//
// Math Template Library
//
// Copyright (c) 2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice, this list of
//      conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright notice, this list of
//      conditions and the following disclaimer in the documentation and/or other materials
//      provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
// WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef MTL_SPARSE_MATRIX_H
#define MTL_SPARSE_MATRIX_H

#include <MTL/DynamicMatrix.h>
#include <MTL/Point2D.h>

namespace MTL
{

// Base class for sparse matrices in various formats such as a coordinate list or
// compressed sparse columns.
template <class T>
class SparseMatrix
{
public:
  SparseMatrix() {}
  SparseMatrix(I32 rows, I32 cols) : Rows_(rows), Cols_(cols) {}
  SparseMatrix(I32 rows, I32 cols,
               const DynamicVector<I32>& Ap,
               const DynamicVector<I32>& Ai,
               const DynamicVector<T>& Ax)
    : Ap_(Ap), Ai_(Ai), Ax_(Ax)
  {
  }

  MTL_INLINE void Set(const DynamicVector<I32>& Ap, const DynamicVector<I32>& Ai)
  {
    Ap_ = Ap;
    Ai_ = Ai;
    Ax_.Resize(Ai);
  }

  MTL_INLINE void Set(const DynamicVector<I32>& Ap, const DynamicVector<I32>& Ai,
                      const DynamicVector<T>& Ax)
  {
    Ap_ = Ap;
    Ai_ = Ai;
    Ax_ = Ax;
  }

  MTL_INLINE I32 Rows() const       { return Rows_;       }
  MTL_INLINE I32 Cols() const       { return Cols_;       }

  MTL_INLINE const I32* Ap() const  { return Ap_.Begin(); }
  MTL_INLINE const I32* Ai() const  { return Ai_.Begin(); }
  MTL_INLINE const T* Ax() const    { return Ax_.Begin(); }

  MTL_INLINE T* Ax()  { return Ax_.Begin(); }

protected:
  I32 Rows_;
  I32 Cols_;
  DynamicVector<I32> Ap_;
  DynamicVector<I32> Ai_;
  DynamicVector<T>   Ax_;

  void Clear()
  {
    Ap_.Clear();
    Ai_.Clear();
    Ax_.Clear();
  }

  // Data cached to optimize computation of At * A
  DynamicVector<DynamicVector<DynamicVector<Point2D<I32>>>> MultiplyTransposeIndices_;
};

// Compressed sparse column matrix.
template <class T>
class CompressedSparseMatrix : public SparseMatrix<T>
{
public:
  CompressedSparseMatrix() {}
  CompressedSparseMatrix(I32 rows, I32 cols) : SparseMatrix<T>(rows, cols) {}
  CompressedSparseMatrix(I32 rows, I32 cols,
                         const DynamicVector<I32>& Ap,
                         const DynamicVector<I32>& Ai,
                         const DynamicVector<T>& Ax)
    : SparseMatrix<T>(rows, cols, Ap, Ai, Ax)
  {
    OptimizeMultiplyTransposeByThis();
  }

  // Assumes all row indices in each column are ordered.
  CompressedSparseMatrix(I32 rows, I32 cols,
                         const DynamicVector<DynamicVector<I32>>& sparseColumns)
  {
    Create(rows, cols, sparseColumns);
  }

  CompressedSparseMatrix(const DynamicMatrix<T>& full)
  {
    Rows_ = full.Rows();
    Cols_ = full.Cols();

    for (I32 col = 0; col < Cols_; col++)
    {
      I32 lastSize = (I32)Ai_.Size();

      for (I32 row = 0; row < Rows_; row++)
      {
        if (full[row][col] != 0)
        {
          Ai_.PushBack(row);
          Ax_.PushBack(full[row][col]);
        }
      }

      Ap_.PushBack(lastSize);
    }

    Ap_.PushBack((I32)Ai_.Size());
  }

  void Create(I32 rows, I32 cols,
              const DynamicVector<DynamicVector<I32>>& sparseColumns)
  {
    Rows_ = rows;
    Cols_ = cols;
    Clear();

    Ap_.PushBack(0);
    for (I32 col = 0; col < Cols_; col++)
    {
      Ai_.AddBack(sparseColumns[col]);
      Ap_.PushBack(I32(Ap_[Ap_.Size()-1] + sparseColumns[col].Size()));
    }
    Ax_.Resize(Ai_.Size());

    OptimizeMultiplyTransposeByThis();
  }

  // x = At * b where At is the transpose of this matrix.
  void MultiplyTransposed(DynamicVector<T>& x, const DynamicVector<T>& b) const
  {
    assert(Rows_ == b.Size());

    x.Resize(Cols_);
    x.Zeros();

    const I32* ap = Ap();
    const I32* ai = Ai();
    const T* ax = Ax();

    for (I32 col = 0 ; col < Cols_; col++)
    {
      for (I32 k = ap[col]; k < ap[col+1]; k++)
      {
        x[col] += b[ai[k]] * ax[k];
      }
    }
  }

  // Returns P = At * A as a dense matrix where A is this matrix.
  void MultiplyTransposeByThis(DynamicMatrix<T>& P) const
  {
    P.Resize(Cols_, Cols_);
    P.Zeros();

    const I32* ap = Ap();
    const I32* ai = Ai();
    const T* ax = Ax();

    if (MultiplyTransposeIndices_.Size() == Cols_)
    {
      for (I32 i = 0; i < Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < Cols_; j++)
        {
          const DynamicVector<Point2D<I32>>& pairs = MultiplyTransposeIndices_[i][j - i - 1];
          const Point2D<I32>* pPairs    = pairs.Begin();
          const Point2D<I32>* pPairsEnd = pairs.End();

          T sum = 0;

          for (; pPairs < pPairsEnd; pPairs++)
            sum += ax[pPairs->x()] * ax[pPairs->y()];

          P[i][j] = sum;
          P[j][i] = sum;
        }
      }
    }
    else
    {
      for (I32 i = 0; i < Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < Cols_; j++)
        {
          T sum = 0;
          I32 p = ap[i];
          I32 q = ap[j];
          I32 pEnd = ap[i+1];
          I32 qEnd = ap[j+1];
          while (p < pEnd && q < qEnd)
          {
            if (ai[p] < ai[q])
            {
              p++;
            }
            else if (ai[p] > ai[q])
            {
              q++;
            }
            else
            {
              sum += ax[p] * ax[q];
              p++;
              q++;
            }
          }

          P[i][j] = sum;
          P[j][i] = sum;
        }
      }
    }
  }
  void MultiplyTransposeByThisParallel(DynamicMatrix<T>& P) const
  {
    P.Resize(Cols_, Cols_);
    P.Zeros();

    const I32* ap = Ap();
    const I32* ai = Ai();
    const T* ax = Ax();

    if (MultiplyTransposeIndices_.Size() == Cols_)
    {
      #pragma omp parallel for
      for (I32 i = 0; i < Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < Cols_; j++)
        {
          const DynamicVector<Point2D<I32>>& pairs = MultiplyTransposeIndices_[i][j - i - 1];
          const Point2D<I32>* pPairs    = pairs.Begin();
          const Point2D<I32>* pPairsEnd = pairs.End();

          T sum = 0;

          for (; pPairs < pPairsEnd; pPairs++)
            sum += ax[pPairs->x()] * ax[pPairs->y()];

          P[i][j] = sum;
          P[j][i] = sum;
        }
      }
    }
    else
    {
      #pragma omp parallel for
      for (I32 i = 0; i < Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < Cols_; j++)
        {
          T sum = 0;
          I32 p = ap[i];
          I32 q = ap[j];
          I32 pEnd = ap[i+1];
          I32 qEnd = ap[j+1];
          while (p < pEnd && q < qEnd)
          {
            if (ai[p] < ai[q])
            {
              p++;
            }
            else if (ai[p] > ai[q])
            {
              q++;
            }
            else
            {
              sum += ax[p] * ax[q];
              p++;
              q++;
            }
          }

          P[i][j] = sum;
          P[j][i] = sum;
        }
      }
    }
  }

  void OptimizeMultiplyTransposeByThis()
  {
    const I32* ap = Ap();
    const I32* ai = Ai();
    const T* ax = Ax();

    MultiplyTransposeIndices_.Resize(Cols_);

    for (I32 i = 0; i < Cols_ - 1; i++)
    {
      DynamicVector<DynamicVector<Point2D<I32>>>& colIndices = MultiplyTransposeIndices_[i];

      I32 offset = i + 1;
      colIndices.Resize(Cols_ - offset);

      for (I32 j = i + 1; j < Cols_; j++)
      {
        DynamicVector<Point2D<I32>>& pairs = colIndices[j - offset];
        pairs.Clear();

        I32 p = ap[i];
        I32 q = ap[j];
        I32 pEnd = ap[i+1];
        I32 qEnd = ap[j+1];
        while (p < pEnd && q < qEnd)
        {
          if (ai[p] < ai[q])
          {
            p++;
          }
          else if (ai[p] > ai[q])
          {
            q++;
          }
          else
          {
            pairs.PushBack(Point2D<I32>(p,q));
            p++;
            q++;
          }
        }
      }
    }
  }
};


}  // namespace MTL

#endif // MTL_SPARSE_MATRIX_H
