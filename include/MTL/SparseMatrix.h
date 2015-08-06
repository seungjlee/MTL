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
    SparseMatrix<T>::Rows_ = full.Rows();
    SparseMatrix<T>::Cols_ = full.Cols();

    for (I32 col = 0; col < SparseMatrix<T>::Cols_; col++)
    {
      I32 lastSize = (I32)SparseMatrix<T>::Ai_.Size();

      for (I32 row = 0; row < SparseMatrix<T>::Rows_; row++)
      {
        if (full[row][col] != 0)
        {
          SparseMatrix<T>::Ai_.PushBack(row);
          SparseMatrix<T>::Ax_.PushBack(full[row][col]);
        }
      }

      SparseMatrix<T>::Ap_.PushBack(lastSize);
    }

    SparseMatrix<T>::Ap_.PushBack((I32)SparseMatrix<T>::Ai_.Size());
  }

  void Create(I32 rows, I32 cols,
              const DynamicVector<DynamicVector<I32>>& sparseColumns)
  {
    SparseMatrix<T>::Rows_ = rows;
    SparseMatrix<T>::Cols_ = cols;
    SparseMatrix<T>::Clear();

    SparseMatrix<T>::Ap_.PushBack(0);
    for (I32 col = 0; col < SparseMatrix<T>::Cols_; col++)
    {
      SparseMatrix<T>::Ai_.AddBack(sparseColumns[col]);
      SparseMatrix<T>::Ap_.PushBack
        (I32(SparseMatrix<T>::Ap_[SparseMatrix<T>::Ap_.Size()-1] + sparseColumns[col].Size()));
    }
    SparseMatrix<T>::Ax_.Resize(SparseMatrix<T>::Ai_.Size());

    OptimizeMultiplyTransposeByThis();
  }

  // x = At * b where At is the transpose of this matrix.
  void MultiplyTransposed(DynamicVector<T>& x, const DynamicVector<T>& b) const
  {
    assert(SparseMatrix<T>::Rows_ == b.Size());

    x.Resize(SparseMatrix<T>::Cols_);
    x.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    for (I32 col = 0 ; col < SparseMatrix<T>::Cols_; col++)
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
    P.Resize(SparseMatrix<T>::Cols_, SparseMatrix<T>::Cols_);
    P.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    if (Mp_.Size() > 0)
    {
      const I32* Mq = Mq_.Begin();
      const I32* Mp = Mp_.Begin();
      const Point2D<I32>* pPairs = MultiplyTransposeIndices_.Begin();

      for (I32 i = 0; i < SparseMatrix<T>::Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
        if (i < SparseMatrix<T>::Cols_ - 1)
        {
          I32 index = Mq[i];
          for (I32 j = i + 1; j < SparseMatrix<T>::Cols_; j++, index++)
          {
            if (Mp[index] < Mp[index + 1])
            {
              const Point2D<I32>* p = pPairs + Mp[index];
              const Point2D<I32>* pEnd = pPairs + Mp[index + 1];

              T sum = 0;

              for (; p < pEnd; p++)
                sum += ax[p->x()] * ax[p->y()];

              P[i][j] = sum;
              P[j][i] = sum;
            }
          }
        }
      }
    }
    else
    {
      for (I32 i = 0; i < SparseMatrix<T>::Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < SparseMatrix<T>::Cols_; j++)
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
    P.Resize(SparseMatrix<T>::Cols_, SparseMatrix<T>::Cols_);
    P.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    if (Mp_.Size() > 0)
    {
      const I32* Mq = Mq_.Begin();
      const I32* Mp = Mp_.Begin();
      const Point2D<I32>* pPairs = MultiplyTransposeIndices_.Begin();

      #pragma omp parallel for
      for (I32 i = 0; i < SparseMatrix<T>::Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
        if (i < SparseMatrix<T>::Cols_ - 1)
        {
          I32 index = Mq[i];
          for (I32 j = i + 1; j < SparseMatrix<T>::Cols_; j++, index++)
          {
            if (Mp[index] < Mp[index + 1])
            {
              const Point2D<I32>* p = pPairs + Mp[index];
              const Point2D<I32>* pEnd = pPairs + Mp[index + 1];

              T sum = 0;

              for (; p < pEnd; p++)
                sum += ax[p->x()] * ax[p->y()];

              P[i][j] = sum;
              P[j][i] = sum;
            }
          }
        }
      }
    }
    else
    {
      #pragma omp parallel for
      for (I32 i = 0; i < SparseMatrix<T>::Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Sequential(ax + ap[i], ap[i+1] - ap[i]);
#else
        P[i][i] = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif

        for (I32 j = i + 1; j < SparseMatrix<T>::Cols_; j++)
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
    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    Mp_.Clear();
    Mq_.Resize(SparseMatrix<T>::Cols_ - 1);
    MultiplyTransposeIndices_.Clear();

    Mp_.PushBack(0);
    for (I32 i = 0; i < SparseMatrix<T>::Cols_ - 1; i++)
    {
      Mq_[i] = (I32)Mp_.Size() - 1; 
      for (I32 j = i + 1; j < SparseMatrix<T>::Cols_; j++)
      {
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
            MultiplyTransposeIndices_.PushBack(Point2D<I32>(p,q));
            p++;
            q++;
          }
        }
        Mp_.PushBack((I32)MultiplyTransposeIndices_.Size());
      }
    }
  }

protected:
  // Data cached to optimize computation of At * A
  DynamicVector<I32> Mq_;
  DynamicVector<I32> Mp_;
  DynamicVector<Point2D<I32>> MultiplyTransposeIndices_;
};


}  // namespace MTL

#endif // MTL_SPARSE_MATRIX_H
