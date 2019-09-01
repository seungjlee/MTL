//
// Math Template Library
//
// Copyright (c) 2015: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/sparse.h>

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
    Ax_.Resize(Ai.Size());
  }

  MTL_INLINE void Set(const DynamicVector<I32>& Ap, const DynamicVector<I32>& Ai,
                      const DynamicVector<T>& Ax)
  {
    Ap_ = Ap;
    Ai_ = Ai;
    Ax_ = Ax;
  }

  MTL_INLINE I32 Rows() const             { return Rows_;       }
  MTL_INLINE I32 Cols() const             { return Cols_;       }

  MTL_INLINE const I32* Ap() const        { return Ap_.Begin(); }
  MTL_INLINE const I32* Ai() const        { return Ai_.Begin(); }
  MTL_INLINE const T* Ax() const          { return Ax_.Begin(); }

  MTL_INLINE I32* Ap()                    { return Ap_.Begin(); }
  MTL_INLINE I32* Ai()                    { return Ai_.Begin(); }
  MTL_INLINE T* Ax()                      { return Ax_.Begin(); }

  MTL_INLINE SizeType NumberOfElements()  { return Ai_.Size();  }

protected:
  I32 Rows_;
  I32 Cols_;
  DynamicVector<I32> Ap_;
  DynamicVector<I32> Ai_;
  DynamicVector<T>   Ax_;

  void Release()
  {
    Ap_.Release();
    Ai_.Release();
    Ax_.Release();
  }
};

// Compressed sparse column matrix.
template <class T = F64>
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

  // Constructor that preallocates the buffers.
  CompressedSparseMatrix(I32 rows, I32 cols, SizeType allocationSize)
    : SparseMatrix<T>(rows, cols)
  {
    SparseMatrix<T>::Ap_.Resize(this->Cols_ + 1);
    Allocate(allocationSize);
  }

  // Assumes all row indices in each column are ordered.
  CompressedSparseMatrix(I32 rows, I32 cols,
                         const DynamicVector<DynamicVector<I32>>& sparseColumns)
  {
    Create(rows, cols, sparseColumns);
  }
  CompressedSparseMatrix(I32 rows, I32 cols,
                         const DynamicVector<DynamicVector<I32>>& sparseColumns,
                         const DynamicVector<DynamicVector<T>>& sparseValues)
  {
    Create(rows, cols, sparseColumns, sparseValues);
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

  void Allocate(SizeType allocationSize)
  {
    SparseMatrix<T>::Ai_.Resize(allocationSize);
    SparseMatrix<T>::Ax_.Resize(allocationSize);
  }

  // Expects row numbers in each column to be sorted.
  void Create(I32 rows, I32 cols,
              const DynamicVector<DynamicVector<I32>>& sparseColumns,
              bool updateOptimizedMultiplyStructure = false)
  {
    SizeType totalDataSize = 0;
    for (const auto& sparseCol : sparseColumns)
      totalDataSize += sparseCol.Size();

    ReleaseCached();
    SparseMatrix<T>::Release();
    SparseMatrix<T>::Ap_.Reserve(cols+1);
    SparseMatrix<T>::Ai_.Reserve(totalDataSize);

    SparseMatrix<T>::Rows_ = rows;
    SparseMatrix<T>::Cols_ = cols;

    SparseMatrix<T>::Ap_.PushBack(0);
    for (I32 col = 0; col < SparseMatrix<T>::Cols_; col++)
    {
      SparseMatrix<T>::Ai_.AddBack(sparseColumns[col]);
      SparseMatrix<T>::Ap_.PushBack
        (I32(SparseMatrix<T>::Ap_[SparseMatrix<T>::Ap_.Size()-1] + sparseColumns[col].Size()));
    }
    SparseMatrix<T>::Ax_.Resize(SparseMatrix<T>::Ai_.Size());

    if (updateOptimizedMultiplyStructure)
      OptimizeMultiplyTransposeByThis();
  }

  // Expects row numbers in each column to be sorted.
  void Create(I32 rows, I32 cols,
              const DynamicVector<DynamicVector<I32>>& sparseColumns,
              const DynamicVector<DynamicVector<T>>& sparseValues,
              bool updateOptimizedMultiplyStructure = false)
  {
    SizeType totalDataSize = 0;
    for (const auto& sparseCol : sparseColumns)
      totalDataSize += sparseCol.Size();

    ReleaseCached();
    SparseMatrix<T>::Release();
    SparseMatrix<T>::Ap_.Reserve(cols+1);
    SparseMatrix<T>::Ai_.Reserve(totalDataSize);
    SparseMatrix<T>::Ax_.Reserve(totalDataSize);

    SparseMatrix<T>::Rows_ = rows;
    SparseMatrix<T>::Cols_ = cols;

    SparseMatrix<T>::Ap_.PushBack(0);
    for (I32 col = 0; col < SparseMatrix<T>::Cols_; col++)
    {
      assert(sparseColumns[col].Size() == sparseValues[col].Size());

      SparseMatrix<T>::Ai_.AddBack(sparseColumns[col]);
      SparseMatrix<T>::Ax_.AddBack(sparseValues[col]);
      SparseMatrix<T>::Ap_.PushBack
        (I32(SparseMatrix<T>::Ap_[SparseMatrix<T>::Ap_.Size()-1] + sparseColumns[col].Size()));
    }

    if (updateOptimizedMultiplyStructure)
      OptimizeMultiplyTransposeByThis();
  }

  T operator()(I32 row, I32 col) const
  {
    assert(col <= SparseMatrix<T>::Cols_);

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    for (I32 k = ap[col]; k < ap[col+1]; k++)
    {
      if (ai[k] == row)
        return ax[k];
    }

    return T(0);
  }

  T MaxDiagonal() const
  {
    const I32 N = SparseMatrix<T>::Cols_;

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    T maximum = -kINF;

    for (I32 col = 0 ; col < N; col++)
    {
      bool found = false;
      for (I32 k = ap[col]; k < ap[col+1]; k++)
      {
        if (ai[k] == col)
        {
          maximum = Max(maximum, ax[k]);
          found = true;
          break;
        }
      }

      if (found == false)
        maximum = Max(maximum, T(0));
    }

    return maximum;
  }

  void AddToDiagonalsIfEntryExists(T val)
  {
    const I32 N = SparseMatrix<T>::Cols_;

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    T* ax = SparseMatrix<T>::Ax();

    for (I32 col = 0 ; col < N; col++)
    {
      for (I32 k = ap[col]; k < ap[col+1]; k++)
      {
        if (ai[k] == col)
        {
          ax[k] += val;
          break;
        }

        if (ai[k] > col)
        {
          break;
        }
      }
    }
  }

  // x = At * b where At is the transpose of this matrix.
  void MultiplyTransposed(DynamicVector<T>& x, const DynamicVector<T>& b) const
  {
    assert(SparseMatrix<T>::Rows_ == b.Size());

    const I32 N = SparseMatrix<T>::Cols_;

    x.Resize(N);
    x.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    for (I32 col = 0 ; col < N; col++)
    {
      for (I32 k = ap[col]; k < ap[col+1]; k++)
      {
        x[col] += b[ai[k]] * ax[k];
      }
    }
  }

  // Returns P = At * A where A is this matrix.
  void MultiplyTransposeByThis(CompressedSparseMatrix<T>& P) const
  {
    DynamicVector<DynamicVector<I32>>&
      sparseRows = const_cast<CompressedSparseMatrix*>(this)->SparseRows_;
    DynamicVector<DynamicVector<T>>&
      sparseValues = const_cast<CompressedSparseMatrix*>(this)->SparseValues_;

    const I32 N = SparseMatrix<T>::Cols_;

    sparseRows.Resize(N);
    sparseValues.Resize(N);
    for (I32 i = 0; i < N; i++)
    {
      sparseRows[i].Clear();
      sparseValues[i].Clear();
    }

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

    if (RowIndexPairs_.Size() > 0)
    {
      for (I32 index0 = 0; index0 < (I32)ColumnPairs_.Size(); index0++)
      {
        I32 pairIndex = 0;
        for (I32 index1 = 0; index1 < (I32)ColumnPairs_[index0].Size(); index1++)
        {
          I32 i = ColumnPairs_[index0][index1].x();
          I32 j = ColumnPairs_[index0][index1].y();

          if (i == j)
          {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
            T sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                          numberOfThreads);
#else
            T sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
            sparseRows[i].PushBack(i);
            sparseValues[i].PushBack(sum);
          }
          else
          {
            const Point2D<I32>*
              p = RowIndexPairs_[index0].Begin() + RowIndexPairsStart_[index0][pairIndex];
            const Point2D<I32>*
              pEnd = RowIndexPairs_[index0].Begin() + RowIndexPairsEnd_[index0][pairIndex];

            T sum = 0;

            for (; p < pEnd; p++)
              sum += ax[p->x()] * ax[p->y()];

            sparseRows[i].PushBack(j);
            sparseValues[i].PushBack(sum);

            pairIndex++;
          }
        }
      }
    }
    else
    {
      for (I32 i = 0; i < N; i++)
      {
        for (I32 j = 0; j < i+1; j++)
        {
          if (i == j)
          {
            if (ap[i+1] > ap[i])
            {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
              T sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                            numberOfThreads);
#else
              T sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
              sparseRows[i].PushBack(i);
              sparseValues[i].PushBack(sum);
            }
          }
          else
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

            sparseRows[i].PushBack(j);
            sparseValues[i].PushBack(sum);
          }
        }
      }
    }

    for (I32 col = 1; col < N; col++)
    {
      for (U32 i = 0; i < sparseRows[col].Size(); i++)
      {
        I32 row = sparseRows[col][i];
        if (row != col)
        {
          sparseRows[row].PushBack(col);
          sparseValues[row].PushBack(sparseValues[col][i]);
        }
      }
    }

    P.Create(N, N, sparseRows, sparseValues, false);
  }

  // Returns P = At * A where A is this matrix.
  void MultiplyTransposeByThisParallel(CompressedSparseMatrix<T>& P) const
  {
    DynamicVector<DynamicVector<I32>>&
      sparseRows = const_cast<CompressedSparseMatrix*>(this)->SparseRows_;
    DynamicVector<DynamicVector<T>>&
      sparseValues = const_cast<CompressedSparseMatrix*>(this)->SparseValues_;

    const I32 N = SparseMatrix<T>::Cols_;

    sparseRows.Resize(N);
    sparseValues.Resize(N);
    for (I32 i = 0; i < N; i++)
    {
      sparseRows[i].Clear();
      sparseValues[i].Clear();
    }

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

    if (RowIndexPairs_.Size() > 0)
    {
      MTL_PARALLEL_FOR_BLOCKS(ColumnPairs_.Size())
      for (I32 index0 = 0; index0 < (I32)ColumnPairs_.Size(); index0++)
      {
        I32 pairIndex = 0;
        for (I32 index1 = 0; index1 < (I32)ColumnPairs_[index0].Size(); index1++)
        {
          I32 i = ColumnPairs_[index0][index1].x();
          I32 j = ColumnPairs_[index0][index1].y();

          if (i == j)
          {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
            T sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                          numberOfThreads);
#else
            T sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
            sparseRows[i].PushBack(i);
            sparseValues[i].PushBack(sum);
          }
          else
          {
            const Point2D<I32>*
              p = RowIndexPairs_[index0].Begin() + RowIndexPairsStart_[index0][pairIndex];
            const Point2D<I32>*
              pEnd = RowIndexPairs_[index0].Begin() + RowIndexPairsEnd_[index0][pairIndex];

            T sum = 0;

            for (; p < pEnd; p++)
              sum += ax[p->x()] * ax[p->y()];

            sparseRows[i].PushBack(j);
            sparseValues[i].PushBack(sum);

            pairIndex++;
          }
        }
      }
    }
    else
    {
      MTL_PARALLEL_FOR_BLOCKS(N)
      for (I32 i = 0; i < N; i++)
      {
        for (I32 j = 0; j < i+1; j++)
        {
          if (i == j)
          {
            if (ap[i+1] > ap[i])
            {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
              T sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                            numberOfThreads);
#else
              T sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
              sparseRows[i].PushBack(i);
              sparseValues[i].PushBack(sum);
            }
          }
          else
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

            sparseRows[i].PushBack(j);
            sparseValues[i].PushBack(sum);
          }
        }
      }
    }

    for (I32 col = 1; col < N; col++)
    {
      for (U32 i = 0; i < sparseRows[col].Size(); i++)
      {
        I32 row = sparseRows[col][i];
        if (row != col)
        {
          sparseRows[row].PushBack(col);
          sparseValues[row].PushBack(sparseValues[col][i]);
        }
      }
    }

    P.Create(N, N, sparseRows, sparseValues, false);
  }

  // Returns P = At * A as a dense matrix where A is this matrix.
  void MultiplyTransposeByThis(DynamicMatrix<T>& P) const
  {
    P.Resize(SparseMatrix<T>::Cols_, SparseMatrix<T>::Cols_);
    P.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

    if (RowIndexPairs_.Size() > 0)
    {
      for (I32 index0 = 0; index0 < (I32)ColumnPairs_.Size(); index0++)
      {
        I32 pairIndex = 0;
        for (I32 index1 = 0; index1 < (I32)ColumnPairs_[index0].Size(); index1++)
        {
          I32 i = ColumnPairs_[index0][index1].x();
          I32 j = ColumnPairs_[index0][index1].y();

          if (i == j)
          {
            T sum = 0;
            if (ap[i+1] > ap[i])
            {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
              sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                          numberOfThreads);
#else
              sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
            }
            P[i][i] = sum;
          }
          else
          {
            const Point2D<I32>*
              p = RowIndexPairs_[index0].Begin() + RowIndexPairsStart_[index0][pairIndex];
            const Point2D<I32>*
              pEnd = RowIndexPairs_[index0].Begin() + RowIndexPairsEnd_[index0][pairIndex];

            T sum = 0;

            for (; p < pEnd; p++)
              sum += ax[p->x()] * ax[p->y()];

            P[i][j] = sum;
            P[j][i] = sum;

            pairIndex++;
          }
        }
      }
    }
    else
    {
      for (I32 i = 0; i < SparseMatrix<T>::Cols_; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                        numberOfThreads);
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

  // Returns P = At * A as a dense matrix where A is this matrix.
  void MultiplyTransposeByThisParallel(DynamicMatrix<T>& P) const
  {
    const I32 N = SparseMatrix<T>::Cols_;

    P.Resize(N, N);
    P.Zeros();

    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();
    const T* ax = SparseMatrix<T>::Ax();

    I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

    if (RowIndexPairs_.Size() > 0)
    {
      MTL_PARALLEL_FOR_BLOCKS(ColumnPairs_.Size())
      for (I32 index0 = 0; index0 < (I32)ColumnPairs_.Size(); index0++)
      {
        I32 pairIndex = 0;
        for (I32 index1 = 0; index1 < (I32)ColumnPairs_[index0].Size(); index1++)
        {
          I32 i = ColumnPairs_[index0][index1].x();
          I32 j = ColumnPairs_[index0][index1].y();

          if (i == j)
          {
            T sum = 0;
            if (ap[i+1] > ap[i])
            {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
              sum = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                          numberOfThreads);
#else
              sum = SumOfSquares_Sequential(ax + ap[i], ax + ap[i+1]);
#endif
            }
            P[i][i] = sum;
          }
          else
          {
            const Point2D<I32>*
              p = RowIndexPairs_[index0].Begin() + RowIndexPairsStart_[index0][pairIndex];
            const Point2D<I32>*
              pEnd = RowIndexPairs_[index0].Begin() + RowIndexPairsEnd_[index0][pairIndex];

            T sum = 0;

            for (; p < pEnd; p++)
              sum += ax[p->x()] * ax[p->y()];

            P[i][j] = sum;
            P[j][i] = sum;

            pairIndex++;
          }
        }
      }
    }
    else
    {
      MTL_PARALLEL_FOR_BLOCKS(N)
      for (I32 i = 0; i < N; i++)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i][i] = SumOfSquares_StreamUnaligned_Parallel(ax + ap[i], ap[i+1] - ap[i],
                                                        numberOfThreads);
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
    I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfCores();

    // Mostly logical operations so we can take advantage of hyperthreading.
    if (MTL::CPU::Instance().IsIntel() && MTL::CPU::Instance().Multithreading())
      numberOfThreads *= 2;

    DynamicVector<DynamicVector<Point2D<I32>>> ThreadRowIndexPairs(numberOfThreads);
    DynamicVector<DynamicVector<I32>> ThreadNonZeroMap(numberOfThreads);

    const I32 N = SparseMatrix<T>::Cols_;

    RowIndexPairs_.Resize(N);
    RowIndexPairsStart_.Resize(N);
    RowIndexPairsEnd_.Resize(N);
    ColumnPairs_.Resize(N);

    for (I32 i = 0; i < N; i++)
    {
      RowIndexPairs_[i].Clear();
      ColumnPairs_[i].Clear();
      RowIndexPairsStart_[i].Clear();
      RowIndexPairsEnd_[i].Clear();
    }
    ColumnPairs_[0].PushBack(Point2D<I32>(0,0));

    //
    // Attempting to balance the load for the threads.
    //
    MTL_PARALLEL_FOR_BLOCKS_THREADS(N/2, numberOfThreads)
    for (I32 i = 1; i < N/2+1; i++)
    {
      I32 threadNumber = omp_get_thread_num();

      ThreadRowIndexPairs[threadNumber].Reserve(SparseMatrix<T>::Rows_);

      I32 j = N - i;

      OptimizeMultiplyTransposeByThis(i, ThreadRowIndexPairs[threadNumber],
                                      ThreadNonZeroMap[threadNumber]);
      if (j > i)
        OptimizeMultiplyTransposeByThis(j, ThreadRowIndexPairs[threadNumber],
                                        ThreadNonZeroMap[threadNumber]);
    }
  }

protected:
  // For computing MultiplyTransposeByThis().
  DynamicVector<DynamicVector<I32>> SparseRows_;
  DynamicVector<DynamicVector<T>> SparseValues_;

  // Data cached to optimize computation of At * A
  DynamicVector<DynamicVector<Point2D<I32>>> ColumnPairs_;
  DynamicVector<DynamicVector<Point2D<I32>>> RowIndexPairs_;
  DynamicVector<DynamicVector<I32>> RowIndexPairsStart_;
  DynamicVector<DynamicVector<I32>> RowIndexPairsEnd_;

  void ReleaseCached()
  {
    SparseRows_.Release();
    SparseValues_.Release();
    ColumnPairs_.Release();
    RowIndexPairs_.Release();
    RowIndexPairsStart_.Release();
    RowIndexPairsEnd_.Release();
  }

private:
  void OptimizeMultiplyTransposeByThis(I32 i, DynamicVector<Point2D<I32>>& rowIndexPairs,
                                       DynamicVector<I32>& nonZeroMap)
  {
    const I32* ap = SparseMatrix<T>::Ap();
    const I32* ai = SparseMatrix<T>::Ai();

    I32 p = ap[i];
    const I32 pEnd = ap[i+1];
    const I32 pRows = pEnd - p;

    if(pRows > 0)
    {
      const I32 rows = SparseMatrix<T>::Rows_;

      const I32 pFirstNonZeroRow = ai[p];
      const I32 pLastNonZeroRow = ai[pEnd - 1];

      I32* pNonZeroMap = NULL;
      if (pRows != rows)
      {
        assert(pLastNonZeroRow > pFirstNonZeroRow);
        nonZeroMap.Resize(pLastNonZeroRow - pFirstNonZeroRow + 1);
        nonZeroMap.Zeros();
        pNonZeroMap = nonZeroMap.Begin();

        for (I32 k = ap[i]; k < ap[i+1]; k++)
          pNonZeroMap[ai[k] - pFirstNonZeroRow] = k+1;
      }

      ColumnPairs_[i].Clear();
      RowIndexPairs_[i].Clear();

      for (I32 j = 0; j < i; j++)
      {
        I32 count = 0;
        I32 q = ap[j];
        const I32 qEnd = ap[j+1];
        const I32 qRows = qEnd - q;

        p = ap[i];

        if (pRows != rows && qRows != rows)
        {
          if (ai[p] <= ai[qEnd - 1] && ai[q] <= pLastNonZeroRow)
          {
            rowIndexPairs.Clear();

            while (q < qEnd && ai[q] <= ai[pEnd - 1])
            {
              if (ai[q] < ai[p])
              {
                if (ai[p] > ai[qEnd - 1])
                  break;

                while (q < qEnd && ai[q] < ai[p])
                  q++;
              }

              if (q >= qEnd)
                break;

              if (ai[q] == ai[p])
              {
                rowIndexPairs.PushBack(Point2D<I32>(p,q));

                if (ai[p] == ai[qEnd - 1])
                  break;
                if (++p == pEnd)
                  break;
              }
              else
              {
                assert(ai[q] >= pFirstNonZeroRow);

                if (ai[q] > pLastNonZeroRow)
                  break;

                U32 pp = pNonZeroMap[ai[q] - pFirstNonZeroRow];
                if (pp)
                {
                  rowIndexPairs.PushBack(Point2D<I32>(pp-1,q));

                  if (pp == pEnd || ai[pp-1] == ai[qEnd - 1])
                    break;

                  p = pp;
                }
              }
              q++;
            }

            if (rowIndexPairs.Size() > 0)
            {
              RowIndexPairsStart_[i].PushBack((I32)RowIndexPairs_[i].Size());
              RowIndexPairs_[i].AddBack(rowIndexPairs);
              RowIndexPairsEnd_[i].PushBack((I32)RowIndexPairs_[i].Size());
              ColumnPairs_[i].PushBack(Point2D<I32>(i,j));
            }
          }
        }
        else if (pRows != rows)
        {
          SizeType addSize = pEnd - p;
          I32 startIndex = (I32)RowIndexPairs_[i].Size();
          RowIndexPairsStart_[i].PushBack(startIndex);
          RowIndexPairs_[i].Resize(startIndex + addSize);

          Point2D<I32>* pPairs = &RowIndexPairs_[i][startIndex];
          for (I32 k = 0; k < addSize; k++)
          {
            pPairs[k].x(p);
            pPairs[k].y(q + ai[p]);
            p++;
          }

          RowIndexPairsEnd_[i].PushBack((I32)RowIndexPairs_[i].Size());
          ColumnPairs_[i].PushBack(Point2D<I32>(i,j));
        }
        else if (qRows != rows)
        {
          SizeType addSize = qEnd - q;
          I32 startIndex = (I32)RowIndexPairs_[i].Size();
          RowIndexPairsStart_[i].PushBack(startIndex);
          RowIndexPairs_[i].Resize(startIndex + addSize);

          Point2D<I32>* pPairs = &RowIndexPairs_[i][startIndex];
          for (I32 k = 0; k < addSize; k++)
          {
            pPairs[k].x(p + ai[q]);
            pPairs[k].y(q);
            q++;
          }

          RowIndexPairsEnd_[i].PushBack((I32)RowIndexPairs_[i].Size());
          ColumnPairs_[i].PushBack(Point2D<I32>(i,j));
        }
        else
        {
          I32 startIndex = (I32)RowIndexPairs_[i].Size();
          RowIndexPairsStart_[i].PushBack(startIndex);
          RowIndexPairs_[i].Resize(startIndex + rows);

          Point2D<I32>* pPairs = &RowIndexPairs_[i][startIndex];
          for (I32 k = 0; k < rows; k++)
          {
            pPairs[k].x(p++);
            pPairs[k].y(q++);
          }

          RowIndexPairsEnd_[i].PushBack((I32)RowIndexPairs_[i].Size());
          ColumnPairs_[i].PushBack(Point2D<I32>(i,j));
        }
      }

      ColumnPairs_[i].PushBack(Point2D<I32>(i,i));
    }
  }

  template<class TT>
  SizeType ComputeTotalSize(const DynamicVector<DynamicVector<DynamicVector<TT>>>& v)
  {
    SizeType total = 0;

    for (U32 i = 0; i < v.Size(); i++)
    {
      for (U32 k = 0; k < v[i].Size(); k++)
      {
        total += v[i][k].Size();
      }
    }

    return total;
  }
};

template <class T>
static I32 SolveLDLt(DynamicVector<T>& x, const CompressedSparseMatrix<T>& A,
                     const DynamicVector<T>& b, SymbolicLDLt<T>& symLDLt,
                     const T& tolerance = NumericalEpsilon<T>())
{
  assert(A.Cols() == A.Rows());

  I32 n = A.Cols();

  const I32* Ap = A.Ap();
  const I32* Ai = A.Ai();
  const T* Ax = A.Ax();

  DynamicVector<T> y(n);
  DynamicVector<T> d(n);
  T* Y = y.Begin();
  T* D = d.Begin();

  const I32* P = symLDLt.P().Begin();

  const I32* Lp = symLDLt.L().Ap();
  I32* Li = symLDLt.L().Ai();
  T*   Lx = symLDLt.L().Ax();

  LDLt_Numeric(n, Ap, Ai, Ax, Lp, symLDLt.Parent(), symLDLt.Lnz(), Li, Lx, D, Y,
               symLDLt.Pattern(), symLDLt.Flag(), P, symLDLt.Pinv());

  I32 rank = 0;
  for (I32 i = 0; i < n; i++)
    if (Abs(D[i]) > tolerance)
      rank++;

  for (int i = 0; i < n; i++)
    Y[i] = b[P[i]];

  // Compute inv(L) * y using forward substitution.
  for (int j = 0; j < n; j++)
  {
    for (int p = Lp[j]; p < Lp[j+1]; p++)
      Y[Li[p]] -= Lx[p] * Y[j];
  }

  // Compute inv(D) * y.
  y /= d;

  // Finally, compute inv(Lt) * y using back substitution.
  for (int j = n-1; j >= 0; j--)
  {
    for (int p = Lp[j]; p < Lp[j+1]; p++)
      Y[j] -= Lx[p] * Y[Li[p]];
  }

  x.Resize(n);
  T* X = x.Begin();
  for (int i = 0; i < n; i++)
    X[P[i]] = Y[i];

  return rank;
}



}  // namespace MTL

#endif // MTL_SPARSE_MATRIX_H
