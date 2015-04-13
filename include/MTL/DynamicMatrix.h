//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_DYNAMIC_MATRIX_H
#define MTL_DYNAMIC_MATRIX_H

#include "DynamicVector.h"

namespace MTL
{

// P = M * Mt.
template<class T>
MTL_INLINE static void MultiplyByTranspose(T* P, const T* M, I32 rows, I32 cols,
                                           I32 rowSizeP, I32 rowSizeM)
{
  for (I32 i = 0; i < rows; i++)
  {
    for (I32 j = i; j < rows; j++)
    {
      if (i == j)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i*rowSizeP + j] = SumOfSquares_StreamAligned_Parallel(M + i*rowSizeM, cols);
#else
        P[i*rowSizeP + j] = SumOfSquares_Sequential(M + i*rowSizeM, M + i*rowSizeM + cols);
#endif
      }
      else
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i*rowSizeP + j] = DotProduct_StreamAligned_Parallel(M + i*rowSizeM,
                                                              M + j*rowSizeM, cols);
#else
        P[i*rowSizeP + j] = DotProduct_Sequential(M + i*rowSizeM,
                                                  M + j*rowSizeM, M + i*rowSizeM + cols);
#endif
        P[j*rowSizeP + i] = P[i*rowSizeP + j];
      }
    }
  }
}

template<class T>
MTL_INLINE static void Multiply(T* P, const T* A, const T* B, I32 rows, I32 N, I32 cols,
                                I32 rowSizeP, I32 rowSizeA, I32 rowSizeB)
{
  for (I32 row = 0; row < rows; row++)
  {
    for (I32 col = 0; col < cols; col++)
    {
      P[row*rowSizeP + col] = A[row*rowSizeA] * B[col];
      for (I32 i = 1; i < N; i++)
        P[row*rowSizeP + col] += A[row*rowSizeA + i] * B[i*rowSizeB + col];
    }
  }
}

template<class T>
MTL_INLINE static void SwapRows(T* a, I32 row1, I32 row2, I32 dataSize, I32 rowSize)
{
  for (I32 i = 0; i < dataSize; i++)
    Swap(a[row1 * rowSize + i], a[row2 * rowSize + i]);
}

template<class T>
class DynamicMatrix
{
public:
  MTL_INLINE DynamicMatrix() : Rows_(0), Cols_(0) {}

  enum Initialize
  {
    eNothing,
    eIdentity,
    eZeros
  };

  MTL_INLINE DynamicMatrix(I32 rows, I32 cols, I32 byteAlignment = MTL_STREAM_BYTES)
  {
    Resize(rows, cols, byteAlignment);
  }

  MTL_INLINE DynamicMatrix(I32 rows, I32 cols, Initialize i, I32 byteAlignment = MTL_STREAM_BYTES)
  {
    Resize(rows, cols, byteAlignment);

    switch (i)
    {
      case eIdentity:  Identity();    break;
      case eZeros:     Zeros();       break;
      default:                        break;
    }
  }

  template<I32 M, I32 N>
  MTL_INLINE DynamicMatrix(const Matrix<M,N,T>& matrix, I32 byteAlignment = MTL_STREAM_BYTES)
  {
    Resize(M, N);
    for (I32 row = 0; row < M; row++)
      memcpy((*this)[row], matrix[row], N*sizeof(T));
  }

  MTL_INLINE void Resize(I32 rows, I32 cols, I32 byteAlignment = MTL_STREAM_BYTES)
  {
    Rows_ = rows;
    Cols_ = cols;

    RowSize_ =
      ((Cols_ * sizeof(T) + byteAlignment - 1) / byteAlignment) * byteAlignment / sizeof(T);

    Data_.Resize(Rows_ * RowSize_);
  }

  MTL_INLINE void Identity()
  {
    Zeros();
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] = T(1);
  }

  MTL_INLINE void Zeros()
  {
    for (I32 i = 0; i < Rows(); i++)
      memset((*this)[i], 0, Cols() * sizeof(T));
  }

  MTL_INLINE void AddToDiagonals(const T& value)
  {
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] += value;
  }

  MTL_INLINE void SetAll(const T& newVal)
  {
    for (I32 i = 0; i < Rows(); i++)
      OptimizedAssignAll((*this)[i], newVal, Cols());
  }

  MTL_INLINE DynamicMatrix operator*(const DynamicMatrix& B) const
  {
    assert(Cols() == B.Rows());

    DynamicMatrix result(Rows(), B.Cols());
    Multiply(result[0], (*this)[0], B[0], Rows(), Cols(), B.Cols(),
             result.RowSize(), RowSize(), B.RowSize());

    return result;
  }

  MTL_INLINE DynamicVector<T> operator*(const DynamicVector<T>& v) const
  {
    assert(Cols() == (I32)v.Size());

    DynamicVector<T> result(Rows());
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      result[i] = DotProduct_StreamAligned_Sequential((*this)[i], v.Begin(), Cols());
#else
      result[i] = DotProduct_Sequential((*this)[i], v.Begin(), (*this)[i] + Cols());
#endif

    return result;
  }

  // Computes (*this) * this->getTranspose().
  MTL_INLINE DynamicMatrix MultiplyByTranspose() const
  {
    DynamicMatrix squareSymmetricMatrix(Rows(), Rows());
    ::MultiplyByTranspose(squareSymmetricMatrix[0], (*this)[0], Rows(), Cols(),
                          squareSymmetricMatrix.RowSize(), RowSize());

    return squareSymmetricMatrix;
  }

  MTL_INLINE DynamicMatrix operator+=(const DynamicMatrix& B)
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      Addition_StreamAligned_Sequential((*this)[i], B[i], Cols());
#else
      Addition_Sequential((*this)[i], B[i], (*this)[i] + Cols());
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator+(const DynamicMatrix& B) const
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    Matrix C = *this;
    C += B;
    return C;
  }

  MTL_INLINE DynamicMatrix& operator-=(const DynamicMatrix& B)
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      Subtraction_StreamAligned_Sequential((*this)[i], B[i], Cols());
#else
      Subtraction_Sequential((*this)[i], B[i], (*this)[i] + Cols());
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator-(const DynamicMatrix& B) const
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    DynamicMatrix C = *this;
    C -= B;
    return C;
  }

  MTL_INLINE DynamicMatrix operator-() const
  {
    DynamicMatrix A = *this;
    for (I32 i = 0; i < A.Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      UnaryMinus_StreamAligned_Sequential(A[i], A.Cols());
#else
      UnaryMinus_Sequential(A[i], A[i] + A.Cols());
#endif

    return A;
  }

  MTL_INLINE DynamicMatrix& operator+=(const T& scalar)
  {
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarAddition_StreamAligned_Sequential(A[i], scalar, A.Cols())
#else
      ScalarAddition_Sequential(A[i], scalar, A[i] + A.Cols())
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator+(const T& scalar) const
  {
    DynamicMatrix A = *this;
    A += scalar;
    return A;
  }

  MTL_INLINE DynamicMatrix& operator-=(const T& scalar)
  {
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarSubtraction_StreamAligned_Sequential(A[i], scalar, A.Cols())
#else
      ScalarSubtraction_Sequential(A[i], scalar, A[i] + A.Cols())
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator-(const T& scalar) const
  {
    DynamicMatrix A = *this;
    A -= scalar;
    return A;
  }

  MTL_INLINE DynamicMatrix& operator*=(const T& scalar)
  {
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarMultiplication_StreamAligned_Sequential(A[i], scalar, A.Cols())
#else
      ScalarMultiplication_Sequential(A[i], scalar, A[i] + A.Cols())
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator*(const T& scalar) const
  {
    DynamicMatrix A = *this;
    A *= scalar;
    return A;
  }

  MTL_INLINE DynamicMatrix& operator/=(const T& scalar)
  {
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarDivision_StreamAligned_Sequential(A[i], scalar, A.Cols())
#else
      ScalarDivision_Sequential(A[i], scalar, A[i] + A.Cols())
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator/(const T& scalar) const
  {
    DynamicMatrix A = *this;
    A /= scalar;
    return A;
  }

  MTL_INLINE bool IsZero() const  { return *this == DynamicMatrix(eZeros); }

  MTL_INLINE T Sum() const
  {
    T sum = T(0);
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      sum += Sum_StreamAligned_Parallel((*this)[i], Cols());
#else
      sum += Sum_Sequential((*this)[i], (*this)[i] + Cols());
#endif

    return sum;
  }

  MTL_INLINE T SumOfSquares() const
  {
    T sum = T(0);
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      sum += SumOfSquares_StreamAligned_Parallel((*this)[i], Cols());
#else
      sum += SumOfSquares_Sequential((*this)[i], (*this)[i] + Cols());
#endif

    return sum;
  }

  MTL_INLINE double Mean() const
  {
    return Sum() / double(Rows()*Cols());
  }

  MTL_INLINE T MaxNorm() const
  {
    T max = T(0);
    for (I32 i = 0; i < Rows(); i++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      T sum = SumOfAbsolutes_StreamAligned_Parallel((*this)[i], Cols());
#else
      T sum = SumOfAbsolutes_Sequential((*this)[i], (*this)[i] + Cols());
#endif
      if (sum > max)
        max = sum;
    }

    return max;
  }

  MTL_INLINE T FrobeniusNorm() const
  {
    return Sqrt(SumOfSquares());
  }

  MTL_INLINE DynamicMatrix ComputeTranspose() const
  {
    DynamicMatrix t(Cols(),Rows());
    for (I32 i = 0; i < Rows(); i++)
      for (I32 j = 0; j < Cols(); j++)
        t[j][i] = (*this)[i][j];

    return t;
  }

  void SetRow(I32 row, const DynamicVector<T>& v)
  {
    assert(row >= 0 && row < Rows());
    assert(v.Size() <= Cols());

    OptimizedCopy((*this)[row], v.Begin(), v.Size());
  }

  MTL_INLINE const T* operator[](I32 row) const
  {
    assert(row >= 0 && row < Rows());
    return Data_.Begin() + row * RowSize();
  }
  MTL_INLINE T* operator[](I32 row)
  {
    assert(row >= 0 && row < Rows());
    return Data_.Begin() + row * RowSize();
  }

  MTL_INLINE I32 Rows() const              { return Rows_;         }
  MTL_INLINE I32 Cols() const              { return Cols_;         }
  MTL_INLINE I32 RowSize() const           { return RowSize_;      }

  MTL_INLINE const T* Data() const         { return Data_.begin(); }
  MTL_INLINE T* Data()                     { return Data_.begin(); }

private:
  I32 Rows_;
  I32 Cols_;
  I32 RowSize_;
  DynamicVector<T> Data_;
};

}  // namespace MTL

#endif // MTL_DYNAMIC_MATRIX_H
