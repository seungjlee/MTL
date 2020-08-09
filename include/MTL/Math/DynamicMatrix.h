//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT
#define MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT MTL_STREAM_BYTES
#endif

#include <MTL/Math/DynamicVector.h>
#include <MTL/Stream/StreamArray.h>

namespace MTL
{

// P = M * Mt.
template<class T>
MTL_INLINE static void MultiplyByTranspose(T* P, const T* M, I32 rows, I32 cols,
                                           I32 rowSizeP, I32 rowSizeM)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 i = 0; i < rows; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    P[i*rowSizeP + i] = SumOfSquares_StreamAligned_Sequential(M + i*rowSizeM, cols);
#else
    P[i*rowSizeP + i] = SumOfSquares_Sequential(M + i*rowSizeM, cols);
#endif

    for (I32 j = i + 1; j < rows; j++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      P[i*rowSizeP + j] = DotProduct_StreamAligned_Sequential(M + i*rowSizeM,
                                                              M + j*rowSizeM, cols);
#else
      P[i*rowSizeP + j] = DotProduct_Sequential(M + i*rowSizeM,
                                                M + j*rowSizeM, cols);
#endif
      P[j*rowSizeP + i] = P[i*rowSizeP + j];
    }
  }
}
// P += M * Mt.
template<class T>
MTL_INLINE static void AddMultiplyByTranspose(T* P, const T* M, I32 rows, I32 cols,
                                              I32 rowSizeP, I32 rowSizeM)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 i = 0; i < rows; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    P[i*rowSizeP + i] += SumOfSquares_StreamAligned_Sequential(M + i*rowSizeM, cols);
#else
    P[i*rowSizeP + i] += SumOfSquares_Sequential(M + i*rowSizeM, M + i*rowSizeM + cols);
#endif

    for (I32 j = i + 1; j < rows; j++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      P[i*rowSizeP + j] += DotProduct_StreamAligned_Sequential(M + i*rowSizeM,
                                                               M + j*rowSizeM, cols);
#else
      P[i*rowSizeP + j] += DotProduct_Sequential(M + i*rowSizeM,
                                                 M + j*rowSizeM, cols);
#endif
    }
  }

  for (I32 i = 0; i < rows; i++)
    for (I32 j = i; j < rows; j++)
      P[j*rowSizeP + i] = P[i*rowSizeP + j];
}
// Only fills lower part of matrix.
template<class T>
MTL_INLINE static void MultiplyByTranspose_LowerMatrix(T* P, const T* M, I32 rows, I32 cols,
                                                       I32 rowSizeP, I32 rowSizeM)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 i = 0; i < rows; i++)
  {
    for (I32 j = i; j < rows; j++)
    {
      if (i == j)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[i*rowSizeP + j] = SumOfSquares_StreamAligned_Sequential(M + i*rowSizeM, cols);
#else
        P[i*rowSizeP + j] = SumOfSquares_Sequential(M + i*rowSizeM, cols);
#endif
      }
      else
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
        P[j*rowSizeP + i] = DotProduct_StreamAligned_Sequential(M + i*rowSizeM, M + j*rowSizeM, cols);
#else
        P[j*rowSizeP + i] = DotProduct_Sequential(M + i*rowSizeM, M + j*rowSizeM, cols);
#endif
      }
    }
  }
}

template<class T>
MTL_INLINE static void Multiply(T* P, const T* A, const T* B, I32 rows, I32 N, I32 cols,
                                I32 rowSizeP, I32 rowSizeA, I32 rowSizeB)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 row = 0; row < rows; row++)
  {
    for (I32 col = 0; col < cols; col++)
    {
      P[row*rowSizeP + col] =
      P[row*rowSizeP + col] = A[row*rowSizeA] * B[col];
      for (I32 i = 1; i < N; i++)
        P[row*rowSizeP + col] += A[row*rowSizeA + i] * B[i*rowSizeB + col];
    }
  }
}
template<class T>
MTL_INLINE static void MultiplyTransposed(T* P, const T* A, const T* Bt,
                                          I32 rows, I32 N, I32 cols,
                                          I32 rowSizeP, I32 rowSizeA, I32 rowSizeBt)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 row = 0; row < rows; row++)
  {
    for (I32 col = 0; col < cols; col++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      P[row*rowSizeP + col] = DotProduct_StreamAligned_Sequential(A + row*rowSizeA, Bt + col*rowSizeBt, N);
#else
      P[row*rowSizeP + col] = DotProduct_Sequential(A + row*rowSizeA, Bt + col*rowSizeBt, N);
#endif
    }
  }
}
template<class T>
MTL_INLINE static void AddMultiplyTransposed(T* P, const T* A, const T* Bt,
                                            I32 rows, I32 N, I32 cols,
                                            I32 rowSizeP, I32 rowSizeA, I32 rowSizeBt)
{
  MTL_PARALLEL_FOR_BLOCKS(rows)
  for (I32 row = 0; row < rows; row++)
  {
    for (I32 col = 0; col < cols; col++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      P[row*rowSizeP + col] += DotProduct_StreamAligned_Sequential(A + row*rowSizeA, Bt + col*rowSizeBt, N);
#else
      P[row*rowSizeP + col] += DotProduct_Sequential(A + row*rowSizeA, Bt + col*rowSizeBt, N);
#endif
    }
  }
}

template<class T>
MTL_INLINE static void SwapRows(T* a, I32 row1, I32 row2, I32 dataSize, I32 rowSize)
{
  for (I32 i = 0; i < dataSize; i++)
    Swap(a[row1 * rowSize + i], a[row2 * rowSize + i]);
}

template<class T> class DynamicMatrix;
template<class T> MTL_INLINE void ComputeTranspose(DynamicMatrix<T>& dst, const DynamicMatrix<T>& src);
template<class T> MTL_INLINE void ComputeTransposeParallel(DynamicMatrix<T>& dst, const DynamicMatrix<T>& src);

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

  MTL_INLINE DynamicMatrix(I32 rows, I32 cols, I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT)
  {
    Resize(rows, cols, byteAlignment);
  }

  MTL_INLINE DynamicMatrix(I32 rows, I32 cols, Initialize i, I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT)
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
  MTL_INLINE DynamicMatrix(const Matrix<M,N,T>& matrix, I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT)
  {
    Resize(M, N, byteAlignment);
    for (I32 row = 0; row < M; row++)
      memcpy((*this)[row], matrix[row], N*sizeof(T));
  }

  MTL_INLINE DynamicMatrix(I32 rows, I32 cols, I32 rowSize, T* pData)
  {
    Rows_    = rows;
    Cols_    = cols;
    RowSize_ = rowSize;
    pData_   = pData;
  }

  // Copy constructor.
  MTL_INLINE DynamicMatrix(const DynamicMatrix& rhs)
  {
    *this = rhs;
  }

  MTL_INLINE DynamicMatrix& operator=(const DynamicMatrix& rhs) 
  {
    if (this != &rhs)
    {
      Rows_    = rhs.Rows_;
      Cols_    = rhs.Cols_;
      RowSize_ = rhs.RowSize_;
      Data_    = rhs.Data_;

      if (Data_.Size() == 0 && Rows_ > 0)
      {
        Data_.Resize(Cols_ * RowSize_);
        OptimizedCopy(Data_.Begin(), rhs.pData_, Data_.Size());
      }
      pData_ = Data_.Begin();
    }

    return *this;
  }


  MTL_INLINE void Resize(I32 rows, I32 cols, I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT)
  {
    Rows_ = rows;
    Cols_ = cols;

    RowSize_ =
      ((Cols_ * (I32)sizeof(T) + byteAlignment - 1) / byteAlignment) * byteAlignment / (I32)sizeof(T);

    Data_.Resize(Rows_ * RowSize_);
    pData_ = Data_.Begin();
  }

  MTL_INLINE void Identity()
  {
    Zeros();
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] = T(1);
  }

  MTL_INLINE void Zeros()
  {
    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
      memset((*this)[i], 0, Cols() * sizeof(T));
  }

  MTL_INLINE void AddToDiagonals(const T& value)
  {
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] += value;
  }

  MTL_INLINE void SetDiagonals(const T& value)
  {
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] = value;
  }

  MTL_INLINE void SetDiagonals(const DynamicVector<T>& values)
  {
    for (I32 i = 0; i < Min(Rows(), Cols()); i++)
      (*this)[i][i] = values[i];
  }

  MTL_INLINE void SetAll(const T& newVal)
  {
    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
      AssignAll_Stream((*this)[i], newVal, Cols());
  }

  MTL_INLINE DynamicVector<T> operator*(const DynamicVector<T>& v) const
  {
    assert(Cols() == (I32)v.Size());

    DynamicVector<T> result(Rows());

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      result[i] = DotProduct_StreamAligned_Sequential((*this)[i], v.Begin(), Cols());
#else
      result[i] = DotProduct_Sequential((*this)[i], v.Begin(), Cols());
#endif

    return result;
  }

  MTL_INLINE DynamicMatrix operator*(const DynamicMatrix& B) const
  {
    assert(Cols() == B.Rows());

    DynamicMatrix result;
    Multiply(result, B);

    return result;
  }

  MTL_INLINE void Multiply(DynamicMatrix& product, const DynamicMatrix& B) const
  {
    DynamicMatrix Bt;
    B.ComputeTranspose(Bt);
    MultiplyTransposed(product, Bt);
  }

  MTL_INLINE void MultiplyTransposed(DynamicMatrix& product, const DynamicMatrix& Bt) const
  {
    assert(Cols() == Bt.Cols());

    product.Resize(Rows(), Bt.Rows());
    //MTL::MultiplyTransposed(product[0], (*this)[0], Bt[0], Rows(), Cols(), Bt.Rows(),
    //                        product.RowSize(), RowSize(), Bt.RowSize());
    product.Zeros();
    int blockSize = ComputeMultiplicationBlockSize();
    int offset = 0;
    while (offset < Cols())
    {
      int columnsToMultiply = Min(blockSize, Cols() - offset);
      AddMultiplyTransposed(product[0], (*this)[0] + offset, Bt[0] + offset, Rows(),
                            columnsToMultiply, Bt.Rows(),
                            product.RowSize(), RowSize(), Bt.RowSize());
      offset += columnsToMultiply;
    }
  }

  // Computes (*this) * this->getTranspose().
  MTL_INLINE DynamicMatrix MultiplyByTranspose() const
  {
    DynamicMatrix squareSymmetricMatrix;
    MultiplyByTranspose(squareSymmetricMatrix);
    return squareSymmetricMatrix;
  }
  MTL_INLINE void MultiplyByTranspose(DynamicMatrix& squareSymmetricMatrix) const
  {
    squareSymmetricMatrix.Resize(Rows(), Rows());
    squareSymmetricMatrix.Zeros();
    int blockSize = ComputeMultiplicationBlockSize();
    int offset = 0;
    while (offset < Cols())
    {
      int columnsToMultiply = Min(blockSize, Cols() - offset);
      AddMultiplyByTranspose(squareSymmetricMatrix[0], Data() + offset,
                             Rows(), columnsToMultiply,
                             squareSymmetricMatrix.RowSize(), RowSize());
      offset += columnsToMultiply;
    }
  }

  MTL_INLINE DynamicMatrix operator+=(const DynamicMatrix& B)
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      Addition_StreamAligned_Sequential((*this)[i], B[i], Cols());
#else
      Addition_Sequential((*this)[i], B[i], Cols());
#endif

    return *this;
  }
  MTL_INLINE DynamicMatrix operator+(const DynamicMatrix& B) const
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    DynamicMatrix C = *this;
    C += B;
    return C;
  }

  MTL_INLINE DynamicMatrix& operator-=(const DynamicMatrix& B)
  {
    assert(Cols() == B.Cols());
    assert(Rows() == B.Rows());

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      Subtraction_StreamAligned_Sequential((*this)[i], B[i], Cols());
#else
      Subtraction_Sequential((*this)[i], B[i], Cols());
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

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < A.Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      UnaryMinus_StreamAligned_Sequential(A[i], A.Cols());
#else
      UnaryMinus_Sequential(A[i], A.Cols());
#endif

    return A;
  }

  MTL_INLINE DynamicMatrix& operator+=(const T& scalar)
  {
    DynamicMatrix& A = *this;

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarAddition_StreamAligned_Sequential(A[i], scalar, A.Cols());
#else
    ScalarAddition_Sequential(A[i], scalar, A.Cols());
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
    DynamicMatrix& A = *this;

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarSubtraction_StreamAligned_Sequential(A[i], scalar, A.Cols());
#else
      ScalarSubtraction_Sequential(A[i], scalar, A.Cols());
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
    DynamicMatrix& A = *this;

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarMultiplication_StreamAligned_Sequential(A[i], scalar, A.Cols());
#else
      ScalarMultiplication_Sequential(A[i], scalar, A.Cols());
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
    DynamicMatrix& A = *this;

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarDivision_StreamAligned_Sequential(A[i], scalar, A.Cols());
#else
      ScalarDivision_Sequential(A[i], scalar, A.Cols());
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

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      sum += Sum_StreamAligned_Sequential((*this)[i], Cols());
#else
      sum += Sum_Sequential((*this)[i], Cols());
#endif

    return sum;
  }

  MTL_INLINE T SumOfSquares() const
  {
    T sum = T(0);

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      sum += SumOfSquares_StreamAligned_Sequential((*this)[i], Cols());
#else
      sum += SumOfSquares_Sequential((*this)[i], Cols());
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

    MTL_PARALLEL_FOR_BLOCKS(Rows())
    for (I32 i = 0; i < Rows(); i++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      T sum = SumOfAbsolutes_StreamAligned_Sequential((*this)[i], Cols());
#else
      T sum = SumOfAbsolutes_Sequential((*this)[i], Cols());
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

  MTL_INLINE DynamicMatrix ComputeTranspose(I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT) const
  {
    DynamicMatrix t;
    ComputeTranspose(t, byteAlignment);

    return t;
  }

  MTL_INLINE void ComputeTranspose(DynamicMatrix& t, I32 byteAlignment = MTL_DYNAMIC_MATRIX_DEFAULT_BYTE_ALIGNMENT) const
  {
    t.Resize(Cols(), Rows(), byteAlignment);
    MTL::ComputeTranspose(t, *this);
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
    return Data() + row * RowSize();
  }
  MTL_INLINE T* operator[](I32 row)
  {
    assert(row >= 0 && row < Rows());
    return Data() + row * RowSize();
  }

  MTL_INLINE I32 Rows() const              { return Rows_;    }
  MTL_INLINE I32 Cols() const              { return Cols_;    }
  MTL_INLINE I32 RowSize() const           { return RowSize_; }

  MTL_INLINE const T* Data() const         { return pData_;   }
  MTL_INLINE T* Data()                     { return pData_;   }

private:
  I32 Rows_;
  I32 Cols_;
  I32 RowSize_;
  DynamicVector<T> Data_;
  T* pData_;

  int ComputeMultiplicationBlockSize() const
  {
    return 256 * (Square(1024)/Square(Rows()) + 1);
  }
};

template<class T> MTL_INLINE void ComputeTranspose(DynamicMatrix<T>& dst, const DynamicMatrix<T>& src)
{
  for (I32 i = 0; i < src.Rows(); i++)
    for (I32 j = 0; j < src.Cols(); j++)
      dst[j][i] = src[i][j];
}
template<class T> MTL_INLINE void ComputeTransposeParallel(DynamicMatrix<T>& dst, const DynamicMatrix<T>& src)
{
  MTL_PARALLEL_FOR_BLOCKS(src.Rows())
  for (I32 i = 0; i < src.Rows(); i++)
    for (I32 j = 0; j < src.Cols(); j++)
      dst[j][i] = src[i][j];
}

}  // namespace MTL

#endif // MTL_DYNAMIC_MATRIX_H
