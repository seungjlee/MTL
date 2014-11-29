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


#ifndef MTL_MATRIX_H
#define MTL_MATRIX_H

#include "StreamMath.h"

namespace MTL
{

template<I32 M, I32 N, class T = F64>
class Matrix
{
  typedef T DataType[M][N];

public:
  MTL_INLINE Matrix() {}

  enum Initialize
  {
    eNothing,
    eIdentity,
    eZeros
  };

  MTL_INLINE Matrix(Initialize i)
  {
    switch (i)
    {
      case eIdentity:  Identity();    break;
      case eZeros:     Zeros();       break;
      default:                        break;
    }
  }

  // Constructor that sets all elements to initial.
  MTL_INLINE explicit Matrix(const T& initial)
  {
    SetAll(initial);
  }

  MTL_INLINE explicit Matrix(const T data[M][N])
  {
    memcpy(Data_, data, DataSizeInBytes());
  }

  // Sets matrix to be the identity matrix.
  MTL_INLINE void Identity()
  {
    Zeros();
    SetDiagonal<kMinimumDimension_>(Data_, T(1));
  }

  MTL_INLINE void Zeros()
  {
    memset(Data_, 0, DataSizeInBytes());
  }

  // Sets all elements in the matrix to newVal
  MTL_INLINE void SetAll(const T& newVal)
  {
    Set<M*N>(Data_[0], newVal);
  }

  MTL_INLINE void SetColumn(I32 col, const Matrix<M,1,T>& v)
  {
    assert(col >= 0 && col < N);
    SetColumn<M>(Data_, v.Data(), col);
  }

  MTL_INLINE void SetRow(I32 row, const Matrix<1,N,T>& v)
  {
    assert(row >= 0 && row < M);
    memcpy((*this)[row], v.data_[0], N*sizeof(T));
  }

  // Multiplication.
  template <I32 P>
  MTL_INLINE Matrix<M,P,T> operator*(const Matrix<N,P,T>& B) const
  {
    Matrix<M,P,T> result;
    for (I32 row = 0; row < M; row++)
      MultiplyRowByMatrix<N,P>::Compute<P>(result[row], Data_[row], B.Data());

    return result;
  }

  // Computes (*this) * this->ComputeTranspose().
  MTL_INLINE Matrix<M,M,T> MultiplyByTranspose() const
  {
    Matrix<M,M,T> squareSymmetricMatrix;

    for (I32 row = 0; row < M; row++)
    {
      for (I32 col = row; col < M; col++)
      {
        if (row == col)
        {
          squareSymmetricMatrix[row][col] = SumOfSquares<N>((*this)[row]);
        }
        else
        {
          squareSymmetricMatrix[row][col] = DotProduct<N>((*this)[row], (*this)[col]);
          squareSymmetricMatrix[col][row] = squareSymmetricMatrix[row][col];
        }
      }
    }

    return squareSymmetricMatrix;
  }

  // Computes this->ComputeTranspose() * (*this).
  MTL_INLINE Matrix<N,N,T> MultiplyTransposeByThis() const
  {
    Matrix<N,N,T> squareSymmetricMatrix;

    for (I32 row = 0; row < N; row++)
    {
      for (I32 col = row; col < N; col++)
      {
        if (row == col)
        {
          squareSymmetricMatrix[row][col] = ColumnSumOfSquares(col);
        }
        else
        {
          squareSymmetricMatrix[row][col] = DotProductOfColumns(row, col);
          squareSymmetricMatrix[col][row] = squareSymmetricMatrix[row][col];
        }
      }
    }

    return squareSymmetricMatrix;
  }

  // Adds matrix below this matrix.
  template <I32 P>
  MTL_INLINE Matrix<M+P,N,T> operator&&(const Matrix<P,N,T>& B) const
  {
    Matrix<M+P,N,T> matrix;
    memcpy(matrix[0],   Data_[0], M*N*sizeof(T));
    memcpy(matrix[M], B.Data_[0], P*N*sizeof(T));
    return matrix;
  }

  // Adds matrix to the right of this matrix.
  template <I32 P>
  MTL_INLINE Matrix<M,N+P,T> operator||(const Matrix<M,P,T>& B) const
  {
    Matrix<M,N+P,T> matrix;
    for (I32 row = 0; row < M; row++)
    {
      memcpy(matrix[row], (*this)[row], N*sizeof(T));
      memcpy(matrix[row] + N, B[row], P*sizeof(T));
    }
    return matrix;
  }

  template <I32 I, I32 J, I32 ROWS, I32 COLS>
  MTL_INLINE Matrix<ROWS,COLS> SubMatrix() const
  {
    assert(I >= 0 && I < M);
    assert(J >= 0 && J < N);
    assert(I + ROWS <= M);
    assert(J + COLS <= N);

    Matrix<ROWS,COLS,T> matrix;

    for (I32 row = 0; row < ROWS; row++)
      memcpy(matrix[row], (*this)[I + row] + J, COLS*sizeof(T));

    return matrix;
  }

  template<I32 COL>
  MTL_INLINE Matrix<M,1> Column() const
  {
    assert(COL >= 0 && COL < N);

    Matrix<M,1,T> columnVector;

    for (I32 row = 0; row < M; row++)
      columnVector[row][0] = (*this)[row][COL];

    return columnVector;
  }

  template<I32 ROW>
  MTL_INLINE Matrix<1,N> Row() const
  {
    assert(ROW >= 0 && ROW < M);

    Matrix<1,N,T> rowVector;

    memcpy(rowVector[0], (*this)[ROW], N*sizeof(T));

    return rowVector;
  }

  MTL_INLINE Matrix& operator+=(const Matrix& B)
  {
    Add<M*N>(Data_[0], B.Data_[0]);
    return *this;
  }
  MTL_INLINE Matrix operator+(const Matrix& B) const
  {
    Matrix C = *this;
    C += B;
    return C;
  }

  MTL_INLINE Matrix& operator-=(const Matrix& B)
  {
    Subtract<M*N>(Data_[0], B.Data_[0]);
    return *this;
  }
  MTL_INLINE Matrix operator-(const Matrix& B) const
  {
    Matrix C = *this;
    C -= B;
    return C;
  }

  MTL_INLINE Matrix operator-() const
  {
    Matrix A;
    UnaryMinus<M*N>(A.Data_[0], Data_[0]);
    return A;
  }

  MTL_INLINE Matrix& operator+=(const T& scalar)
  {
    ScalarAdd<M*N>(Data_[0], scalar);
    return *this;
  }
  MTL_INLINE Matrix operator+(const T& scalar) const
  {
    Matrix A = *this;
    A += scalar;
    return A;
  }

  MTL_INLINE Matrix& operator-=(const T& scalar)
  {
    ScalarSubtract<M*N>(Data_[0], scalar);
    return *this;
  }
  MTL_INLINE Matrix operator-(const T& scalar) const
  {
    Matrix A = *this;
    A -= scalar;
    return A;
  }

  MTL_INLINE Matrix& operator*=(const T& scalar)
  {
    ScalarMultiply<M*N>(Data_[0], scalar);
    return *this;
  }
  MTL_INLINE Matrix operator*(const T& scalar) const
  {
    Matrix A = *this;
    A *= scalar;
    return A;
  }

  MTL_INLINE Matrix& operator/=(const T& scalar)
  {
    ScalarDivide<M*N>(Data_[0], scalar);
    return *this;
  }
  MTL_INLINE Matrix operator/(const T& scalar) const
  {
    Matrix A = *this;
    A /= scalar;
    return A;
  }

  MTL_INLINE T Sum() const
  {
    return MTL::Sum<M*N>(Data_[0]);
  }

  // Sum of the squares of all elements.
  MTL_INLINE T SumOfSquares() const
  {
    return SumOfSquares<M*N>(Data_[0]);
  }

  // Returns maximum of all elements.
  MTL_INLINE T Max() const
  {
    return Maximum<M*N>(Data_[0]);
  }

  // Returns minimum of all elements.
  MTL_INLINE T Min() const
  {
    return Minimum<M*N>(Data_[0]);
  }

  MTL_INLINE T MaxOfAbsolutes() const
  {
    return Maximum<M*N>(Data_[0]);
  }
  MTL_INLINE T MinOfAbsolutes() const
  {
    return Minimum<M*N>(Data_[0]);
  }


  MTL_INLINE void RowMultiply(I32 row, const T& scalar)  { ScalarMultiply<N>(Data_[row], scalar); }
  MTL_INLINE void RowDivide  (I32 row, const T& scalar)  { ScalarDivide<N>  (Data_[row], scalar); }

  MTL_INLINE void ColumnMultiply(I32 col, const T& scalar)
  { ColumnScalarMultiply<M>(Data_[0] + col, scalar); }
  MTL_INLINE void ColumnDivide(I32 col, const T& scalar)
  { ColumnScalarDivide<M>(Data_[0] + col, scalar); }

  MTL_INLINE T ColumnSum(I32 col) const           { return ColumnSum<M>(col);          }
  MTL_INLINE T ColumnSumOfSquares(I32 col) const  { return ColumnSumOfSquares<M>(col); }
  MTL_INLINE T DotProductOfColumns(I32 col1, I32 col2) const
  { return DotProductOfColumns<M>(col1, col2); }

  MTL_INLINE F64 Mean() const
  {
    return Sum() / F64(M*N);
  }

  MTL_INLINE T FrobeniusNorm() const
  {
    return Sqrt(SumOfSquares());
  }

  MTL_INLINE Matrix<N,M,T> ComputeTranspose() const
  {
    Matrix<N,M,T> t;
    for (I32 i = 0; i < M; i++)
      for (I32 j = 0; j < N; j++)
        t[j][i] = Data_[i][j];

    return t;
  }

  MTL_INLINE bool operator==(const Matrix& B) const
  {
    return !memcmp(Data_, B.Data_, DataSizeInBytes());
  }
  MTL_INLINE bool operator!=(const Matrix& B) const
  {
    return !operator==(B);
  }

  MTL_INLINE bool IsZero() const      { return *this == Matrix(eZeros);    }
  MTL_INLINE bool IsIdentity() const  { return *this == Matrix(eIdentity); }

  MTL_INLINE const T* operator[](I32 row) const  { assert(row >= 0 && row < M); return Data_[row]; }
  MTL_INLINE T*       operator[](I32 row)        { assert(row >= 0 && row < M); return Data_[row]; }

  MTL_INLINE I32 Rows() const              { return M;             }
  MTL_INLINE I32 Cols() const              { return N;             }

  MTL_INLINE const DataType& Data() const  { return Data_;         }
  MTL_INLINE DataType& Data()              { return Data_;         }

  MTL_INLINE I32 DataSizeInBytes() const   { return sizeof(Data_); }

  template<I32 N> static T Maximum(const T* a)
  {
    return MTL::Max(Maximum<N-1>(a), a[N-1]);
  }
  template<> T static Maximum<1>(const T* a)
  {
    return a[0];
  }

  template<I32 N> static T Minimum(const T* a)
  {
    return MTL::Min(Minimum<N-1>(a), a[N-1]);
  }
  template<> T static Minimum<1>(const T* a)
  {
    return a[0];
  }

  template<I32 N> static T MaximumOfAbsolutes(const T* a)
  {
    return MTL::Max(MaximumOfAbsolutes<N-1>(a), Abs(a[N-1]));
  }
  template<> T static MaximumOfAbsolutes<1>(const T* a)
  {
    return Abs(a[0]);
  }

  template<I32 N> static T MinimumOfAbsolutes(const T* a)
  {
    return MTL::Min(MinimumOfAbsolutes<N-1>(a), Abs(a[N-1]));
  }
  template<> T static MinimumOfAbsolutes<1>(const T* a)
  {
    return Abs(a[0]);
  }

  template<I32 Q> static T SumOfSquares(const T* ptr)
  {
    return SumOfSquares<Q-1>(ptr) + Pow<2>(ptr[Q-1]);
  }
  template<> static T SumOfSquares<1>(const T* ptr)
  {
    return Pow<2>(ptr[0]);
  }

  template<I32 Q> static void UnaryMinus(T* a, const T* b)
  {
    UnaryMinus<Q-1>(a, b);
    a[Q-1] = -b[Q-1];
  }
  template<> static void UnaryMinus<1>(T* a, const T* b)
  {
    a[0] = -b[0];
  }

  template<I32 Q> static void Add(T* a, const T* b)
  {
    Add<Q-1>(a, b);
    a[Q-1] += b[Q-1];
  }
  template<> static void Add<1>(T* a, const T* b)
  {
    a[0] += b[0];
  }

  template<I32 Q> static void Subtract(T* a, const T* b)
  {
    Subtract<Q-1>(a, b);
    a[Q-1] -= b[Q-1];
  }
  template<> static void Subtract<1>(T* a, const T* b)
  {
    a[0] -= b[0];
  }

  template<I32 N> static void Multiply(T* a, const T* b)
  {
    Multiply<N-1>(a, b);
    a[N-1] *= b[N-1];
  }
  template<> static void Multiply<1>(T* a, const T* b)
  {
    a[0] *= b[0];
  }

  template<I32 N> static void Divide(T* a, const T* b)
  {
    Divide<N-1>(a, b);
    a[N-1] /= b[N-1];
  }
  template<> static void Divide<1>(T* a, const T* b)
  {
    a[0] /= b[0];
  }

  template<I32 Q> static void ScalarAdd(T* a, const T& scalar)
  {
    ScalarAdd<Q-1>(a, scalar);
    a[Q-1] += scalar;
  }
  template<> static void ScalarAdd<1>(T* a, const T& scalar)
  {
    a[0] += scalar;
  }

  template<I32 Q> static void ScalarSubtract(T* a, const T& scalar)
  {
    ScalarSubtract<Q-1>(a, scalar);
    a[Q-1] -= scalar;
  }
  template<> static void ScalarSubtract<1>(T* a, const T& scalar)
  {
    a[0] -= scalar;
  }

  template<I32 Q> static void ScalarMultiply(T* a, const T& scalar)
  {
    ScalarMultiply<Q-1>(a, scalar);
    a[Q-1] *= scalar;
  }
  template<> static void ScalarMultiply<1>(T* a, const T& scalar)
  {
    a[0] *= scalar;
  }

  template<I32 Q> static void ScalarDivide(T* a, const T& scalar)
  {
    ScalarDivide<Q-1>(a, scalar);
    a[Q-1] /= scalar;
  }
  template<> static void ScalarDivide<1>(T* a, const T& scalar)
  {
    a[0] /= scalar;
  }

protected:
  DataType Data_;

  template<I32 Q> void SetDiagonal(T ptr[M][N], const T& newVal)
  {
    SetDiagonal<Q-1>(ptr, newVal);
    ptr[Q-1][Q-1] = newVal;
  }
  template<> void SetDiagonal<1>(T ptr[M][N], const T& newVal)
  {
    ptr[0][0] = newVal;
  }

  template<I32 Q> void SetDiagonal(T ptr[M][N], const T newVals[])
  {
    SetDiagonal<Q-1>(ptr, newVals);
    ptr[Q-1][Q-1] = newVals[Q-1];
  }
  template<> void SetDiagonal<1>(T ptr[M][N], const T newVals[])
  {
    ptr[0][0] = newVals[0];
  }

  template<I32 Q> static T DotProduct(const T* a, const T* b)
  {
    return DotProduct<Q-1>(a, b) + a[Q-1] * b[Q-1];
  }
  template<> T static DotProduct<1>(const T* a, const T* b)
  {
    return a[0] * b[0];
  }

  template <I32 Cols>
  struct MultiplyRowByCol
  {
    template<I32 Q> static T Compute(const T* pRow, const T colData[][Cols], I32 col)
    {
      return Compute<Q-1>(pRow, colData, col) + pRow[Q-1] * colData[Q-1][col];
    }
    template<> static T Compute<1>(const T* pRow, const T colData[][Cols], I32 col)
    {
      return pRow[0] * colData[0][col];
    }
  };

  template <I32 Cols1, I32 Cols2>
  struct MultiplyRowByMatrix
  {
    template<I32 Q> static void Compute(T* pDst, const T* pRow, const T colData[][Cols2])
    {
      MultiplyRowByMatrix<Cols1,Cols2>::Compute<Q-1>(pDst, pRow, colData);
      pDst[Q-1] = MultiplyRowByCol<Cols2>::Compute<Cols1>(pRow, colData, Q-1);
    }
    template<> static void Compute<1>(T* pDst, const T* pRow, const T colData[][Cols2])
    {
      pDst[0] = MultiplyRowByCol<Cols2>::Compute<Cols1>(pRow, colData, 0);
    }
  };

  template<I32 N> static T Max_Norm(const T* a)
  {
    return Max(Max_Norm<N-1>(a), Abs(a[N-1]));
  }
  template<> T static Max_Norm<1>(const T* a)
  {
    return Abs(a[0]);
  }

private:
  enum { kMinimumDimension_ = M < N ? M : N };

  template<I32 Q> static void Set(T* ptr, const T& newVal)
  {
    Set<Q-1>(ptr, newVal);
    ptr[Q-1] = newVal;
  }
  template<> static void Set<1>(T* ptr, const T& newVal)
  {
    ptr[0] = newVal;
  }

  template<I32 Q> static void SetColumn(T dst[M][N], const T src[M][1], I32 col)
  {
    SetColumn<Q-1>(dst, src, col);
    dst[Q-1][col] = src[Q-1][0];
  }
  template<> static void SetColumn<1>(T dst[M][N], const T src[M][1], I32 col)
  {
    dst[0][col] = src[0][0];
  }

  template<I32 Q> static void ColumnScalarMultiply(T* a, const T& scalar)
  {
    *a *= scalar;
    ColumnScalarMultiply<Q-1>(a + N, scalar);
  }
  template<> static void ColumnScalarMultiply<1>(T* a, const T& scalar)
  {
    *a *= scalar;
  }

  template<I32 Q> static void ColumnScalarDivide(T* a, const T& scalar)
  {
    *a /= scalar;
    ColumnScalarDivide<Q-1>(a + N, scalar);
  }
  template<> static void ColumnScalarDivide<1>(T* a, const T& scalar)
  {
    *a /= scalar;
  }

  template<I32 Q> T ColumnSum(I32 col) const
  {
    return ColumnSum<Q-1>(col) + Data_[Q-1][col];
  }
  template<> T ColumnSum<1>(I32 col) const
  {
    return Data_[0][col];
  }

  template<I32 Q> T ColumnSumOfSquares(I32 col) const
  {
    return ColumnSumOfSquares<Q-1>(col) + Pow<2>(Data_[Q-1][col]);
  }
  template<> T ColumnSumOfSquares<1>(I32 col) const
  {
    return Pow<2>(Data_[0][col]);
  }

  template<I32 Q> T DotProductOfColumns(I32 col1, I32 col2) const
  {
    return DotProductOfColumns<Q-1>(col1, col2) + Data_[Q-1][col1] * Data_[Q-1][col2];
  }
  template<> T DotProductOfColumns<1>(I32 col1, I32 col2) const
  {
    return Data_[0][col1] * Data_[0][col2];
  }
};


#define MTL_MATRIX_COMMON_DEFINITIONS(Class, Base, M, N, T)                                   \
  MTL_INLINE Class(const Base& v) : Base(v) {}                                                \
  MTL_INLINE Class& operator=(const Base& v)  { Base::operator=(v); return *this; }           \
  MTL_INLINE Class(Initialize i) : Base(i) {}                                                 \
  MTL_INLINE explicit Class(const T& initialVal) : Base(initialVal) {}                        \
  MTL_INLINE explicit Class(const T data[M][N]) : Base(data) {}                               \
  MTL_INLINE Class& operator+=(const Base& b)  { Base::operator+=(b); return *this; }         \
  MTL_INLINE Class& operator-=(const Base& b)  { Base::operator-=(b); return *this; }         \
  MTL_INLINE Class& operator+=(const Class& b) { Base::operator+=((Base&)b); return *this; }  \
  MTL_INLINE Class& operator-=(const Class& b) { Base::operator-=((Base&)b); return *this; }  \
  MTL_INLINE Class operator-() const  { return Base::operator-(); }                           \
  MTL_INLINE Class operator+(const Base& b) const  { return Base::operator+(b); }             \
  MTL_INLINE Class operator-(const Base& b) const  { return Base::operator-(b); }             \
  MTL_INLINE Class operator+(const Class& b) const  { return Base::operator+((Base&)b); }     \
  MTL_INLINE Class operator-(const Class& b) const  { return Base::operator-((Base&)b); }     \
  MTL_INLINE Class& operator+=(const T& scalar)  { Base::operator+=(scalar); return *this; }  \
  MTL_INLINE Class& operator-=(const T& scalar)  { Base::operator-=(scalar); return *this; }  \
  MTL_INLINE Class& operator*=(const T& scalar)  { Base::operator*=(scalar); return *this; }  \
  MTL_INLINE Class& operator/=(const T& scalar)  { Base::operator/=(scalar); return *this; }  \
  MTL_INLINE Class operator+(const T& scalar) const  { return Base::operator+(scalar); }      \
  MTL_INLINE Class operator-(const T& scalar) const  { return Base::operator-(scalar); }      \
  MTL_INLINE Class operator*(const T& scalar) const  { return Base::operator*(scalar); }      \
  MTL_INLINE Class operator/(const T& scalar) const  { return Base::operator/(scalar); }      \
  MTL_INLINE bool operator==(const Base& b) const  { return Base::operator==(b); }            \
  MTL_INLINE bool operator!=(const Base& b) const  { return Base::operator!=(b); }            \
  MTL_INLINE bool operator==(const Class& b) const  { return Base::operator==((Base&)b); }    \
  MTL_INLINE bool operator!=(const Class& b) const  { return Base::operator!=((Base&)b); }    \
  template <I32 P> MTL_INLINE Matrix<M,P,T> operator*(const Matrix<N,P,T>& b) const           \
  { return Base::operator*(b); }                                                              \

#define MTL_SQUARE_MATRIX_COMMON_DEFINITIONS(Class, Base, N, T)                               \
  MTL_MATRIX_COMMON_DEFINITIONS(Class, Base, N, N, T)                                         \
  MTL_INLINE Class operator*(const Base& b) const   { return Base::operator*(b);        }     \
  MTL_INLINE Class operator*(const Class& b) const  { return Base::operator*((Base&)b); }     \
  MTL_INLINE Class& operator*=(const Base& b)       { return (*this = *this * b);       }     \
  MTL_INLINE Class& operator*=(const Class& b)      { return (*this = *this * b);       }     \

#define MTL_COLUMN_VECTOR_COMMON_DEFINITIONS(Class, Base, M, T)                               \
  MTL_MATRIX_COMMON_DEFINITIONS(Class, Base, M, 1, T)                                         \
  MTL_INLINE explicit Class(const T data[M]) : Base(data) {}                                  \
  MTL_INLINE Class& operator*=(const Base& v)  { Base::operator*=(v); return *this; }         \
  MTL_INLINE Class& operator/=(const Base& v)  { Base::operator/=(v); return *this; }         \
  MTL_INLINE Class operator*(const Base& v) const  { return Base::operator*(v); }             \
  MTL_INLINE Class operator/(const Base& v) const  { return Base::operator/(v); }             \
  MTL_INLINE Class& operator*=(const Class& v)  {Base::operator*=((Base&)v); return *this;}   \
  MTL_INLINE Class& operator/=(const Class& v)  {Base::operator/=((Base&)v); return *this;}   \
  MTL_INLINE Class operator*(const Class& v) const  { return Base::operator*((Base&)v); }     \
  MTL_INLINE Class operator/(const Class& v) const  { return Base::operator/((Base&)v); }     \


template<I32 M, class T = F64>
class ColumnVector : public Matrix<M,1,T>
{
  typedef Matrix<M,1,T> Base;

public:
  MTL_MATRIX_COMMON_DEFINITIONS(ColumnVector, Base, M, 1, T);

  MTL_INLINE ColumnVector() : Base() {}

  MTL_INLINE explicit ColumnVector(const T data[M])
  {
    memcpy(Data_, data, DataSizeInBytes());
  }

  MTL_INLINE ColumnVector& operator*=(const ColumnVector& v)
  {
    Multiply<M>(Data_[0], v.Data_[0]);
    return *this;
  }
  MTL_INLINE ColumnVector& operator*=(const Base& v)
  {
    return operator*=(ColumnVector(v));
  }
  MTL_INLINE ColumnVector operator*(const Base& v2) const
  {
    ColumnVector<M,T> v3 = *this;
    v3 *= v2;
    return v3;
  }

  MTL_INLINE ColumnVector& operator/=(const ColumnVector& v)
  {
    Divide<M>(Data_[0], v.Data_[0]);
    return *this;
  }
  MTL_INLINE ColumnVector& operator/=(const Base& v)
  {
    return operator/=(ColumnVector(v));
  }
  MTL_INLINE ColumnVector<M,T> operator/(const Base& v2) const
  {
    ColumnVector<M,T> v3 = *this;
    v3 /= v2;
    return v3;
  }

  MTL_INLINE T Dot(const Base& v2) const
  {
    return DotProduct<M>(Data_[0], v2.Data_[0]);
  }

  MTL_INLINE T MaxNorm() const
  {
    return Max_Norm<M>(Data_[0]);
  }

  MTL_INLINE void Normalize()
  {
    *this /= Sqrt(SumOfSquares());
  }

  MTL_INLINE T& operator[](I32 i)               { assert(i >= 0 && i < M);  return Data_[i][0]; }
  MTL_INLINE const T& operator[](I32 i) const   { assert(i >= 0 && i < M);  return Data_[i][0]; }

  MTL_INLINE I32 Size() const  { return Rows();}
};

typedef ColumnVector<1,F64> ColumnVector1D;
typedef ColumnVector<2,F64> ColumnVector2D;
typedef ColumnVector<3,F64> ColumnVector3D;
typedef ColumnVector<4,F64> ColumnVector4D;
typedef ColumnVector<5,F64> ColumnVector5D;
typedef ColumnVector<6,F64> ColumnVector6D;


MTL_INLINE static ColumnVector3D Cross(const ColumnVector3D& u, const ColumnVector3D& v)
{
  ColumnVector<3,F64> crossProduct;
  crossProduct[0] = u[1] * v[2] - u[2] * v[1];
  crossProduct[1] = u[2] * v[0] - u[0] * v[2];
  crossProduct[2] = u[0] * v[1] - u[1] * v[0];
  return crossProduct; 
}


// Solves A*X = B where LU is the matrix containg the L and U factors and P is the
// permutation vector from the LUP factorization of A.
// The solution X is returned in B.
template <I32 N, I32 C, class T>
MTL_INLINE static void SolveLUP(Matrix<N,C,T>& B,
                                const Matrix<N,N,T>& LU, const ColumnVector<N,I32>& P)
{
  ColumnVector<N,T> Y;

  // Solve for each column.
  for (I32 col = 0; col < C; col++)
  {
    // Solve L*Y = B first.
    for (I32 row = 0; row < N; row++)
    {
      Y[row] = B.Data()[P[row]][col];

      for (I32 i = 0; i < row; i++)
        Y[row] -= LU.Data()[P[row]][i] * Y[i];
    }

    // Now solve U*B = Y,
    for (I32 row = N-1; row >= 0; row--)
    {
      for (I32 i = row + 1; i < N; i++)
        Y[row] -= LU.Data()[P[row]][i] * Y[i];

      Y[row] /= LU.Data()[P[row]][row];
      B.Data()[row][col] = Y[row];
    }
  }
}


template<I32 N, class T = F64>
class SquareMatrix : public Matrix<N,N,T>
{
  typedef Matrix<N,N,T> Base;

public:
  MTL_SQUARE_MATRIX_COMMON_DEFINITIONS(SquareMatrix, Base, N, T);

  MTL_INLINE SquareMatrix() : Base() {}
  MTL_INLINE explicit SquareMatrix(const T diagonalValues[N]) : Base()
  {
    Zeros();
    setDiagonals(diagonalValues);
  }

  MTL_INLINE ColumnVector<N,T> operator*(const ColumnVector<N,T>& b) const
  { return Base::operator*(b); }

  MTL_INLINE void SetDiagonals(const T& value)
  {
    SetDiagonal<N>(Data_, value);
  }

  MTL_INLINE void SetDiagonals(const T values[N])
  {
    SetDiagonal<N>(Data_, values);
  }

  // Transposes this matrix.
  MTL_INLINE void Transpose()
  {
    for (I32 i = 0; i < N; i++)
      for (I32 j = i + 1; j < N; j++)
        Swap(Data_[i][j], Data_[j][i]);
  }

  MTL_INLINE bool IsSingular() const  { return determinant() == 0; }

  // Returns determinant of the matrix. Computes determinant using LUP factorization.
  MTL_INLINE T Determinant() const
  {
    SquareMatrix LU;
    ColumnVector<N,I32> P;
    T sign;

    if (LUP(LU, P, sign))
      return Determinant(LU, P, sign);
    else
      return 0;
  }
  MTL_INLINE static T Determinant(const SquareMatrix& LU, const ColumnVector<N,I32>& P,
                                  const T& sign)
  {
    T product = LU.Data_[P[0]][0];
    for (I32 i = 1; i < N; i++)
      product *= LU.Data_[P[i]][i];

    return sign * product;
  }

  MTL_INLINE SquareMatrix ComputeTranspose() const
  {
    return Base::ComputeTranspose();
  }

  MTL_INLINE SquareMatrix Inverse() const
  {
    SquareMatrix inverse;
    inverse.Identity();

    if (SolveLUP(inverse))
      return inverse;
    else
      return SquareMatrix(T(kNAN));
  }

  // Returns L and U matrices compressed in one matrix. The diagonals are all U values since
  // L's diagonals are all 1's. Returns true if successful. Returns false if matrix is singular
  // for the computational precision. Returns the permutation vector P and its corresponding sign.
  MTL_INLINE bool LUP(SquareMatrix& LU, ColumnVector<N,I32>& P, T& sign) const
  {
    LU = *this;
    sign = 1;

    ColumnVector<N,T> scale;
    for (I32 row = 0; row < N; row++)
    {
      P[row] = row;

      scale[row] = Abs(LU.Data_[row][1]);

      for (I32 col = 2; col < N; col++)
      {
        T temp = Abs(LU.Data_[row][col]);
        if (temp > scale[row])
          scale[row] = temp;
      }
    }

    for (I32 row = 0; row < N-1; row++)
    {
      I32 pivot = row;

      // Find max for pivot row.
      T max = Abs(LU.Data_[P[pivot]][row]) / scale[P[pivot]];
      for (I32 i = row + 1; i < N; i++)
      {
        T temp = Abs(LU.Data_[P[i]][row]) / scale[P[i]];
        if (temp > max)
        {
          max = temp;
          pivot = i;
        }
      }

      if (pivot != row)
      {
        Swap(P[pivot], P[row]);
        sign = -sign;
      }

      if (Abs(LU.Data_[P[row]][row]) < Epsilon<T>())
        return false;

      // Compute L and U.
      for (I32 i = row + 1; i < N; i++)
      {
        T temp = LU.Data_[P[i]][row] / LU.Data_[P[row]][row];
        LU.Data_[P[i]][row] = temp;

        for (I32 j = row + 1; j < N; j++)
          LU.Data_[P[i]][j] -= temp * LU.Data_[P[row]][j];
      }
    }

    return true;
  }

  // Solves A * X = B where A is this matrix. Uses LUP factorization to solve the system of linear
  // equations. Replaces B with the result X if successful. Returns true if successful.
  // Returns false if matrix is singular.
  template<I32 C>
  MTL_INLINE bool SolveLUP(Matrix<N,C,T>& B) const
  {
    SquareMatrix LU;
    ColumnVector<N,I32> P;
    T sign;

    if (LUP(LU, P, sign))
    {
      MTL::SolveLUP<N,C,T>(B, LU, P);
      return true;
    }
    else
      return false;
  }

  // Computes the condition number.
  MTL_INLINE T ConditionNumber() const
  {
    SquareMatrix LU;
    ColumnVector<N,I32> P;
    T sign;

    if (factorizeLUP(LU, P, sign))
    {
      return conditionNumber(*this, LU, P);
    }
    else
      return T(kINF);  // Return infinity.
  }
  MTL_INLINE static T ConditionNumber(const SquareMatrix<N,T>& M, const SquareMatrix<N,T>& LU,
                                      const ColumnVector<N,I32>& P)
  {
    SquareMatrix inverse;
    inverse.Identity();
    SolveLUP<N,N,T>(inverse, LU, P);

    return M.FrobeniusNorm() * inverse.FrobeniusNorm();
  }

private:
  template<I32 Q> static T ProductOfDiagonals(T ptr[N][N])
  {
    return ProductOfDiagonals<Q-1>(ptr) * ptr[Q-1][Q-1];
  }
  template<> static T ProductOfDiagonals<1>(T ptr[N][N])
  {
    return ptr[0][0];
  }
};

typedef SquareMatrix<1,F64> SquareMatrix1x1;
typedef SquareMatrix<2,F64> SquareMatrix2x2;
typedef SquareMatrix<3,F64> SquareMatrix3x3;
typedef SquareMatrix<4,F64> SquareMatrix4x4;
typedef SquareMatrix<5,F64> SquareMatrix5x5;
typedef SquareMatrix<6,F64> SquareMatrix6x6;

MTL_INLINE void SquareMatrix1x1::Transpose()
{
}
MTL_INLINE void SquareMatrix2x2::Transpose()
{
  Swap(Data_[0][1], Data_[1][0]);
}
MTL_INLINE void SquareMatrix3x3::Transpose()
{
  Swap(Data_[0][1], Data_[1][0]);
  Swap(Data_[0][2], Data_[2][0]);
  Swap(Data_[1][2], Data_[2][1]);
}
MTL_INLINE void SquareMatrix4x4::Transpose()
{
  Swap(Data_[0][1], Data_[1][0]);
  Swap(Data_[0][2], Data_[2][0]);
  Swap(Data_[0][3], Data_[3][0]);
  Swap(Data_[1][2], Data_[2][1]);
  Swap(Data_[1][3], Data_[3][1]);
  Swap(Data_[2][3], Data_[3][2]);
}

MTL_INLINE F64 SquareMatrix1x1::Determinant() const
{
  return Data_[0][0];
}
MTL_INLINE F64 SquareMatrix2x2::Determinant() const
{
  return Data_[0][0] * Data_[1][1] - Data_[0][1] * Data_[1][0];
}
MTL_INLINE F64 SquareMatrix3x3::Determinant() const
{
  return (Data_[0][0] * (Data_[1][1] * Data_[2][2] - Data_[1][2] * Data_[2][1]) +
          Data_[0][1] * (Data_[1][2] * Data_[2][0] - Data_[1][0] * Data_[2][2]) +
          Data_[0][2] * (Data_[1][0] * Data_[2][1] - Data_[1][1] * Data_[2][0]));
}

MTL_INLINE SquareMatrix1x1 SquareMatrix1x1::ComputeTranspose() const
{
  return *this;
}
MTL_INLINE SquareMatrix2x2 SquareMatrix2x2::ComputeTranspose() const
{
  SquareMatrix2x2 t;
  t[0][0] = Data_[0][0];
  t[0][1] = Data_[1][0];
  t[1][0] = Data_[0][1];
  t[1][1] = Data_[1][1];
  return t;
}
MTL_INLINE SquareMatrix3x3 SquareMatrix3x3::ComputeTranspose() const
{
  SquareMatrix3x3 t;
  t[0][0] = Data_[0][0];
  t[0][1] = Data_[1][0];
  t[0][2] = Data_[2][0];
  t[1][0] = Data_[0][1];
  t[1][1] = Data_[1][1];
  t[1][2] = Data_[2][1];
  t[2][0] = Data_[0][2];
  t[2][1] = Data_[1][2];
  t[2][2] = Data_[2][2];
  return t;
}

MTL_INLINE SquareMatrix4x4 SquareMatrix4x4::ComputeTranspose() const
{
  SquareMatrix4x4 t;
  t[0][0] = Data_[0][0];
  t[0][1] = Data_[1][0];
  t[0][2] = Data_[2][0];
  t[0][3] = Data_[3][0];
  t[1][0] = Data_[0][1];
  t[1][1] = Data_[1][1];
  t[1][2] = Data_[2][1];
  t[1][3] = Data_[3][1];
  t[2][0] = Data_[0][2];
  t[2][1] = Data_[1][2];
  t[2][2] = Data_[2][2];
  t[2][3] = Data_[3][2];
  t[3][0] = Data_[0][3];
  t[3][1] = Data_[1][3];
  t[3][2] = Data_[2][3];
  t[3][3] = Data_[3][3];
  return t;
}

template <I32 N>
MTL_INLINE static SquareMatrix<N,F64> Inverse(const SquareMatrix<N,F64>& M, F64)
{
  return M.inverse();
}

template<>
MTL_INLINE static SquareMatrix1x1 Inverse(const SquareMatrix1x1& M, F64 determinant)
{
  return SquareMatrix1x1(1. / determinant);
}

MTL_INLINE static SquareMatrix2x2 Inverse(const SquareMatrix2x2& M, F64 determinant)
{
  F64 reciprocalDet = 1. / determinant;
  SquareMatrix2x2 inv;
  inv[0][0] =  M[1][1] * reciprocalDet;
  inv[0][1] = -M[0][1] * reciprocalDet;
  inv[1][0] = -M[1][0] * reciprocalDet;
  inv[1][1] =  M[0][0] * reciprocalDet;
  return inv;
}

MTL_INLINE static SquareMatrix3x3 Inverse(const SquareMatrix3x3& M, F64 determinant)
{
  F64 reciprocalDet = 1. / determinant;
  SquareMatrix3x3 inv;

  inv[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[2][1]) * reciprocalDet;
  inv[0][1] = (M[0][2]*M[2][1]-M[0][1]*M[2][2]) * reciprocalDet;
  inv[0][2] = (M[0][1]*M[1][2]-M[0][2]*M[1][1]) * reciprocalDet;

  inv[1][0] = (M[1][2]*M[2][0]-M[1][0]*M[2][2]) * reciprocalDet;
  inv[1][1] = (M[0][0]*M[2][2]-M[0][2]*M[2][0]) * reciprocalDet;
  inv[1][2] = (M[0][2]*M[1][0]-M[0][0]*M[1][2]) * reciprocalDet;

  inv[2][0] = (M[1][0]*M[2][1]-M[1][1]*M[2][0]) * reciprocalDet;
  inv[2][1] = (M[0][1]*M[2][0]-M[0][0]*M[2][1]) * reciprocalDet;
  inv[2][2] = (M[0][0]*M[1][1]-M[0][1]*M[1][0]) * reciprocalDet;

  return inv;
}

MTL_INLINE SquareMatrix1x1 SquareMatrix1x1::Inverse() const
{
  return SquareMatrix1x1(1./Data_[0][0]);
}
MTL_INLINE SquareMatrix2x2 SquareMatrix2x2::Inverse() const
{
  return MTL::Inverse(*this, Determinant());
}
MTL_INLINE SquareMatrix3x3 SquareMatrix3x3::Inverse() const
{
  return MTL::Inverse(*this, Determinant());
}

MTL_INLINE SquareMatrix3x3 SquareMatrix3x3::operator*(const SquareMatrix3x3& B) const
{
  SquareMatrix3x3 product;
  MultiplyRowByMatrix<3,3>::Compute<3>(product[0], Data_[0], B.Data());
  MultiplyRowByMatrix<3,3>::Compute<3>(product[1], Data_[1], B.Data());
  MultiplyRowByMatrix<3,3>::Compute<3>(product[2], Data_[2], B.Data());

  return product;
}

MTL_INLINE ColumnVector3D SquareMatrix3x3::operator*(const ColumnVector3D& B) const
{
  ColumnVector3D product;
  product[0] = MultiplyRowByCol<1>::Compute<3>(Data_[0], B.Data(), 0);
  product[1] = MultiplyRowByCol<1>::Compute<3>(Data_[1], B.Data(), 0);
  product[2] = MultiplyRowByCol<1>::Compute<3>(Data_[2], B.Data(), 0);

  return product;
}

MTL_INLINE SquareMatrix2x2 SquareMatrix2x2::operator*(const SquareMatrix2x2& B) const
{
  SquareMatrix2x2 product;
  MultiplyRowByMatrix<2,2>::Compute<2>(product[0], Data_[0], B.Data());
  MultiplyRowByMatrix<2,2>::Compute<2>(product[1], Data_[1], B.Data());

  return product;
}

MTL_INLINE ColumnVector2D SquareMatrix2x2::operator*(const ColumnVector2D& B) const
{
  ColumnVector2D product;
  product[0] = MultiplyRowByCol<1>::Compute<2>(Data_[0], B.Data(), 0);
  product[1] = MultiplyRowByCol<1>::Compute<2>(Data_[1], B.Data(), 0);

  return product;
}

MTL_INLINE SquareMatrix3x3 ComputeCrossProductMatrix(const ColumnVector3D& v)
{
  SquareMatrix3x3 m;
  m[0][0] =     0;  m[0][1] = -v[2];  m[0][2] =   v[1];
  m[1][0] =  v[2];  m[1][1] =     0;  m[1][2] =  -v[0];
  m[2][0] = -v[1];  m[2][1] =  v[0];  m[2][2] =      0;

  return m;
}

// Solves A * x = b by computing the inverse of A. Returns true if the inverse can be computed,
// false otherwise. The inverse is computed using the determinant of the matrix.
template<I32 N>
MTL_INLINE static bool SolveWithInverse(ColumnVector<N>& x,
                                        const SquareMatrix<N> A, const ColumnVector<N>& b)
{
  T determinant = A.determinant();
  if (determinant == 0)
    return false;

  x = Inverse(A, determinant) * b;

  return true;
}

// Recursive implementation of matrix determinant computation.
template<I32 N, class T> MTL_INLINE static T DeterminantRecursive(const SquareMatrix<N,T>& A)
{
  T sum = T(0);
  T sign = T(1);

  for (I32 i = 0; i < N; i++)
  {
    T m[N-1][N-1];
    for (I32 row = 0; row < N-1; row++)
      for (I32 j = 0, k = 0; j < N; j++)
        if (i != j)
          m[row][k++] = A[row+1][j];

    SquareMatrix<N-1,T> M(m);
    sum += sign * A[0][i] * DeterminantRecursive(M);
    sign *= T(-1);
  }

  return sum;
}
template<class T> MTL_INLINE static T DeterminantRecursive(const SquareMatrix<1,T>& A)
{
  return A[0][0];
}

// Recursive implementation of matrix inverse computation.
template<I32 N, class T>
MTL_INLINE static SquareMatrix<N,T> InverseRecursive(const SquareMatrix<N,T>& A)
{
  SquareMatrix<N,T> inv;
  T reciprocalDet = 1. / DeterminantRecursive(A);

  for (I32 i = 0; i < N; i++)
  {
    T sign = i & 1 ? T(-1) : T(1);
    for (I32 j = 0; j < N; j++)
    {
      T m[N-1][N-1];
      for (I32 k = 0, p = 0; k < N; k++)
      {
        if (i != k)
        {
          for (I32 l = 0, q = 0; l < N; l++)
            if (j != l)
              m[p][q++] = A.Data()[k][l];
          p++;
        }
      }

      SquareMatrix<N-1,T> M(m);
      inv.Data()[j][i] = sign * DeterminantRecursive(M) * reciprocalDet;
      sign *= T(-1);
    }
  }

  return inv;
}

}  // namespace MTL

#endif // MTL_MATRIX_H
