//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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

#include <assert.h>
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
    Array2D<M,N,T,kMinimumDimension_>::SetDiagonal(Data_, T(1));
  }

  MTL_INLINE void Zeros()
  {
    memset(Data_, 0, DataSizeInBytes());
  }

  MTL_INLINE void SetDiagonals(const T& value)
  {
    Array2D<M,N,T,kMinimumDimension_>::SetDiagonal(Data_, value);
  }

  MTL_INLINE void SetDiagonals(const T values[])
  {
    Array2D<M,N,T,kMinimumDimension_>::SetDiagonal(Data_, values);
  }

  MTL_INLINE void AddToDiagonals(const T& value)
  {
    Array2D<M,N,T,kMinimumDimension_>::AddToDiagonal(Data_, value);
  }

  // Sets all elements in the matrix to newVal
  MTL_INLINE void SetAll(const T& newVal)
  {
    Array<M*N,T>::Set(Data_[0], newVal);
  }

  MTL_INLINE void SetColumn(I32 col, const Matrix<M,1,T>& v)
  {
    assert(col >= 0 && col < N);
    Array2D<M,N,T,M>::SetColumn(Data(), v.Data(), col);
  }

  MTL_INLINE void SetRow(I32 row, const Matrix<1,N,T>& v)
  {
    assert(row >= 0 && row < M);
    memcpy((*this)[row], v.Data()[0], N*sizeof(T));
  }

  // Multiplication.
  template <I32 P>
  MTL_INLINE Matrix<M,P,T> operator*(const Matrix<N,P,T>& B) const
  {
    Matrix<M,P,T> result;
    for (I32 row = 0; row < M; row++)
      Array2D<N,P,T,P>::MultiplyRowByMatrix(result[row], Data_[row], B.Data());

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
          squareSymmetricMatrix[row][col] = MTL::SumOfSquares<N>((*this)[row]);
        }
        else
        {
          squareSymmetricMatrix[row][col] = MTL::Dot<N>((*this)[row], (*this)[col]);
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
    memcpy(matrix[0],   Data()[0], M*N*sizeof(T));
    memcpy(matrix[M], B.Data()[0], P*N*sizeof(T));
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
    Array<M*N,T>::Add(Data()[0], B.Data()[0]);
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
    Array<M*N,T>::Subtract(Data()[0], B.Data()[0]);
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
    Array<M*N,T>::UnaryMinus(A.Data()[0], Data()[0]);
    return A;
  }

  MTL_INLINE Matrix& operator+=(const T& scalar)
  {
    Array<M*N,T>::AddScalar(Data()[0], scalar);
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
    Array<M*N,T>::SubtractScalar<M*N>(Data()[0], scalar);
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
    Array<M*N,T>::MultiplyScalar(Data()[0], scalar);
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
    Array<M*N,T>::DivideScalar(Data()[0], scalar);
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
    return MTL::Sum<M*N>(Data()[0]);
  }

  // Sum of the squares of all elements.
  MTL_INLINE T SumOfSquares() const
  {
    return MTL::SumOfSquares<M*N>(Data()[0]);
  }

  // Returns maximum of all elements.
  MTL_INLINE T Max() const
  {
    return MTL::Maximum<M*N>(Data()[0]);
  }

  // Returns minimum of all elements.
  MTL_INLINE T Min() const
  {
    return MTL::Minimum<M*N>(Data()[0]);
  }

  MTL_INLINE T MaxOfAbsolutes() const
  {
    return MTL::MaximumOfAbsolutes<M*N>(Data()[0]);
  }
  MTL_INLINE T MinOfAbsolutes() const
  {
    return MTL::MinimumOfAbsolutes<M*N>(Data()[0]);
  }

  MTL_INLINE void RowMultiply(I32 row, const T& scalar)  { ScalarMultiply<N>(Data()[row], scalar); }
  MTL_INLINE void RowDivide  (I32 row, const T& scalar)  { ScalarDivide<N>  (Data()[row], scalar); }

  MTL_INLINE void ColumnMultiply(I32 col, const T& scalar)
  { Array_2D<N,T,M>::ColumnScalarMultiply(Data()[0] + col, scalar); }
  MTL_INLINE void ColumnDivide(I32 col, const T& scalar)
  { Array_2D<N,T,M>::ColumnScalarDivide(Data()[0] + col, scalar); }

  MTL_INLINE T ColumnSum(I32 col) const  { return Array_2D<N,T,M>::ColumnSum(Data()[0] + col); }
  MTL_INLINE T ColumnSumOfSquares(I32 col) const
  { return Array_2D<N,T,M>::ColumnSumOfSquares(Data()[0] + col); }
  MTL_INLINE T DotProductOfColumns(I32 col1, I32 col2) const
  { return Array_2D<N,T,M>::DotProductOfColumns(Data()[0] + col1, Data()[0] + col2); }

  MTL_INLINE T Mean() const
  {
    return Sum() / T(M*N);
  }

  MTL_INLINE T FrobeniusNorm() const
  {
    return Sqrt(SumOfSquares());
  }

  MTL_INLINE T RMS() const
  {
    return Sqrt(SumOfSquares() / T(M*N));
  }

  MTL_INLINE Matrix<N,M,T> ComputeTranspose() const
  {
    Matrix<N,M,T> t;
    for (I32 i = 0; i < M; i++)
      for (I32 j = 0; j < N; j++)
        t[j][i] = Data()[i][j];

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

protected:
  DataType Data_;

private:
  enum { kMinimumDimension_ = M < N ? M : N };
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
  template <MTL::I32 P> MTL_INLINE Matrix<M,P,T> operator*(const Matrix<N,P,T>& b) const      \
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
    Multiply<M>(Data()[0], v.Data()[0]);
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
    Divide<M>(Data()[0], v.Data()[0]);
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

  MTL_INLINE T Dot(const Base& v) const
  {
    return MTL::Dot<M>(Data()[0], v.Data()[0]);
  }

  MTL_INLINE T MaxNorm() const
  {
    return Array<M,T>::MaxNorm(Data()[0]);
  }

  MTL_INLINE void Normalize()
  {
    *this /= Sqrt(SumOfSquares());
  }

  template <I32 OFFSET, I32 ROWS>
  MTL_INLINE ColumnVector<ROWS,T> SubColumn() const
  {
    return SubMatrix<OFFSET,0,ROWS,1>();
  }  

  MTL_INLINE T& operator[](I32 i)               { assert(i >= 0 && i < M);  return Data_[i][0]; }
  MTL_INLINE const T& operator[](I32 i) const   { assert(i >= 0 && i < M);  return Data_[i][0]; }

  MTL_INLINE I32 Size() const  { return Rows();}
};

template<I32 N, class T = F64>
class RowVector : public Matrix<1,N,T>
{
  typedef Matrix<1,N,T> Base;

public:
  MTL_MATRIX_COMMON_DEFINITIONS(RowVector, Base, 1, N, T);

  MTL_INLINE RowVector() : Base() {}

  MTL_INLINE explicit RowVector(const T data[N])
  {
    memcpy(Data_, data, DataSizeInBytes());
  }

  MTL_INLINE RowVector& operator*=(const RowVector& v)
  {
    Multiply<N>(Data()[0], v.Data()[0]);
    return *this;
  }
  MTL_INLINE RowVector& operator*=(const Base& v)
  {
    return operator*=(RowVector(v));
  }
  MTL_INLINE RowVector operator*(const Base& v2) const
  {
    RowVector<N,T> v3 = *this;
    v3 *= v2;
    return v3;
  }

  MTL_INLINE RowVector& operator/=(const RowVector& v)
  {
    Divide<N>(Data()[0], v.Data()[0]);
    return *this;
  }
  MTL_INLINE RowVector& operator/=(const Base& v)
  {
    return operator/=(RowVector(v));
  }
  MTL_INLINE RowVector<N,T> operator/(const Base& v2) const
  {
    RowVector<N,T> v3 = *this;
    v3 /= v2;
    return v3;
  }

  MTL_INLINE T Dot(const Base& v) const
  {
    return DotProduct<N>(Data()[0], v.Data()[0]);
  }

  MTL_INLINE T MaxNorm() const
  {
    return Max_Norm<N>(Data()[0]);
  }

  MTL_INLINE void Normalize()
  {
    *this /= Sqrt(SumOfSquares());
  }

  MTL_INLINE T& operator[](I32 i)               { assert(i >= 0 && i < N);  return Data_[0][i]; }
  MTL_INLINE const T& operator[](I32 i) const   { assert(i >= 0 && i < N);  return Data_[0][i]; }

  MTL_INLINE I32 Size() const  { return Cols();}
};

typedef ColumnVector<1,F64> ColumnVector1D;
typedef ColumnVector<2,F64> ColumnVector2D;
typedef ColumnVector<3,F64> ColumnVector3D;
typedef ColumnVector<4,F64> ColumnVector4D;
typedef ColumnVector<5,F64> ColumnVector5D;
typedef ColumnVector<6,F64> ColumnVector6D;

typedef RowVector<1,F64> RowVector1D;
typedef RowVector<2,F64> RowVector2D;
typedef RowVector<3,F64> RowVector3D;
typedef RowVector<4,F64> RowVector4D;
typedef RowVector<5,F64> RowVector5D;
typedef RowVector<6,F64> RowVector6D;

template <class T>
MTL_INLINE static ColumnVector<3,T> Cross(const ColumnVector<3,T>& u, const ColumnVector<3,T>& v)
{
  ColumnVector<3,T> crossProduct;
  crossProduct[0] = u[1] * v[2] - u[2] * v[1];
  crossProduct[1] = u[2] * v[0] - u[0] * v[2];
  crossProduct[2] = u[0] * v[1] - u[1] * v[0];
  return crossProduct; 
}
template <class T>
MTL_INLINE static RowVector<3,T> Cross(const RowVector<3,T>& u, const RowVector<3,T>& v)
{
  RowVector<3,T> crossProduct;
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
  Array2D<3,3,F64,3>::MultiplyRowByMatrix(product[0], Data_[0], B.Data());
  Array2D<3,3,F64,3>::MultiplyRowByMatrix(product[1], Data_[1], B.Data());
  Array2D<3,3,F64,3>::MultiplyRowByMatrix(product[2], Data_[2], B.Data());

  return product;
}

MTL_INLINE ColumnVector3D SquareMatrix3x3::operator*(const ColumnVector3D& B) const
{
  ColumnVector3D product;
  product[0] = Array_2D<1,F64,3>::MultiplyRowByCol(Data_[0], B.Data(), 0);
  product[1] = Array_2D<1,F64,3>::MultiplyRowByCol(Data_[1], B.Data(), 0);
  product[2] = Array_2D<1,F64,3>::MultiplyRowByCol(Data_[2], B.Data(), 0);

  return product;
}

MTL_INLINE SquareMatrix2x2 SquareMatrix2x2::operator*(const SquareMatrix2x2& B) const
{
  SquareMatrix2x2 product;
  Array2D<2,2,F64,2>::MultiplyRowByMatrix(product[0], Data_[0], B.Data());
  Array2D<2,2,F64,2>::MultiplyRowByMatrix(product[1], Data_[1], B.Data());

  return product;
}

MTL_INLINE ColumnVector2D SquareMatrix2x2::operator*(const ColumnVector2D& B) const
{
  ColumnVector2D product;
  product[0] = Array_2D<1,F64,2>::MultiplyRowByCol(Data_[0], B.Data(), 0);
  product[1] = Array_2D<1,F64,2>::MultiplyRowByCol(Data_[1], B.Data(), 0);

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
