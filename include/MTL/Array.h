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

#ifndef MTL_ARRAY_H
#define MTL_ARRAY_H

#include <MTL/Math.h>

namespace MTL
{

template <int N, class T> class Array
{
public:
  MTL_INLINE static T Sum(const T* p)  { return Array<N-1,T>::Sum(p) + p[N-1];             }
  MTL_INLINE static T Min(const T* p)  { return MTL::Min<T>(Array<N-1,T>::Min(p), p[N-1]); }
  MTL_INLINE static T Max(const T* p)  { return MTL::Max<T>(Array<N-1,T>::Max(p), p[N-1]); }

  MTL_INLINE static T MinOfAbsolutes(const T* p)
  {
    return MTL::Min(Array<N-1,T>::MinOfAbsolutes(p), Abs(p[N-1]));
  }
  MTL_INLINE static T MaxOfAbsolutes(const T* p)
  {
    return MTL::Max(Array<N-1,T>::MaxOfAbsolutes(p), Abs(p[N-1]));
  }

  MTL_INLINE static T SumOfSquares(const T* p)
  {
    return Array<N-1,T>::SumOfSquares(p) + Square(p[N-1]);
  }
  MTL_INLINE static T Dot(const T* p1, const T* p2)
  {
    return Array<N-1,T>::Dot(p1, p2) + p1[N-1] * p2[N-1];
  }

  MTL_INLINE static void UnaryMinus(T* a, const T* b)
  {
    Array<N-1,T>::UnaryMinus(a, b);
    a[N-1] = -b[N-1];
  }

  MTL_INLINE static void Add(T* a, const T* b)
  {
    Array<N-1,T>::Add(a, b);
    a[N-1] += b[N-1];
  }
  MTL_INLINE static void Subtract(T* a, const T* b)
  {
    Array<N-1,T>::Subtract(a, b);
    a[N-1] -= b[N-1];
  }
  MTL_INLINE static void Multiply(T* a, const T* b)
  {
    Array<N-1,T>::Multiply(a, b);
    a[N-1] *= b[N-1];
  }
  MTL_INLINE static void Divide(T* a, const T* b)
  {
    Array<N-1,T>::Divide(a, b);
    a[N-1] /= b[N-1];
  }

  MTL_INLINE static void AddScalar(T* a, const T& scalar)
  {
    Array<N-1,T>::AddScalar(a, scalar);
    a[N-1] += scalar;
  }
  MTL_INLINE static void SubtractScalar(T* a, const T& scalar)
  {
    Array<N-1,T>::SubtractScalar(a, scalar);
    a[N-1] -= scalar;
  }
  MTL_INLINE static void MultiplyScalar(T* a, const T& scalar)
  {
    Array<N-1,T>::MultiplyScalar(a, scalar);
    a[N-1] *= scalar;
  }
  MTL_INLINE static void DivideScalar(T* a, const T& scalar)
  {
    Array<N-1,T>::DivideScalar(a, scalar);
    a[N-1] /= scalar;
  }

  MTL_INLINE static T MaxNorm(const T* a)
  {
    return MTL::Max(Array<N-1,T>::MaxNorm(a), Abs(a[N-1]));
  }

  MTL_INLINE static void Set(T* ptr, const T& newVal)
  {
    Array<N-1,T>::Set(ptr, newVal);
    ptr[N-1] = newVal;
  }

  // Copies or converts.
  template<class SrcT>
  MTL_INLINE static void Set(T* dst, const SrcT* src)
  {
    Array<N-1,T>::Set(dst, src);
    dst[N-1] = (T)src[N-1];
  }
};
template <class T> class Array<1,T>
{
public:
  MTL_INLINE static T Sum(const T* p)                           { return p[0];          }
  MTL_INLINE static T Min(const T* p)                           { return p[0];          }
  MTL_INLINE static T Max(const T* p)                           { return p[0];          }
  MTL_INLINE static T MinOfAbsolutes(const T* p)                { return Abs(p[0]);     }
  MTL_INLINE static T MaxOfAbsolutes(const T* p)                { return Abs(p[0]);     }
  MTL_INLINE static T SumOfSquares(const T* p)                  { return Square(p[0]);  }
  MTL_INLINE static T Dot(const T* p1, const T* p2)             { return p1[0] * p2[0]; }
  MTL_INLINE static void UnaryMinus(T* a, const T* b)           { a[0] = -b[0];         }

  MTL_INLINE static void Add(T* a, const T* b)                  { a[0] += b[0];         }
  MTL_INLINE static void Subtract(T* a, const T* b)             { a[0] -= b[0];         }
  MTL_INLINE static void Multiply(T* a, const T* b)             { a[0] *= b[0];         }
  MTL_INLINE static void Divide(T* a, const T* b)               { a[0] /= b[0];         }
  MTL_INLINE static void AddScalar(T* a, const T& scalar)       { a[0] += scalar;       }
  MTL_INLINE static void SubtractScalar(T* a, const T& scalar)  { a[0] -= scalar;       }
  MTL_INLINE static void MultiplyScalar(T* a, const T& scalar)  { a[0] *= scalar;       }
  MTL_INLINE static void DivideScalar(T* a, const T& scalar)    { a[0] /= scalar;       }

  MTL_INLINE static T MaxNorm(const T* a)                       { return Abs(a[0]);     }

  MTL_INLINE static void Set(T* ptr, const T& newVal)           { ptr[0] = newVal;      }

  template<class SrcT>
  MTL_INLINE static void Set(T* dst, const SrcT* src)           { dst[0] = (T)src[0];   }
};

template <int N, class T, int Q> class Array_2D
{
public:
  MTL_INLINE static T MultiplyRowByCol(const T* pRow, const T colData[][N], I32 col)
  {
    return Array_2D<N,T,Q-1>::MultiplyRowByCol(pRow, colData, col) + pRow[Q-1] * colData[Q-1][col];
  }
  MTL_INLINE static void ColumnScalarMultiply(T* a, const T& scalar)
  {
    *a *= scalar;
    Array_2D<N,T,Q-1>::ColumnScalarMultiply(a + N, scalar);
  }
  MTL_INLINE static void ColumnScalarDivide(T* a, const T& scalar)
  {
    *a /= scalar;
    Array_2D<N,T,Q-1>::ColumnScalarDivide(a + N, scalar);
  }

  MTL_INLINE static T ColumnSum(const T* a)
  {
    return *a + Array_2D<N,T,Q-1>::ColumnSumOfSquares(a + N);
  }
  MTL_INLINE static T ColumnSumOfSquares(const T* a)
  {
    return Square(*a) + Array_2D<N,T,Q-1>::ColumnSumOfSquares(a + N);
  }
  MTL_INLINE static T DotProductOfColumns(const T* a, const T* b)
  {
    return *a * *b + Array_2D<N,T,Q-1>::DotProductOfColumns(a + N, b + N);
  }
};
template <int N, class T> class Array_2D<N,T,1>
{
public:
  MTL_INLINE static T MultiplyRowByCol(const T* pRow, const T colData[][N], I32 col)
  {
    return pRow[0] * colData[0][col];
  }
  MTL_INLINE static void ColumnScalarMultiply(T* a, const T& scalar)
  {
    *a *= scalar;
  }
  MTL_INLINE static void ColumnScalarDivide(T* a, const T& scalar)
  {
    *a /= scalar;
  }

  MTL_INLINE static T ColumnSum(const T* a)
  {
    return *a;
  }
  MTL_INLINE static T ColumnSumOfSquares(const T* a)
  {
    return Square(*a);
  }
  MTL_INLINE static T DotProductOfColumns(const T* a, const T* b)
  {
    return *a * *b;
  }
};

template <int M, int N, class T, int Q> class Array2D
{
public:
  MTL_INLINE static void SetDiagonal(T ptr[M][N], const T& newVal)
  {
    Array2D<M,N,T,Q-1>::SetDiagonal(ptr, newVal);
    ptr[Q-1][Q-1] = newVal;
  }
  MTL_INLINE static void SetDiagonal(T ptr[M][N], const T newVals[])
  {
    Array2D<M,N,T,Q-1>::SetDiagonal(ptr, newVals);
    ptr[Q-1][Q-1] = newVals[Q-1];
  }
  MTL_INLINE static void AddToDiagonal(T ptr[M][N], const T& newVal)
  {
    Array2D<M,N,T,Q-1>::AddToDiagonal(ptr, newVal);
    ptr[Q-1][Q-1] += newVal;
  }

  MTL_INLINE static void MultiplyRowByMatrix(T* pDst, const T* pRow, const T colData[][N])
  {
    Array2D<M,N,T,Q-1>::MultiplyRowByMatrix(pDst, pRow, colData);
    pDst[Q-1] = Array_2D<N,T,M>::MultiplyRowByCol(pRow, colData, Q-1);
  }

  MTL_INLINE static void SetColumn(T dst[M][N], const T src[M][1], I32 col)
  {
    Array2D<M,N,T,Q-1>::SetColumn(dst, src, col);
    dst[Q-1][col] = src[Q-1][0];
  }
};
template <int M, int N, class T> class Array2D<M,N,T,1>
{
public:
  MTL_INLINE static void SetDiagonal(T ptr[M][N], const T& newVal)    { ptr[0][0] = newVal;     }
  MTL_INLINE static void SetDiagonal(T ptr[M][N], const T newVals[])  { ptr[0][0] = newVals[0]; }
  MTL_INLINE static void AddToDiagonal(T ptr[M][N], const T& newVal)  { ptr[0][0] += newVal;    }

  MTL_INLINE static void MultiplyRowByMatrix(T* pDst, const T* pRow, const T colData[][N])
  {
    pDst[0] = Array_2D<N,T,M>::MultiplyRowByCol(pRow, colData, 0);
  }

  MTL_INLINE static void SetColumn(T dst[M][N], const T src[M][1], I32 col)
  {
    dst[0][col] = src[0][0];
  }
};

// Sum & mean of fixed sized array.
template <int N, class T> MTL_INLINE static T Sum(const T* p)
{
  return Array<N,T>::Sum(p);
}
template <int N, class T> MTL_INLINE static T Mean(const T* p)
{
  return Sum<N>(p) / T(N);
}

// Minimum and maximum.
template <int N, class T> MTL_INLINE static T Minimum(const T* p)
{
  return Array<N,T>::Min(p);
}
template <int N, class T> MTL_INLINE static T Maximum(const T* p)
{
  return Array<N,T>::Max(p);
}
template <int N, class T> MTL_INLINE static T MinimumOfAbsolutes(const T* p)
{
  return Array<N,T>::MinOfAbsolutes(p);
}
template <int N, class T> MTL_INLINE static T MaximumOfAbsolutes(const T* p)
{
  return Array<N,T>::MaxOfAbsolutes(p);
}

template <int N, class T> MTL_INLINE static T SumOfSquares(const T* p)
{
  return Array<N,T>::SumOfSquares(p);
}

// Dot product.
template <int N, class T> MTL_INLINE static T Dot(const T* p1, const T* p2)
{
  return Array<N,T>::Dot(p1, p2);
}

}  // namespace MTL

#endif  // MTL_ARRAY_H
