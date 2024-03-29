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


#ifndef MTL_LDLT_H
#define MTL_LDLT_H

#include "DynamicMatrix.h"

namespace MTL
{

// LDLt decomposition of matrix A so that A = L * D * Lt
// Matrix A is passed in as matrix L.
template<I32 N, class T>
static void LDLt_JustForSolver(SquareMatrix<N,T>& L, T D[])
{
  D[0] = L[0][0];

  for (I32 row = 1; row < N; row++)
  {
    for (I32 col = 0; col < row; col++)
    {
      for (I32 k = 0; k < col; k++)
        L[row][col] -= L[row][k] * L[col][k] * D[k];

      L[row][col] /= D[col];
    }

    D[row] = L[row][row];

    for (I32 col = 0; col < row; col++)
      D[row] -= Pow<2>(L[row][col]) * D[col];
  }
}
template<I32 N, class T>
static void LDLt(SquareMatrix<N,T>& L, T D[])
{
  LDLt_JustForSolver(L, D);

  // Fill ones and zeros in L.
  for (I32 row = 0; row < N; row++)
  {
    L[row][row] = T(1);

    for (I32 col = row + 1; col < N; col++)
      L[row][col] = T(0);
  }
}

template<I32 N, class T>
static void SolveLDLt(ColumnVector<N,T>& x, const SquareMatrix<N,T>& L, const T D[])
{
  // Compute Inv(L) * x using forward substitution.
  for (I32 row = 1; row < N; row++)
  {
    for (I32 col = 0; col < row; col++)
      x[row] -= L[row][col] * x[col];
  }

  // Compute Inv(D) * x.
  Array<N,T>::Divide(&x[0], D);

  // Finally, compute Inv(Lt) * x using back substitution.
  for (I32 row = N-2; row >= 0; row--)
  {
    for (I32 col = row + 1; col < N; col++)
      x[row] -= L[col][row] * x[col];
  }
}

//
// Solves A * x = b using LDLt Decomposition. b is passed as x.
// Returns the effective rank of matrix A.
//
// Note that A will return an incomplete and dirty L.
//
template<I32 N, class T>
static I32 SolveLDLt(ColumnVector<N,T>& x, SquareMatrix<N,T>& A,
                     const T& tolerance = NumericalEpsilon<T>())
{
  T D[N];
  LDLt_JustForSolver(A, D);

  I32 rank = 0;
  for (I32 i = 0; i < N; i++)
    if (Abs(D[i]) > tolerance)
      rank++;

  if (rank != N)
  {
    // Rank deficient.
    x.Zeros();
    return rank;
  }

  SolveLDLt(x, A, D);

  return rank;
}

template<class T>
static void LDLt_JustForSolver(T* L, T* D, I32 N, I32 rowSize)
{
  D[0] = L[0];

  for (I32 row = 1; row < N; row++)
  {
    for (I32 col = 0; col < row; col++)
    {
      for (I32 k = 0; k < col; k++)
        L[row * rowSize + col] -= L[row * rowSize + k] * L[col * rowSize + k] * D[k];

      L[row * rowSize + col] /= D[col];
    }

    D[row] = L[row * rowSize + row];

    for (I32 col = 0; col < row; col++)
      D[row] -= Pow<2>(L[row * rowSize + col]) * D[col];
  }
}

template<class T>
static void SolveLDLt(T* x, T* L, const T* D, I32 N, I32 rowSize)
{
  // Compute Inv(L) * x using forward substitution.
  for (I32 row = 1; row < N; row++)
  {
    for (I32 col = 0; col < row; col++)
      x[row] -= L[row * rowSize + col] * x[col];
  }

  // Compute Inv(D) * x.
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  Division_StreamAligned_Sequential(x, D, N);
#else
  Division_Sequential(x, D, N);
#endif

  // Finally, compute Inv(Lt) * x using back substitution.
  for (I32 row = N-2; row >= 0; row--)
  {
    for (I32 col = row + 1; col < N; col++)
      x[row] -= L[col * rowSize + row] * x[col];
  }
}

template<class T>
static I32 SolveLDLt(T* x, T* A, I32 N, I32 rowSize,
                     const T& tolerance = NumericalEpsilon<T>())
{
  DynamicVector<T> D(N);
  LDLt_JustForSolver(A, D.Begin(), N, rowSize);

  I32 rank = 0;
  for (I32 i = 0; i < N; i++)
    if (Abs(D[i]) > tolerance)
      rank++;

  if (rank != N)
  {
    // Rank deficient.
    OptimizedZeros(x, N);
    return rank;
  }

  SolveLDLt(x, A, D.Begin(), N, rowSize);

  return rank;
}

//
// Solves A * x = b using LDLt Decomposition. b is passed as x.
// Returns the effective rank of matrix A.
//
// Note that A will return an incomplete and dirty L.
//
template<class T>
static I32 SolveLDLt(DynamicVector<T>& x, DynamicMatrix<T>& A,
                     const T& tolerance = NumericalEpsilon<T>())
{
  I32 N = A.Rows();

  // Maybe we want to also check here if A is symmetric.
  assert(N == A.Cols());
  assert(N == (I32)x.Size());

  return SolveLDLt(x.Begin(), A[0], N, A.RowSize(), tolerance);
}

}  // namespace MTL

#endif // MTL_LDLT_H
