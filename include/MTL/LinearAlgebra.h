//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_LINEAR_ALGEBRA_H
#define MTL_LINEAR_ALGEBRA_H

#include "SVD.h"

namespace MTL
{

// Returns rows of eigen vectors that map to the null space.
// Assumes the input matrix is symmetric.
// Return true if decomposition of the matrix converged.
template<class T>
static bool NullSpaceSymmetric(DynamicMatrix<T>& null, I32 maxIterations = 20)
{
  assert(null.Rows() == null.Cols());

  DynamicVector<T> D;
  bool converged = JacobiSVDTransposed(null, D, maxIterations);

  T tolerance = null.Rows() * D[0] * Epsilon<T>();
  I32 rank = (I32)D.Size();
  for (; rank > 0 && D[rank-1] < tolerance; rank--);

  I32 numberOfRows = null.Rows() - rank;
  for (I32 i = 0; i < numberOfRows; i++)
    OptimizedCopy_Sequential(null[i], null[i + rank], null.Cols());

  null.Resize(numberOfRows, null.Cols());

  return converged;
}
// This is the more stable version in case A is singular or nearly singular.
template<class T>
static bool NullSpaceSymmetric(DynamicMatrix<T>& null, const DynamicMatrix<T>& A, I32 maxIterations = 20)
{
  assert(A.Rows() == A.Cols());

  DynamicMatrix<T> U = A;
  DynamicVector<T> D;
  bool converged = JacobiSVDTransposed(U, D, null, maxIterations);

  T tolerance = null.Rows() * D[0] * Epsilon<T>();
  I32 rank = (I32)D.Size();
  for (; rank > 0 && D[rank-1] < tolerance; rank--);

  I32 numberOfRows = null.Rows() - rank;
  for (I32 i = 0; i < numberOfRows; i++)
    OptimizedCopy_Sequential(null[i], null[i + rank], null.Cols());

  null.Resize(numberOfRows, null.Cols());

  return converged;
}

// Compute null space transform of a row vector.
template<I32 N, class T> static Matrix<N,N-1,T> NullSpace(const RowVector<N,T>& v, I32 maxIterations = 20)
{
  T S[N];
  SquareMatrix<N,T> V;
  RowVector<N,T> temp(v);
  JacobiSVD<1,N>(temp, S, V, maxIterations);

  return V.template SubMatrix<0,1,N,(N-1)>();
}

}  // namespace MTL

#endif // MTL_LINEAR_ALGEBRA_H
