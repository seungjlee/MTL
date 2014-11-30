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


#ifndef MTL_SVD_H
#define MTL_SVD_H

#include "Givens.h"

namespace MTL
{

//
// Computes Singular Value Decomposition so that A = U * S * Vt where U and V are orthonormal
// matrices and S is a diagonal matrix with ordered singular values.
//
template<I32 M, I32 N, class T>
static bool JacobiSVD(Matrix<M,N,T>& A, T W[N], SquareMatrix<N,T>& V)
{
  enum { kMaxIterations = M > 30 ? M : 30 };

  T epsilon = Epsilon<T>();
  T epsilon2 = Square(epsilon);

  V.Identity();

  for (I32 i = 0; i < N; i++)
    W[i] = A.ColumnSumOfSquares(i);

  I32 iteration;
  for (iteration = 0; iteration < kMaxIterations; iteration++)
  {
    bool changed = false;

    for (I32 i = 0; i < N-1; i++)
    {
      for (I32 j = i+1; j < N; j++)
      {
        T a = W[i];
        T b = W[j];

        T p = A.DotProductOfColumns(i, j);
        T pp = p*p;

        if (pp <= epsilon2*a*b)
          continue;

        T beta = a - b;
        T gamma = Hypotenuse(T(2.0)*p, beta);
        T delta;
        T c, s;
        if (beta < 0)
        {
          delta = (gamma - beta) * T(0.5);
          s = Sqrt(delta/gamma);
          c = p / (gamma*s);
        }
        else
        {
          T meanGammaBeta = T(0.5) * (gamma + beta);
          delta = pp / meanGammaBeta;
          c = Sqrt(meanGammaBeta / gamma);
          s = p / (gamma*c);
        }

        if (delta <= 0)
          continue;

        changed = true;

        if (iteration < 10)
        {
          GivensRotation<M,N>(A, i, j, c, s);
          W[i] += delta;
          W[j] -= delta;
        }
        else
        {
          GivensRotation<M,N>(A, i, j, c, s, W[i], W[j]);
          W[i] += delta * T(0.5);
          W[j] -= delta * T(0.5);
        }

        GivensRotation<N,N>(V, i, j, c, s);
      }
    }

    if (!changed)
      break;
  }

  for (I32 i = 0; i < N; i++)
    W[i] = Sqrt(A.ColumnSumOfSquares(i));

  // Sort with respect to singular values.
  for (I32 i = 0; i < N-1; i++)
  {
    I32 maxIndex = i;
    for (I32 k = i+1; k < N; k++)
    {
      if (W[maxIndex] < W[k])
        maxIndex = k;
    }

    if (i != maxIndex)
    {
      Swap(W[i], W[maxIndex]);
      for (I32 k = 0; k < M; k++)
        Swap(A[k][i], A[k][maxIndex]);
      for (I32 k = 0; k < N; k++)
        Swap(V[k][i], V[k][maxIndex]);
    }
  }

  for (I32 i = 0; i < N; i++)
  {
    if (W[i] > 0)
      A.ColumnMultiply(i, T(1)/W[i]);
  }

  return iteration < kMaxIterations;
}

template <I32 N, class T>
I32 ComputeRankFromSingularValues(const T D[N], const T& tolerance = Epsilon<T>())
{
  I32 rank = N;
  for (; rank > 0 && D[rank-1] < tolerance; rank--);
  return rank;
}

template <I32 M, I32 N, class T>
MTL_INLINE static T SolveSVD(ColumnVector<N,T>& x, const Matrix<M,N,T>& U, const T D[N],
                             const SquareMatrix<N,T>& V, I32& rank,
                             const ColumnVector<M,T>& b, const T& tol = T(-1.0))
{
  T tolerance = tol;
  if (tolerance < 0)
    tolerance = M * Epsilon<T>() * D[0];

  rank = ComputeRankFromSingularValues<N,T>(D, tolerance);

  assert(rank > 0);

  T Ub[N];
  for (I32 i = 0; i < rank; i++)
  {
    Ub[i] = U[0][i] * b[0];
    for (I32 k = 1; k < M; k++)
      Ub[i] += U[k][i] * b[k];

    Ub[i] /= D[i];
  }

  x.Zeros();
  for (I32 i = 0; i < N; i++)
  {
    for (I32 k = 0; k < rank; k++)
      x[i] += V[i][k] * Ub[k];
  }

  return D[0] / D[rank-1];
}

template <I32 M, I32 N, class T>
MTL_INLINE static bool SolveJacobiSVD(Matrix<M,N,T>& A, ColumnVector<N,T>& x,
                                      I32& rank, T& conditionNumber,
                                      const ColumnVector<M,T>& b, const T& tolerance = T(-1.0))
{
  SquareMatrix<N,T> V;
  T D[N];

  bool fullyConverged = JacobiSVD<M,N>(A, D, V);
  conditionNumber = SolveSVD<M,N>(x, A, D, V, rank, b, tolerance);

  return fullyConverged;
}
template <I32 N, class T>
MTL_INLINE static bool SolveJacobiSVD(SquareMatrix<N,T>& A, ColumnVector<N,T>& x,
                                      I32& rank, T& conditionNumber,
                                      const T& tol = T(-1.0))
{
  return SolveJacobiSVD(A, x, rank, conditionNumber, x, tol);
}

template<I32 M, I32 N, class T>
static Matrix<N,M,T> ComputePseudoinverseJacobiSVD(const Matrix<M,N,T>& A,
                                                   const T& tol = T(-1.0))
{
  if (N > M)
  {
    Matrix<M,N,T> pinv = ComputePseudoinverseJacobiSVD(A.ComputeTranspose(), tol);
    return pinv.ComputeTranspose();
  }

  Matrix<M,N,T> U = A;
  SquareMatrix<N,T> V;
  T D[N];
  JacobiSVD<M,N,T>(U, D, V);

  T tolerance = tol;
  if (tolerance < 0)
    tolerance = A.Rows() * Epsilon<T>() * D[0];

  I32 rank = ComputeRankFromSingularValues<N,T>(D, tolerance);

  Matrix<N,M,T> pinv;
  for (I32 n = 0; n < N; n++)
  {
    for (I32 m = 0; m < M; m++)
    {
      pinv[n][m] = U[m][0] * V[n][0] / D[0];
      for (I32 i = 1; i < rank; i++)
        pinv[n][m] += U[m][i] * V[n][i] / D[i];
    }
  }

  return pinv;
}
}  // namespace MTL

#endif // MTL_SVD_H
