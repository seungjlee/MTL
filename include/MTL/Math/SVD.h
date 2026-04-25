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


#ifndef MTL_SVD_H
#define MTL_SVD_H

#include "DynamicMatrix.h"
#include "Givens.h"
#include "Point2D.h"

#ifndef SVD_MAX_ITERATIONS
#define SVD_MAX_ITERATIONS 20
#endif

namespace MTL
{

//
// Computes Singular Value Decomposition so that A = U * S * Vt where U and V are orthonormal
// matrices and S is a diagonal matrix with ordered singular values.
//
template<I32 M, I32 N, class T>
static bool JacobiSVD(Matrix<M,N,T>& A, T W[N], SquareMatrix<N,T>& V)
{
  const T epsilon = NumericalEpsilon<T>();
  const T epsilon2 = Square(epsilon);

  V.Identity();

  for (I32 i = 0; i < N; i++)
    W[i] = A.ColumnSumOfSquares(i);

  I32 iteration;
  for (iteration = 0; iteration < SVD_MAX_ITERATIONS; iteration++)
  {
    T sumOfDeltasSquared = 0;

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

        sumOfDeltasSquared += Square(delta);

        GivensRotation<M,N>(A, i, j, c, s);
        W[i] += delta;
        W[j] -= delta;

        GivensRotation<N,N>(V, i, j, c, s);
      }
    }

#ifdef DEBUG_SVD_CONVERGENCE
    wprintf(L"sumOfDeltasSquared = %.12e\n", sumOfDeltasSquared);
#endif
    if (sumOfDeltasSquared < epsilon2)
      break;
  }

#ifdef DEBUG_SVD_CONVERGENCE
  ConsoleOut << L"(1) Iterations: " << iteration << L", Max: " << SVD_MAX_ITERATIONS << std::endl;
#endif

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

  return iteration < SVD_MAX_ITERATIONS;
}

template <I32 N, class T>
I32 ComputeRankFromSingularValues(const T D[N], const T& tolerance = NumericalEpsilon<T>())
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

  bool converged = JacobiSVD<M,N>(A, D, V);
  conditionNumber = SolveSVD<M,N>(x, A, D, V, rank, b, tolerance);

  return converged;
}
template <I32 N, class T>
MTL_INLINE static bool SolveJacobiSVD(SquareMatrix<N,T>& A, ColumnVector<N,T>& x,
                                      I32& rank, T& conditionNumber,
                                      const T& tol = T(-1.0))
{
  return SolveJacobiSVD(A, x, rank, conditionNumber, x, tol);
}

template <I32 M, I32 N, class T>
MTL_INLINE static bool SolveJacobiSVDHomogeneous(Matrix<M,N,T>& A, ColumnVector<N,T>& x,
                                                 I32& rank, T& conditionNumber,
                                                 const T& tol = T(-1.0))
{
  SquareMatrix<N,T> V;
  T D[N];

  bool converged = JacobiSVD<M,N>(A, D, V);

  T tolerance = tol;
  if (tolerance < 0)
    tolerance = A.Rows() * Epsilon<T>() * D[0];

  rank = ComputeRankFromSingularValues<N,T>(D, tolerance);

  assert(rank > 0);

  for (I32 i = 0; i < N; i++)
    x[i] = V[i][N-1];

  conditionNumber = D[0] / D[rank-1];

  return converged;
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

//
// Computes Singular Value Decomposition by first reducing A to upper bidiagonal form using
// Householder reflections, then applying one-sided Jacobi rotations on the (much smaller)
// bidiagonal matrix. Faster than plain JacobiSVD when M >> N because the Jacobi sweeps then
// operate on an N x N matrix instead of an M x N matrix.
//
// Output: A holds U (M x N orthonormal columns), W holds singular values, V is right factor.
// Requires M >= N.
//
template<I32 M, I32 N, class T>
static bool HouseholderJacobiSVD(Matrix<M,N,T>& A, T W[N], SquareMatrix<N,T>& V)
{
  static_assert(M >= N, "HouseholderJacobiSVD requires M >= N");

  // Zero-initialized so that the trailing super[N-1] (never written, never read by the
  // bidiagonal builder below) does not trigger a -Wmaybe-uninitialized false positive that
  // GCC raises through the inlined Array_2D::ColumnSumOfSquares template. Cost is N stores
  // on small stack arrays; SVD itself dominates by orders of magnitude.
  T diag[N] = {};
  T super[N] = {};

  // Bidiagonalize A in place. Reflectors are stored in the lower triangle (left) and the
  // strict upper triangle right of the superdiagonal (right). All reflectors are normalized so
  // that the reflection is H = I - 2 v v^T.
  for (I32 k = 0; k < N; k++)
  {
    // Left Householder zeroing column k below diagonal.
    T sigma = T(0);
    for (I32 i = k; i < M; i++)
      sigma += Square(A[i][k]);

    T alpha = Sqrt(sigma);
    if (A[k][k] > T(0))
      alpha = -alpha;

    T v0 = A[k][k] - alpha;
    T vnorm2 = sigma - Square(A[k][k]) + Square(v0);
    diag[k] = alpha;

    if (vnorm2 > T(0))
    {
      T invNorm = T(1) / Sqrt(vnorm2);
      A[k][k] = v0 * invNorm;
      for (I32 i = k+1; i < M; i++)
        A[i][k] *= invNorm;

      for (I32 j = k+1; j < N; j++)
      {
        T dot = T(0);
        for (I32 i = k; i < M; i++)
          dot += A[i][k] * A[i][j];
        dot *= T(2);
        for (I32 i = k; i < M; i++)
          A[i][j] -= dot * A[i][k];
      }
    }
    else
    {
      A[k][k] = T(0);
    }

    // Right Householder zeroing row k right of superdiagonal.
    if (k + 2 < N)
    {
      T sigmaR = T(0);
      for (I32 j = k+1; j < N; j++)
        sigmaR += Square(A[k][j]);

      T alphaR = Sqrt(sigmaR);
      if (A[k][k+1] > T(0))
        alphaR = -alphaR;

      T u0 = A[k][k+1] - alphaR;
      T unorm2 = sigmaR - Square(A[k][k+1]) + Square(u0);
      super[k] = alphaR;

      if (unorm2 > T(0))
      {
        T invNorm = T(1) / Sqrt(unorm2);
        A[k][k+1] = u0 * invNorm;
        for (I32 j = k+2; j < N; j++)
          A[k][j] *= invNorm;

        for (I32 i = k+1; i < M; i++)
        {
          T dot = T(0);
          for (I32 j = k+1; j < N; j++)
            dot += A[k][j] * A[i][j];
          dot *= T(2);
          for (I32 j = k+1; j < N; j++)
            A[i][j] -= dot * A[k][j];
        }
      }
      else
      {
        A[k][k+1] = T(0);
      }
    }
    else if (k + 1 < N)
    {
      super[k] = A[k][k+1];
    }
  }

  // Build explicit bidiagonal matrix B (N x N).
  SquareMatrix<N,T> B;
  B.Zeros();
  for (I32 k = 0; k < N; k++)
  {
    B[k][k] = diag[k];
    if (k + 1 < N)
      B[k][k+1] = super[k];
  }

  // Apply Jacobi SVD to the bidiagonal matrix.
  bool converged = JacobiSVD<N,N,T>(B, W, V);
  // After JacobiSVD: B's columns are U_b (orthonormal), V is V_b.

  // Form U = U_h * U_b. Start with U holding U_b in the top N rows (as an M x N matrix),
  // then apply each stored left reflector L_k from k = N-1 down to 0.
  Matrix<M,N,T> U;
  U.Zeros();
  for (I32 i = 0; i < N; i++)
    for (I32 j = 0; j < N; j++)
      U[i][j] = B[i][j];

  for (I32 k = N-1; k >= 0; k--)
  {
    // Reflector v stored in column k of A from row k to M-1.
    bool hasReflector = (A[k][k] != T(0));
    if (!hasReflector)
    {
      for (I32 i = k+1; i < M && !hasReflector; i++)
        if (A[i][k] != T(0))
          hasReflector = true;
      if (!hasReflector)
        continue;
    }

    for (I32 j = 0; j < N; j++)
    {
      T dot = T(0);
      for (I32 i = k; i < M; i++)
        dot += A[i][k] * U[i][j];
      dot *= T(2);
      for (I32 i = k; i < M; i++)
        U[i][j] -= dot * A[i][k];
    }
  }

  // Form V = V_h * V_b. V already holds V_b; apply right reflectors R_k from k = N-3 down to 0.
  for (I32 k = N-3; k >= 0; k--)
  {
    bool hasReflector = (A[k][k+1] != T(0));
    if (!hasReflector)
    {
      for (I32 j = k+2; j < N && !hasReflector; j++)
        if (A[k][j] != T(0))
          hasReflector = true;
      if (!hasReflector)
        continue;
    }

    for (I32 j = 0; j < N; j++)
    {
      T dot = T(0);
      for (I32 i = k+1; i < N; i++)
        dot += A[k][i] * V[i][j];
      dot *= T(2);
      for (I32 i = k+1; i < N; i++)
        V[i][j] -= dot * A[k][i];
    }
  }

  // Copy U back into A.
  for (I32 i = 0; i < M; i++)
    for (I32 j = 0; j < N; j++)
      A[i][j] = U[i][j];

  return converged;
}

template <I32 M, I32 N, class T>
MTL_INLINE static bool SolveHouseholderJacobiSVD(Matrix<M,N,T>& A, ColumnVector<N,T>& x,
                                                 I32& rank, T& conditionNumber,
                                                 const ColumnVector<M,T>& b,
                                                 const T& tolerance = T(-1.0))
{
  SquareMatrix<N,T> V;
  T D[N];

  bool converged = HouseholderJacobiSVD<M,N>(A, D, V);
  conditionNumber = SolveSVD<M,N>(x, A, D, V, rank, b, tolerance);

  return converged;
}
template <I32 N, class T>
MTL_INLINE static bool SolveHouseholderJacobiSVD(SquareMatrix<N,T>& A, ColumnVector<N,T>& x,
                                                 I32& rank, T& conditionNumber,
                                                 const T& tol = T(-1.0))
{
  return SolveHouseholderJacobiSVD(A, x, rank, conditionNumber, x, tol);
}

template <I32 M, I32 N, class T>
MTL_INLINE static bool SolveHouseholderJacobiSVDHomogeneous(Matrix<M,N,T>& A,
                                                            ColumnVector<N,T>& x,
                                                            I32& rank, T& conditionNumber,
                                                            const T& tol = T(-1.0))
{
  SquareMatrix<N,T> V;
  T D[N];

  bool converged = HouseholderJacobiSVD<M,N>(A, D, V);

  T tolerance = tol;
  if (tolerance < 0)
    tolerance = A.Rows() * Epsilon<T>() * D[0];

  rank = ComputeRankFromSingularValues<N,T>(D, tolerance);

  assert(rank > 0);

  for (I32 i = 0; i < N; i++)
    x[i] = V[i][N-1];

  conditionNumber = D[0] / D[rank-1];

  return converged;
}

template<I32 M, I32 N, class T>
static Matrix<N,M,T> ComputePseudoinverseHouseholderJacobiSVD(const Matrix<M,N,T>& A,
                                                              const T& tol = T(-1.0))
{
  if constexpr (N > M)
  {
    Matrix<M,N,T> pinv = ComputePseudoinverseHouseholderJacobiSVD(A.ComputeTranspose(), tol);
    return pinv.ComputeTranspose();
  }
  else
  {
    Matrix<M,N,T> U = A;
    SquareMatrix<N,T> V;
    T D[N];
    HouseholderJacobiSVD<M,N,T>(U, D, V);

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
}

template<class T>
static T JacobiRotationsTransposed(T& c, T& s, I32 i, I32 j, T* At, T* W, I32 M, I32 rowSizeA)
{
  const T epsilon = NumericalEpsilon<T>();
  const T epsilon2 = Square(epsilon);

  T* Ai = At + i*rowSizeA;
  T* Aj = At + j*rowSizeA;
  T a = W[i];
  T b = W[j];
        
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  T p = DotProduct_StreamAligned_Sequential(Ai, Aj, M);
#else
  T p = DotProduct_Sequential(Ai, Aj, M);
#endif
  T pp = p*p;

  if (pp <= epsilon2*a*b)
    return T(0);

  T beta = a - b;
  T gamma = Sqrt(T(4)*pp + beta*beta);
  T delta;
  if (beta < 0)
  {
    delta = (gamma - beta) * 0.5;
    s = Sqrt(delta/gamma);
    c = p / (gamma*s);
  }
  else
  {
    T meanGammaBeta = 0.5 * (gamma + beta);
    delta = pp / meanGammaBeta;
    c = Sqrt(meanGammaBeta / gamma);
    s = p / (gamma*c);
  }

  if (delta <= 0)
    return T(0);

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  GivensRotation_StreamAligned_Sequential(Ai, Aj, c, s, M);
#else
  GivensRotation_Sequential(Ai, Aj, c, s, M);
#endif
  W[i] += delta;
  W[j] -= delta;

  return delta;
}

static void ComputeJacobiParallelPairs(DynamicVector<DynamicVector<Point2D<I32>>>& pairs, I32 N)
{
  DynamicVector<Point2D<I32>> allPairs;
  allPairs.Reserve((N-1)*N/2);
  for (I32 i = 0; i < N-1; i++)
  {
    for (I32 j = i+1; j < N; j++)
    {
      allPairs.PushBack(Point2D<I32>(i,j));
    }
  }

  DynamicVector<U8> used(N);

  U32 count = 0;
  while (count < allPairs.Size())
  {
    used.Zeros();
    pairs.Resize(pairs.Size() + 1);
    DynamicVector<Point2D<I32>>& set = pairs.Back();

    for (Point2D<I32>& pair : allPairs)
    {
      if (pair.x() >= 0 && !used[pair.x()] && !used[pair.y()])
      {
        set.PushBack(pair);
        used[pair.x()] = 1;
        used[pair.y()] = 1;
        pair.x(-1);
        count++;
      }
    }
  }
}

// Dynamic matrix version of SVD. Note that it takes A transposed. It expects At and Vt to be
// memory aligned.
template<class T>
static bool JacobiSVDTransposedParallel(T* At, T* W, T* Vt, I32 M, I32 N, I32 rowSizeA, I32 rowSizeV)
{
  I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

  OptimizedZeros(Vt, N*rowSizeV);

  for (I32 i = 0; i < N; i++)
  {
    Vt[i*rowSizeV + i] = T(1);
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M);
#else
    W[i] = SumOfSquares_Sequential(At + i*rowSizeA, M);
#endif
  }

  DynamicVector<DynamicVector<Point2D<I32>>> parallelPairs;
  ComputeJacobiParallelPairs(parallelPairs, N);
  DynamicVector<T> deltas;

  I32 iteration;
  for (iteration = 0; iteration < SVD_MAX_ITERATIONS; iteration++)
  {
    T sumOfDeltasSquared = 0;
    
    for (const auto& set : parallelPairs)
    {
      if (set.Size() > 1)
      {
        deltas.Resize(set.Size());
        {
          I32 blockSize = (I32)MTL::ComputeParallelSubSizesBlockSize(set.Size(), numberOfThreads);

          #pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic, blockSize)
          for (I32 k = 0; k < (I32)set.Size(); k++)
          {
            int i = set[k].x();
            int j = set[k].y();
            T c, s;

            deltas[k] = JacobiRotationsTransposed(c, s, i, j, At, W, M, rowSizeA);
            if (deltas[k] > 0)
            {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
              GivensRotation_StreamAligned_Sequential(Vt + i * rowSizeV, Vt + j * rowSizeV, c, s, N);
#else
              GivensRotation_Sequential(Vt + i * rowSizeV, Vt + j * rowSizeV, c, s, N);
#endif
            }
          }
        }

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
        sumOfDeltasSquared += SumOfSquares_StreamAligned_Sequential(deltas.Begin(), deltas.Size());
#else
        sumOfDeltasSquared += SumOfSquares_Sequential(deltas.Begin(), deltas.Size());
#endif
      }
      else
      {
        int i = set[0].x();
        int j = set[0].y();
        T c, s;

        T delta = JacobiRotationsTransposed(c, s, i, j, At, W, M, rowSizeA);
        if (delta > 0)
        {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
          GivensRotation_StreamAligned_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#else
          GivensRotation_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#endif
          sumOfDeltasSquared += Square(delta);
        }
      }
    }

#ifdef DEBUG_SVD_CONVERGENCE
    wprintf(L"sumOfDeltasSquared = %.12e\n", sumOfDeltasSquared);
#endif
    if (sumOfDeltasSquared < NumericalEpsilonSquared<T>())
      break;
  }

#ifdef DEBUG_SVD_CONVERGENCE
  ConsoleOut << L"(2) Iterations: " << iteration << L", Max: " << SVD_MAX_ITERATIONS << std::endl;
#endif

  for (I32 i = 0; i < N; i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = Sqrt(SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M));
#else
    W[i] = Sqrt(SumOfSquares_Sequential(At + i*rowSizeA, M));
#endif

  for (I32 i = 0; i < N; i++)
  {
    if (W[i] > 0)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
      ScalarMultiplication_StreamAligned_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#else
      ScalarMultiplication_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#endif
  }

  return iteration < SVD_MAX_ITERATIONS;
}
template<class T>
static bool JacobiSVDTransposed(T* At, T* W, T* Vt, I32 M, I32 N, I32 rowSizeA, I32 rowSizeV)
{
  OptimizedZeros(Vt, N*rowSizeV);

  for (I32 i = 0; i < N; i++)
  {
    Vt[i*rowSizeV + i] = T(1);
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M);
#else
    W[i] = SumOfSquares_Sequential(At + i*rowSizeA, M);
#endif
  }

  I32 iteration;
  for (iteration = 0; iteration < SVD_MAX_ITERATIONS; iteration++)
  {
    T sumOfDeltasSquared = 0;

    for (I32 i = 0; i < N-1; i++)
    {
      for (I32 j = i+1; j < N; j++)
      {
        T c, s;
        T delta = JacobiRotationsTransposed(c, s, i, j, At, W, M, rowSizeA);
        if (delta > 0)
        {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
          GivensRotation_StreamAligned_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#else
          GivensRotation_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#endif
          sumOfDeltasSquared += Square(delta);
        }
      }
    }

#ifdef DEBUG_SVD_CONVERGENCE
    wprintf(L"sumOfDeltasSquared = %.12e\n", sumOfDeltasSquared);
#endif
    if (sumOfDeltasSquared < Epsilon<T>())
      break;
  }

#ifdef DEBUG_SVD_CONVERGENCE
  ConsoleOut << L"Iterations: " << iteration << L", Max: " << SVD_MAX_ITERATIONS << std::endl;
#endif

  for (I32 i = 0; i < N; i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = Sqrt(SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M));
#else
    W[i] = Sqrt(SumOfSquares_Sequential(At + i*rowSizeA, M));
#endif

  for (I32 i = 0; i < N; i++)
  {
    if (W[i] > 0)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
      ScalarMultiplication_StreamAligned_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#else
      ScalarMultiplication_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#endif
  }

  return iteration < SVD_MAX_ITERATIONS;
}

template<class T>
static bool JacobiSVDTransposedParallel(T* At, T* W, I32 M, I32 N, I32 rowSizeA)
{
  I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();
  
  for (I32 i = 0; i < N; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M);
#else
    W[i] = SumOfSquares_Sequential(At + i*rowSizeA, M);
#endif
  }

  DynamicVector<DynamicVector<Point2D<I32>>> parallelPairs;
  ComputeJacobiParallelPairs(parallelPairs, N);
  DynamicVector<T> deltas;

  I32 iteration;
  for (iteration = 0; iteration < SVD_MAX_ITERATIONS; iteration++)
  {
    T sumOfDeltasSquared = 0;

    for (const auto& set : parallelPairs)
    {
      if (set.Size() > 1)
      {
        deltas.Resize(set.Size());
        {
          I32 blockSize = (I32)MTL::ComputeParallelSubSizesBlockSize(set.Size(), numberOfThreads);

          #pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic, blockSize)
          for (I32 k = 0; k < (I32)set.Size(); k++)
          {
            int i = set[k].x();
            int j = set[k].y();
            T c, s;

            deltas[k] = JacobiRotationsTransposed(c, s, i, j, At, W, M, rowSizeA);
          }
        }

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
        sumOfDeltasSquared += SumOfSquares_StreamAligned_Sequential(deltas.Begin(), deltas.Size());
#else
        sumOfDeltasSquared += SumOfSquares_Sequential(deltas.Begin(), deltas.Size());
#endif
      }
      else
      {
        T c, s;
        T delta = JacobiRotationsTransposed(c, s, set[0].x(), set[0].y(), At, W, M, rowSizeA);
        sumOfDeltasSquared += Square(delta);
      }
    }

#ifdef DEBUG_SVD_CONVERGENCE
    wprintf(L"sumOfDeltasSquared = %.12e\n", sumOfDeltasSquared);
#endif
    if (sumOfDeltasSquared < NumericalEpsilonSquared<T>())
      break;
  }

#ifdef DEBUG_SVD_CONVERGENCE
  ConsoleOut << L"(3) Iterations: " << iteration << L", Max: " << SVD_MAX_ITERATIONS << std::endl;
#endif

  for (I32 i = 0; i < N; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    W[i] = Sqrt(SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M));
#else
    W[i] = Sqrt(SumOfSquares_Sequential(At + i*rowSizeA, M));
#endif

    if (W[i] > 0)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
      ScalarMultiplication_StreamAligned_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#else
      ScalarMultiplication_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#endif
  }

  return iteration < SVD_MAX_ITERATIONS;
}

template<class T>
static void SortSingularValuesAscending(T* At, T* W, T* Vt, I32 M, I32 N,
                                        I32 rowSizeA, I32 rowSizeV)
{
  // Selection sort of singular values.
  for (I32 i = 0; i < N-1; i++)
  {
    I32 minIndex = i;
    for (I32 k = i+1; k < N; k++)
    {
      if (W[minIndex] > W[k])
        minIndex = k;
    }

    if (i != minIndex)
    {
      Swap(W[i], W[minIndex]);
      SwapRows(At, i, minIndex, M, rowSizeA);
      SwapRows(Vt, i, minIndex, N, rowSizeV);
    }
  }
}
template<class T>
static void SortSingularValuesDescending(T* At, T* W, T* Vt, I32 M, I32 N,
                                         I32 rowSizeA, I32 rowSizeV)
{
  // Selection sort of singular values.
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
      SwapRows(At, i, maxIndex, M, rowSizeA);
      SwapRows(Vt, i, maxIndex, N, rowSizeV);
    }
  }
}

template<class T>
static void SortSingularValuesAscending(T* At, T* W, I32 M, I32 N, I32 rowSizeA)
{
  // Selection sort of singular values.
  for (I32 i = 0; i < N-1; i++)
  {
    I32 minIndex = i;
    for (I32 k = i+1; k < N; k++)
    {
      if (W[minIndex] > W[k])
        minIndex = k;
    }

    if (i != minIndex)
    {
      Swap(W[i], W[minIndex]);
      SwapRows(At, i, minIndex, M, rowSizeA);
    }
  }
}
template<class T>
static void SortSingularValuesDescending(T* At, T* W, I32 M, I32 N, I32 rowSizeA)
{
  // Selection sort of singular values.
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
      SwapRows(At, i, maxIndex, M, rowSizeA);
    }
  }
}

template <class T>
static bool JacobiSVDTransposed(DynamicMatrix<T>& Ut,
                                DynamicVector<T>& D,
                                DynamicMatrix<T>& Vt,
                                bool sortAscending = false)
{
  Vt.Resize(Ut.Rows(), Ut.Rows());
  D.Resize(Ut.Rows());
  bool converged = JacobiSVDTransposedParallel(Ut[0], D.Begin(), Vt[0], Ut.Cols(), Ut.Rows(),
                                               Ut.RowSize(), Vt.RowSize());
  if (sortAscending)
  {
    SortSingularValuesAscending(Ut[0], D.Begin(), Vt[0], Ut.Cols(), Ut.Rows(),
                                Ut.RowSize(), Vt.RowSize());
  }
  else
  {
    SortSingularValuesDescending(Ut[0], D.Begin(), Vt[0], Ut.Cols(), Ut.Rows(),
                                 Ut.RowSize(), Vt.RowSize());
  }
  return converged;
}

template <class T>
static bool JacobiSVDTransposed(DynamicMatrix<T>& Ut,
                                DynamicVector<T>& D,
                                bool sortAscending = false)
{
  D.Resize(Ut.Rows());
  bool converged = JacobiSVDTransposedParallel(Ut[0], D.Begin(), Ut.Cols(), Ut.Rows(),
                                               Ut.RowSize());

  if (sortAscending)
  {
    SortSingularValuesAscending(Ut[0], D.Begin(), Ut.Cols(), Ut.Rows(), Ut.RowSize());
  }
  else
  {
    SortSingularValuesDescending(Ut[0], D.Begin(), Ut.Cols(), Ut.Rows(), Ut.RowSize());
  }
  return converged;
}

// Eigen value decomposition for positive definite matrices.
template <class T>
static bool JacobiEigen(DynamicMatrix<T>& U, DynamicVector<T>& D, bool sortAscending = true)
{
  return JacobiSVDTransposed(U, D, sortAscending);
}

template <class T>
MTL_INLINE static T SolveSVDTransposed(DynamicVector<T>& x,
                                       DynamicMatrix<T>& Ut,
                                       DynamicVector<T>& D,
                                       DynamicMatrix<T>& Vt,
                                       I32& rank,
                                       const DynamicVector<T>& b,
                                       const T& tol = T(-1.0))
{
  T tolerance = tol;
  if (tolerance < 0)
    tolerance = Ut.Cols() * Epsilon<T>() * D[0];

  rank = (I32)D.Size();
  for (; rank > 0 && D[rank-1] < tolerance; rank--);

  Ut.Resize(rank, Ut.Cols());
  Vt.Resize(rank, Vt.Cols());
  D.Resize(rank);

  DynamicVector<T> X = Ut * b;
  X /= D;
  x = Vt.ComputeTranspose() * X;

  return D[0] / D[rank-1];
}

template <class T>
MTL_INLINE static T SolveEigen(DynamicVector<T>& x,
                               DynamicMatrix<T>& U,
                               DynamicVector<T>& D,
                               I32& rank,
                               const DynamicVector<T>& b,
                               const T& tol = T(-1.0))
{
  T tolerance = tol;
  if (tolerance < 0)
    tolerance = U.Cols() * Epsilon<T>() * D[0];

  rank = (I32)D.Size();
  for (; rank > 0 && D[rank-1] < tolerance; rank--);

  U.Resize(rank, U.Cols());
  D.Resize(rank);

  DynamicVector<T> X = U * b;
  X /= D;
  x = U.ComputeTranspose() * X;

  return D[0] / D[rank-1];
}

//
// Solves system of linear equations using SVD (Singular Value Decomposition).
// Solves A * x = b. A is an MxN matrix. At is A transposed.
// Returns true if SVD fully converged before reaching the maximum number of iterations.
// Note that even if it returns false, the answer might be good enough.
//
template <class T>
MTL_INLINE static bool SolveJacobiSVDTransposed(DynamicVector<T>& x, I32& rank,
                                                T& conditionNumber,
                                                const DynamicMatrix<T>& At,
                                                const DynamicVector<T>& b,
                                                const T& tolerance = T(-1.0))
{
  assert(At.Cols() == (I32)b.Size());

  DynamicMatrix<T> Ut = At;
  DynamicMatrix<T> Vt;
  DynamicVector<T> D;

  bool fullyConverged = JacobiSVDTransposed(Ut, D, Vt);
  conditionNumber = SolveSVDTransposed(x, Ut, D, Vt, rank, b, tolerance);

  return fullyConverged;
}
template <class T>
MTL_INLINE static bool SolveJacobiEigen(DynamicVector<T>& x, I32& rank,
                                        T& conditionNumber,
                                        const DynamicMatrix<T>& A,
                                        const DynamicVector<T>& b,
                                        const T& tolerance = T(-1.0))
{
  // A should be symmetric. For now, it will work even for non-symmetric matrices since calls SVD.
  assert(A.Cols() == (I32)A.Rows());
  assert(A.Cols() == (I32)b.Size());

  DynamicMatrix<T> U = A;
  DynamicVector<T> D;

  //
  // Should probably implement a more efficient routine for symmetric matrices.
  //
  bool fullyConverged = JacobiSVDTransposed(U, D);
  conditionNumber = SolveEigen(x, U, D, rank, b, tolerance);

  return fullyConverged;
}

//
// Dynamic-matrix version of HouseholderJacobiSVD.
// A is M x N with M >= N. On output A holds U (M x N orthonormal columns), W has the N
// singular values, V is N x N. Performs Householder bidiagonalization first, then applies
// JacobiSVDTransposed on the small N x N bidiagonal block. Faster than plain JacobiSVD when
// M >> N.
//
template <class T>
static bool HouseholderJacobiSVD(DynamicMatrix<T>& A, DynamicVector<T>& W, DynamicMatrix<T>& V)
{
  const I32 M = A.Rows();
  const I32 N = A.Cols();
  assert(M >= N);

  W.Resize(N);
  V.Resize(N, N);

  DynamicVector<T> diag(N);
  DynamicVector<T> super(N);

  // Bidiagonalize A in place. Left reflectors stored in column k of A from row k to M-1.
  // Right reflectors stored in row k of A from column k+1 to N-1.
  for (I32 k = 0; k < N; k++)
  {
    // Left Householder zeroing column k below diagonal.
    T sigma = T(0);
    for (I32 i = k; i < M; i++)
      sigma += Square(A[i][k]);

    T alpha = Sqrt(sigma);
    if (A[k][k] > T(0))
      alpha = -alpha;

    T v0 = A[k][k] - alpha;
    T vnorm2 = sigma - Square(A[k][k]) + Square(v0);
    diag[k] = alpha;

    if (vnorm2 > T(0))
    {
      T invNorm = T(1) / Sqrt(vnorm2);
      A[k][k] = v0 * invNorm;
      for (I32 i = k+1; i < M; i++)
        A[i][k] *= invNorm;

      for (I32 j = k+1; j < N; j++)
      {
        T dot = T(0);
        for (I32 i = k; i < M; i++)
          dot += A[i][k] * A[i][j];
        dot *= T(2);
        for (I32 i = k; i < M; i++)
          A[i][j] -= dot * A[i][k];
      }
    }
    else
    {
      A[k][k] = T(0);
    }

    // Right Householder zeroing row k right of superdiagonal.
    if (k + 2 < N)
    {
      T sigmaR = T(0);
      for (I32 j = k+1; j < N; j++)
        sigmaR += Square(A[k][j]);

      T alphaR = Sqrt(sigmaR);
      if (A[k][k+1] > T(0))
        alphaR = -alphaR;

      T u0 = A[k][k+1] - alphaR;
      T unorm2 = sigmaR - Square(A[k][k+1]) + Square(u0);
      super[k] = alphaR;

      if (unorm2 > T(0))
      {
        T invNorm = T(1) / Sqrt(unorm2);
        A[k][k+1] = u0 * invNorm;
        for (I32 j = k+2; j < N; j++)
          A[k][j] *= invNorm;

        for (I32 i = k+1; i < M; i++)
        {
          T dot = T(0);
          for (I32 j = k+1; j < N; j++)
            dot += A[k][j] * A[i][j];
          dot *= T(2);
          for (I32 j = k+1; j < N; j++)
            A[i][j] -= dot * A[k][j];
        }
      }
      else
      {
        A[k][k+1] = T(0);
      }
    }
    else if (k + 1 < N)
    {
      super[k] = A[k][k+1];
    }
  }

  // Build the bidiagonal B (N x N) and SVD it via JacobiSVDTransposed.
  // JacobiSVDTransposed expects B^T as input and on output rows of Bt are columns of U_b.
  DynamicMatrix<T> Bt(N, N);
  Bt.Zeros();
  for (I32 k = 0; k < N; k++)
  {
    Bt[k][k] = diag[k];           // (B^T)[k][k] = B[k][k]
    if (k + 1 < N)
      Bt[k+1][k] = super[k];      // (B^T)[k+1][k] = B[k][k+1]
  }

  DynamicMatrix<T> Vb;
  bool converged = JacobiSVDTransposed(Bt, W, Vb);  // descending sort by default

  // After the call: Bt[i] (length N) is the i-th column of U_b, and Vb is V_b^T.
  // Form U = U_h * U_b, where U_b extends to M x N by zero-padding rows N..M-1.
  // We store U back into A. Use a temporary U buffer (M x N).
  DynamicMatrix<T> U(M, N);
  U.Zeros();
  for (I32 i = 0; i < N; i++)
    for (I32 j = 0; j < N; j++)
      U[i][j] = Bt[j][i];   // (U_b)[i][j] = (Bt[j])[i]

  for (I32 k = N-1; k >= 0; k--)
  {
    bool hasReflector = (A[k][k] != T(0));
    if (!hasReflector)
    {
      for (I32 i = k+1; i < M && !hasReflector; i++)
        if (A[i][k] != T(0))
          hasReflector = true;
      if (!hasReflector)
        continue;
    }

    for (I32 j = 0; j < N; j++)
    {
      T dot = T(0);
      for (I32 i = k; i < M; i++)
        dot += A[i][k] * U[i][j];
      dot *= T(2);
      for (I32 i = k; i < M; i++)
        U[i][j] -= dot * A[i][k];
    }
  }

  // Form V = V_h * V_b. Vb on output is V_b^T, so store V^T in V (so V holds V_b^T initially)
  // then apply right reflectors from the left to V (since (V_h * V_b)^T = V_b^T * V_h^T and
  // V_h^T = V_h for Householder reflections).
  // Easiest: convert Vb (which is V_b^T, N x N) into Vfull = V_b (N x N), then apply
  // V = V_h * V_b column-by-column (or row-by-row), and finally store V.
  for (I32 i = 0; i < N; i++)
    for (I32 j = 0; j < N; j++)
      V[i][j] = Vb[j][i];   // V starts as V_b

  for (I32 k = N-3; k >= 0; k--)
  {
    bool hasReflector = (A[k][k+1] != T(0));
    if (!hasReflector)
    {
      for (I32 j = k+2; j < N && !hasReflector; j++)
        if (A[k][j] != T(0))
          hasReflector = true;
      if (!hasReflector)
        continue;
    }

    for (I32 j = 0; j < N; j++)
    {
      T dot = T(0);
      for (I32 i = k+1; i < N; i++)
        dot += A[k][i] * V[i][j];
      dot *= T(2);
      for (I32 i = k+1; i < N; i++)
        V[i][j] -= dot * A[k][i];
    }
  }

  // Copy U back into A.
  for (I32 i = 0; i < M; i++)
    for (I32 j = 0; j < N; j++)
      A[i][j] = U[i][j];

  return converged;
}

//
// Transposed-form HouseholderJacobiSVD that operates in place on At (N x M, M >= N).
// All heavy inner loops run over rows of At (length M) which are contiguous in memory, so this
// is much faster than the column-oriented HouseholderJacobiSVD on tall matrices.
//
// On output:
//   At[i] (length M) is the i-th column of U (rows of At = columns of U), filling only the first
//        (effective) M positions; rows beyond N are not part of U_b but get zeroed by the
//        algorithm where needed.
//   W[i] = i-th singular value (descending).
//   Vt is V^T (N x N).
//
template <class T>
static bool HouseholderJacobiSVDTransposed(DynamicMatrix<T>& At,
                                           DynamicVector<T>& W,
                                           DynamicMatrix<T>& Vt)
{
  const I32 N = At.Rows();
  const I32 M = At.Cols();
  assert(M >= N);

  W.Resize(N);
  Vt.Resize(N, N);

  DynamicVector<T> diag(N);
  DynamicVector<T> super(N);
  DynamicVector<T> dotPerCol(M);   // workspace for right reflector update

  const I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

  // Bidiagonalize. For each k:
  //   Left reflector  L_k zeros A[*][k] below row k. Vector v_k stored contiguously in At[k][k..M-1].
  //   Right reflector R_k zeros A[k][*] right of column k+1. Vector u_k stored strided in
  //                  At[k+1..N-1][k].
  for (I32 k = 0; k < N; k++)
  {
    // ---- Left reflector on column k of A (= row k of At), entries k..M-1 ----
    const I32 lenL = M - k;
    T sigma;
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
    sigma = SumOfSquares_StreamUnaligned_Sequential(At[k] + k, lenL);
#else
    sigma = SumOfSquares_Sequential(At[k] + k, lenL);
#endif

    T alpha = Sqrt(sigma);
    if (At[k][k] > T(0))
      alpha = -alpha;

    T v0 = At[k][k] - alpha;
    T vnorm2 = sigma - Square(At[k][k]) + Square(v0);
    diag[k] = alpha;

    if (vnorm2 > T(0))
    {
      T invNorm = T(1) / Sqrt(vnorm2);
      At[k][k] = v0 * invNorm;
      if (lenL > 1)
      {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
        ScalarMultiplication_StreamUnaligned_Sequential(At[k] + k + 1, invNorm, lenL - 1);
#else
        ScalarMultiplication_Sequential(At[k] + k + 1, invNorm, lenL - 1);
#endif
      }

      // Update At[i][j] for i > k, j >= k. Each row i is independent — parallelize.
      const I32 rowsToUpdate = N - 1 - k;
      if (rowsToUpdate > 0)
      {
        const T* vk = At[k] + k;
        #pragma omp parallel for num_threads(numberOfThreads) schedule(static) if(rowsToUpdate > 1)
        for (I32 i = k + 1; i < N; i++)
        {
          T* row = At[i] + k;
          T dot;
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
          dot = DotProduct_StreamUnaligned_Sequential(vk, row, lenL);
          AdditionScaled_StreamUnaligned_Sequential(row, vk, T(-2) * dot, lenL);
#else
          dot = DotProduct_Sequential(vk, row, lenL);
          for (I32 j = 0; j < lenL; j++)
            row[j] -= T(2) * dot * vk[j];
#endif
        }
      }
    }
    else
    {
      At[k][k] = T(0);
    }

    // ---- Right reflector on row k of A (= column k of At), entries rows k+1..N-1 ----
    if (k + 2 < N)
    {
      T sigmaR = T(0);
      for (I32 i = k+1; i < N; i++)
        sigmaR += Square(At[i][k]);

      T alphaR = Sqrt(sigmaR);
      if (At[k+1][k] > T(0))
        alphaR = -alphaR;

      T u0 = At[k+1][k] - alphaR;
      T unorm2 = sigmaR - Square(At[k+1][k]) + Square(u0);
      super[k] = alphaR;

      if (unorm2 > T(0))
      {
        T invNorm = T(1) / Sqrt(unorm2);
        At[k+1][k] = u0 * invNorm;
        for (I32 i = k+2; i < N; i++)
          At[i][k] *= invNorm;

        // Update At[i][j] for i = k+1..N-1, j = k+1..M-1.
        // Two-pass: dotPerCol[j] = sum_{i} u[i] * At[i][j], then At[i][j] -= 2 u[i] dotPerCol[j].
        const I32 lenR = M - k - 1;
        memset(dotPerCol.Begin() + k + 1, 0, lenR * sizeof(T));
        for (I32 i = k+1; i < N; i++)
        {
          T u = At[i][k];
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
          AdditionScaled_StreamUnaligned_Sequential(dotPerCol.Begin() + k + 1, At[i] + k + 1,
                                                    u, lenR);
#else
          for (I32 j = k+1; j < M; j++)
            dotPerCol[j] += u * At[i][j];
#endif
        }
        for (I32 i = k+1; i < N; i++)
        {
          T u2 = T(-2) * At[i][k];
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
          AdditionScaled_StreamUnaligned_Sequential(At[i] + k + 1, dotPerCol.Begin() + k + 1,
                                                    u2, lenR);
#else
          for (I32 j = k+1; j < M; j++)
            At[i][j] += u2 * dotPerCol[j];
#endif
        }
      }
      else
      {
        At[k+1][k] = T(0);
      }
    }
    else if (k + 1 < N)
    {
      super[k] = At[k+1][k];
    }
  }

  // SVD the bidiagonal block via JacobiSVDTransposed on Bt (N x N).
  DynamicMatrix<T> Bt(N, N);
  Bt.Zeros();
  for (I32 k = 0; k < N; k++)
  {
    Bt[k][k] = diag[k];
    if (k + 1 < N)
      Bt[k+1][k] = super[k];
  }

  DynamicMatrix<T> Vb;
  bool converged = JacobiSVDTransposed(Bt, W, Vb);   // descending sort

  // After call: Bt[i] (N-vector) is the i-th column of U_b. Vb is V_b^T.

  // Snapshot reflectors before we overwrite At.
  DynamicMatrix<T> L(N, M);   // L[k] holds left reflector v_k (entries k..M-1 are nonzero)
  for (I32 k = 0; k < N; k++)
    memcpy(L[k], At[k], M * sizeof(T));

  DynamicMatrix<T> R(N, N);   // R[k] holds right reflector u_k (entries k+1..N-1)
  for (I32 k = 0; k < N-2; k++)
  {
    for (I32 i = 0; i < N; i++)
      R[k][i] = T(0);
    for (I32 i = k+1; i < N; i++)
      R[k][i] = At[i][k];
  }

  // Build initial Ut from U_b extended with zeros: Ut[i][j] = (U_b)[j][i] = Bt[i][j] for j<N,
  // and Ut[i][j] = 0 for j>=N.
  for (I32 i = 0; i < N; i++)
  {
    memcpy(At[i], Bt[i], N * sizeof(T));
    if (M > N)
      memset(At[i] + N, 0, (M - N) * sizeof(T));
  }

  // Apply left reflectors back: Ut <- Ut * L_k for k = N-1 down to 0.
  // For each row i of Ut: dot = <v_k, Ut[i]>; Ut[i] -= 2*dot*v_k (over indices k..M-1).
  for (I32 k = N-1; k >= 0; k--)
  {
    const I32 lenL = M - k;
    const T* vk = L[k] + k;

    bool any = false;
    for (I32 j = 0; j < lenL && !any; j++)
      if (vk[j] != T(0)) any = true;
    if (!any) continue;

    #pragma omp parallel for num_threads(numberOfThreads) schedule(static) if(N > 1)
    for (I32 i = 0; i < N; i++)
    {
      T* row = At[i] + k;
      T dot;
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
      dot = DotProduct_StreamUnaligned_Sequential(vk, row, lenL);
      AdditionScaled_StreamUnaligned_Sequential(row, vk, T(-2) * dot, lenL);
#else
      dot = DotProduct_Sequential(vk, row, lenL);
      for (I32 j = 0; j < lenL; j++)
        row[j] -= T(2) * dot * vk[j];
#endif
    }
  }

  // Build Vt = V^T = V_b^T * V_h. Start from Vb (= V_b^T) then apply each R_k from the right.
  for (I32 i = 0; i < N; i++)
    memcpy(Vt[i], Vb[i], N * sizeof(T));

  for (I32 k = N-3; k >= 0; k--)
  {
    const T* u = R[k];
    bool any = false;
    for (I32 j = k+1; j < N && !any; j++)
      if (u[j] != T(0)) any = true;
    if (!any) continue;

    for (I32 i = 0; i < N; i++)
    {
      T* row = Vt[i];
      T dot = T(0);
      for (I32 j = k+1; j < N; j++)
        dot += u[j] * row[j];
      dot *= T(2);
      for (I32 j = k+1; j < N; j++)
        row[j] -= dot * u[j];
    }
  }

  return converged;
}

//
// Solves A * x = b using HouseholderJacobiSVD. At is A transposed (N x M, with M >= N).
//
template <class T>
MTL_INLINE static bool SolveHouseholderJacobiSVDTransposed(DynamicVector<T>& x, I32& rank,
                                                           T& conditionNumber,
                                                           const DynamicMatrix<T>& At,
                                                           const DynamicVector<T>& b,
                                                           const T& tolerance = T(-1.0))
{
  assert(At.Cols() == (I32)b.Size());

  const I32 M = At.Cols();
  const I32 N = At.Rows();
  assert(M >= N);

  DynamicMatrix<T> Ut = At;
  DynamicMatrix<T> Vt;
  DynamicVector<T> D;

  bool fullyConverged = HouseholderJacobiSVDTransposed(Ut, D, Vt);

  T tol = tolerance;
  if (tol < 0)
    tol = M * Epsilon<T>() * D[0];

  rank = (I32)D.Size();
  for (; rank > 0 && D[rank-1] < tol; rank--);

  assert(rank > 0);

  // x = V * diag(1/D) * U^T * b. Ut[i] (length M) = i-th column of U => (U^T b)[i] = <Ut[i], b>.
  DynamicVector<T> Ub(rank);
  for (I32 i = 0; i < rank; i++)
  {
    const T* row = Ut[i];
    T sum = T(0);
    for (I32 j = 0; j < M; j++)
      sum += row[j] * b[j];
    Ub[i] = sum / D[i];
  }

  // V = Vt^T, so x[i] = sum_k Vt[k][i] * Ub[k].
  x.Resize(N);
  for (I32 i = 0; i < N; i++)
    x[i] = T(0);
  for (I32 k = 0; k < rank; k++)
  {
    const T* vrow = Vt[k];
    T uk = Ub[k];
    for (I32 i = 0; i < N; i++)
      x[i] += vrow[i] * uk;
  }

  conditionNumber = D[0] / D[rank-1];
  return fullyConverged;
}

template <class T>
MTL_INLINE static bool SolveJacobiSVDTransposedHomogeneous
(DynamicVector<T>& x, DynamicMatrix<T>& At, I32& rank, T& conditionNumber, const T& tol = T(-1.0))
{
  DynamicMatrix<T> Vt(At.Rows(), At.Rows());
  DynamicVector<T> D(At.Rows());

  bool fullyConverged = JacobiSVDTransposed(At[0], D.Begin(), Vt[0], At.Cols(), At.Rows(),
                                            At.RowSize(), Vt.RowSize());

  T tolerance = tol;
  if (tolerance < 0)
    tolerance = At.Cols() * Epsilon<T>() * D[0];

  rank = (I32)D.Size();
  for(; rank > 0 && D[rank-1] < tolerance; rank--);

  x.Assign(Vt[Vt.Rows()-1], Vt.Cols());

  conditionNumber = D[0] / D[rank-1];

  return fullyConverged;
}
template <int N, class T>
MTL_INLINE static bool SolveJacobiSVDTransposedHomogeneous
(ColumnVector<N,T>& x, DynamicMatrix<T>& At, I32& rank, T& conditionNumber, const T& tol = T(-1.0))
{
  assert(N == At.Rows());

  SquareMatrix<N,T> Vt;
  T D[N];

  bool fullyConverged = JacobiSVDTransposed(At[0], D, Vt[0], At.Cols(), N,
                                            At.RowSize(), Vt.RowSize());

  T tolerance = tol;
  if (tolerance < 0)
    tolerance = At.Cols() * Epsilon<T>() * D[0];

  rank = ComputeRankFromSingularValues<N,T>(D, tolerance);

  memcpy(&x[0], Vt[N-1], sizeof(x));

  conditionNumber = D[0] / D[rank-1];

  return fullyConverged;
}

}  // namespace MTL

#endif // MTL_SVD_H
