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

#include "DynamicMatrix.h"
#include "Givens.h"
#include "Point2D.h"

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

template<class T>
static bool JacobiRotationsTransposed(T& c, T& s, I32 i, I32 j, I32 iteration, T* At, T* W,
                                      I32 M, I32 N, I32 rowSizeA)
{
  const T epsilon = Epsilon<T>();
  const T epsilon2 = Square(epsilon);

  T* Ai = At + i*rowSizeA;
  T* Aj = At + j*rowSizeA;
  T a = W[i];
  T b = W[j];
        
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
  T p = DotProduct_StreamAligned_Sequential(Ai, Aj, M);
#else
  T p = DotProduct_Sequential(Ai, Aj, Ai + M);
#endif
  T pp = p*p;

  if (pp <= epsilon2*a*b)
    return false;

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
    return false;

  if (iteration < 10)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    GivensRotation_StreamAligned_Sequential(Ai, Aj, c, s, M);
#else
    GivensRotation_Sequential(Ai, Aj, c, s, Ai + M);
#endif
    W[i] += delta;
    W[j] -= delta;
  }
  else
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    GivensRotation_StreamAligned_Sequential(Ai, Aj, c, s, M, W[i], W[j]);
#else
    GivensRotation_Sequential(Ai, Aj, c, s, Ai + M, W[i], W[j]);
#endif
    W[i] += delta * T(0.5);
    W[j] -= delta * T(0.5);
  }

  return true;
}

static void ComputeJacobiParallelPairs(DynamicVector<DynamicVector<Point2D<I32>>>& pairs, I32 N)
{
  DynamicVector<Point2D<I32>> allPairs;
  allPairs.Reserve((N-1)*N/2);
  I32 index = 0;
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
        //break;
      }
    }
  }
}

// Dynamic matrix version of SVD. Note that it takes A transposed. It expects At and Vt to be
// memory aligned.
template<class T>
static bool JacobiSVDTransposed(T* At, T* W, T* Vt, I32 M, I32 N, I32 rowSizeA, I32 rowSizeV)
{
  I32 maxIterations = Max(M, 30);
  I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

  OptimizedZeros(Vt, N*rowSizeV);

  for (I32 i = 0; i < N; i++)
  {
    Vt[i*rowSizeV + i] = T(1);
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    W[i] = SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M);
#else
    W[i] = SumOfSquares_Sequential(At + i*rowSizeA, At + i*rowSizeA + M);
#endif
  }

  DynamicVector<DynamicVector<Point2D<I32>>> parallelPairs;
  ComputeJacobiParallelPairs(parallelPairs, N);
  DynamicVector<I32> setChanged;

  I32 iteration;
  for (iteration = 0; iteration < maxIterations; iteration++)
  {
    bool changed = false;
    
    for (const auto& set : parallelPairs)
    {
      if (set.Size() > 1)
      {
        setChanged.Resize(set.Size());

        I32 blockSize = (I32)MTL::ComputeParallelSubSizesBlockSize(set.Size(), numberOfThreads);

        #pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic, blockSize)
        for (I32 k = 0; k < (I32)set.Size(); k++)
        {
          int i = set[k].x();
          int j = set[k].y();
          T c, s;

          setChanged[k] = JacobiRotationsTransposed(c, s, i, j, iteration, At, W, M, N, rowSizeA);
          if (setChanged[k])
          {

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
            GivensRotation_StreamAligned_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#else
            GivensRotation_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, Vt + i*rowSizeV + N);
#endif
            changed = true;
          }
        }

        changed |= Sum_StreamAligned_Sequential(setChanged.Begin(), setChanged.Size()) > 0;
      }
      else
      {
        int i = set[0].x();
        int j = set[0].y();
        T c, s;

        if (JacobiRotationsTransposed(c, s, i, j, iteration, At, W, M, N, rowSizeA))
        {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
          GivensRotation_StreamAligned_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, N);
#else
          GivensRotation_Sequential(Vt + i*rowSizeV, Vt + j*rowSizeV, c, s, Vt + i*rowSizeV + N);
#endif
          changed = true;
        }
      }
    }

    if (!changed)
      break;
  }

  for (I32 i = 0; i < N; i++)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    W[i] = Sqrt(SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M));
#else
    W[i] = Sqrt(SumOfSquares_Sequential(At + i*rowSizeA, At + i*rowSizeA + M));
#endif

  for (I32 i = 0; i < N; i++)
  {
    if (W[i] > 0)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarMultiplication_StreamAligned_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#else
      ScalarMultiplication_Sequential(At + i*rowSizeA, T(1)/W[i], At + i*rowSizeA + M);
#endif
  }

  return iteration < maxIterations;
}
template<class T>
static bool JacobiSVDTransposed(T* At, T* W, I32 M, I32 N, I32 rowSizeA)
{
  I32 maxIterations = Max(M, 30);
  I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();
  
  for (I32 i = 0; i < N; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    W[i] = SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M);
#else
    W[i] = SumOfSquares_Sequential(At + i*rowSizeA, At + i*rowSizeA + M);
#endif
  }

  DynamicVector<DynamicVector<Point2D<I32>>> parallelPairs;
  ComputeJacobiParallelPairs(parallelPairs, N);
  DynamicVector<I32> setChanged;

  I32 iteration;
  for (iteration = 0; iteration < maxIterations; iteration++)
  {
    bool changed = false;
    
    for (const auto& set : parallelPairs)
    {
      if (set.Size() > 1)
      {
        setChanged.Resize(set.Size());

        I32 blockSize = (I32)MTL::ComputeParallelSubSizesBlockSize(set.Size(), numberOfThreads);

        #pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic, blockSize)
        for (I32 i = 0; i < (I32)set.Size(); i++)
        {
          T c, s;
          setChanged[i] = JacobiRotationsTransposed(c, s, set[i].x(), set[i].y(), iteration,
                                                    At, W, M, N, rowSizeA);
        }

        changed |= Sum_StreamAligned_Sequential(setChanged.Begin(), setChanged.Size()) > 0;
      }
      else
      {
        T c, s;
        if (JacobiRotationsTransposed(c, s, set[0].x(), set[0].y(), iteration,
                                      At, W, M, N, rowSizeA))
          changed = true;
      }
    }

    if (!changed)
      break;
  }

  //#pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic, blockSize)
  for (I32 i = 0; i < N; i++)
  {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    W[i] = Sqrt(SumOfSquares_StreamAligned_Sequential(At + i*rowSizeA, M));
#else
    W[i] = Sqrt(SumOfSquares_Sequential(At + i*rowSizeA, At + i*rowSizeA + M));
#endif

    if (W[i] > 0)
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      ScalarMultiplication_StreamAligned_Sequential(At + i*rowSizeA, T(1)/W[i], M);
#else
      ScalarMultiplication_Sequential(At + i*rowSizeA, T(1)/W[i], At + i*rowSizeA + M);
#endif
  }

  return iteration < maxIterations;
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
MTL_INLINE static bool JacobiSVDTransposed(DynamicMatrix<T>& Ut,
                                           DynamicVector<T>& D,
                                           DynamicMatrix<T>& Vt,
                                           bool sortAscending = false)
{
  Vt.Resize(Ut.Rows(), Ut.Rows());
  D.Resize(Ut.Rows());
  bool converged = JacobiSVDTransposed(Ut[0], D.Begin(), Vt[0], Ut.Cols(), Ut.Rows(),
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
MTL_INLINE static bool JacobiSVDTransposed(DynamicMatrix<T>& Ut,
                                           DynamicVector<T>& D,
                                           bool sortAscending = false)
{
  D.Resize(Ut.Rows());
  bool converged = JacobiSVDTransposed(Ut[0], D.Begin(), Ut.Cols(), Ut.Rows(), Ut.RowSize());

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
MTL_INLINE static bool JacobiEigenTransposed(DynamicMatrix<T>& Ut,
                                             DynamicVector<T>& D,
                                             bool sortAscending = true)
{
  return JacobiSVDTransposed(Ut, D, sortAscending);
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
MTL_INLINE static T SolveEigenTransposed(DynamicVector<T>& x,
                                         DynamicMatrix<T>& Ut,
                                         DynamicVector<T>& D,
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
  D.Resize(rank);

  DynamicVector<T> X = Ut * b;
  X /= D;
  x = Ut.ComputeTranspose() * X;

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
  assert(At.Cols() == b.Size());

  DynamicMatrix<T> Ut = At;
  DynamicMatrix<T> Vt;
  DynamicVector<T> D;

  bool fullyConverged = JacobiSVDTransposed(Ut, D, Vt);
  conditionNumber = SolveSVDTransposed(x, Ut, D, Vt, rank, b, tolerance);

  return fullyConverged;
}
template <class T>
MTL_INLINE static bool SolveJacobiEigenTransposed(DynamicVector<T>& x, I32& rank,
                                                  T& conditionNumber,
                                                  const DynamicMatrix<T>& At,
                                                  const DynamicVector<T>& b,
                                                  const T& tolerance = T(-1.0))
{
  assert(At.Cols() == b.Size());

  DynamicMatrix<T> Ut = At;
  DynamicVector<T> D;

  bool fullyConverged = JacobiSVDTransposed(Ut, D);
  conditionNumber = SolveEigenTransposed(x, Ut, D, rank, b, tolerance);

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
(ColumnVector<N,T>& x, DynamicMatrix<T>& At,
 I32& rank, T& conditionNumber, const T& tol = T(-1.0))
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
