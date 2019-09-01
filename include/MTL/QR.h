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


#ifndef MTL_QR_H
#define MTL_QR_H

#include "DynamicMatrix.h"

namespace MTL
{

// Solves A * x = b. b is passed in as x. Note that x is returned in the top N of the vector.
// Returns the rank of matrix A.
template <I32 M, I32 N, class T>
static I32 SolveHouseholderQR(ColumnVector<M,T>& x, Matrix<M,N,T>& A,
                              const T& tolerance = M * Epsilon<T>())
{
  assert(M >= N);

  I32 P[N];
  for (I32 i = 0; i < N; i++)
    P[i] = i;

  I32 Q = Min(M-1, N);

  for (I32 i = 0; i < Q; i++)
  {
    // Column pivoting.
    I32 maxIndex = i;
    T max = Abs(A[i][P[maxIndex]]);
    for (I32 currentIndex = i+1; currentIndex < N; currentIndex++)
    {
      T current = Abs(A[i][P[currentIndex]]);
      maxIndex = current > max ? currentIndex : maxIndex;
      max = current > max ? current : max;
    }
    Swap(P[i], P[maxIndex]);

    I32 Mi = M - i;
    I32 Ni = N - i;

    T* pV = A[i] + P[i];
    T Aii = *pV;

    T sumOfSquares = Pow<2>(pV[N]);
    for (I32 k = 2; k < Mi; k++)
      sumOfSquares += Pow<2>(pV[k*N]);

    T normAii = Sqrt(Pow<2>(Aii) + sumOfSquares);
    *pV -= normAii;

    T norm = Sqrt(Pow<2>(*pV) + sumOfSquares);
    if (norm < tolerance)
    {
      *pV += normAii;  // Restore.
      continue;
    }

    T div = T(1)/norm;
    for (I32 k = 0; k < Mi; k++)
      pV[k*N] *= div;

    T dotAii = (Aii * (Aii - normAii) + sumOfSquares) * div;
    for (I32 j = i+1; j < N; j++)
    {
      T dotA = pV[0] * A[i][P[j]];
      for (I32 k = 1; k < Mi; k++)
        dotA += pV[k*N] * A[i+k][P[j]];

      for (I32 k = 0; k < Mi; k++)
        A[i+k][P[j]] += T(-2) * dotA * pV[k*N];
    }
    T dotB = pV[0] * x[i];
    for (I32 k = 1; k < Mi; k++)
      dotB += pV[k*N] * x[i+k];

    for (I32 k = 0; k < Mi; k++)
      x[i+k] += T(-2) * dotB * pV[k*N];

    *pV = Aii + *pV * T(-2) * dotAii;
  }

  I32 rank = 0;

  for (I32 j = N - 1; j >= 0; j--)
  {
    for (I32 k = j + 1; k < N; k++)
      x[j] -= A[j][P[k]] * x[k];

    T d = A[j][P[j]];

    if (Abs(d) < tolerance)
    {
      x[j] = 0;
    }
    else
    {
      x[j] /= d;
      rank++;
    }
  }

  T xx[N];
  memcpy(xx, &x[0], sizeof(xx));
  for (I32 i = 0; i < N; i++)
    x[P[i]] = xx[i];

  return rank;
}

// Solves A * x = b. b is passed in as x. Note that x is returned in the top N of the vector.
// Returns the rank of matrix A.
// Note that A is expected to be transposed as input for computational efficiency.
template <class T>
static I32 SolveHouseholderQRTransposed(T* x, T* At, I32 M, I32 N, I32 rowSize,
                                        const T& tolerance = Epsilon<T>())
{
  DynamicVector<I32> permutation(N);
  DynamicVector<T> permuteBuffer(N);
  return SolveHouseholderQRTransposed(x, At, M, N, rowSize, permutation.Begin(),
                                      permuteBuffer.Begin(), tolerance);
}
template <class T>
static I32 SolveHouseholderQRTransposed(T* x, T* At, I32 M, I32 N, I32 rowSize,
                                        I32* P, T* pQtb,
                                        const T& tolerance = Epsilon<T>())
{
  assert(M >= N);  // Not supporting or testing M < N cases.

  I32 numberOfThreads = (I32)MTL::CPU::Instance().NumberOfThreads();

  for (I32 i = 0; i < N; i++)
    P[i] = i;

  I32 Q = Min(M-1, N);

  for (I32 i = 0; i < Q; i++)
  {
    // Column pivoting.
    I32 maxIndex = i;
    T max = Abs(At[P[maxIndex]*rowSize + i]);
    for (I32 currentIndex = i+1; currentIndex < N; currentIndex++)
    {
      T current = Abs(At[P[currentIndex]*rowSize + i]);
      maxIndex = current > max ? currentIndex : maxIndex;
      max = current > max ? current : max;
    }
    Swap(P[i], P[maxIndex]);

    I32 Mi = M - i;

    T* pV = At + P[i]*rowSize + i;
    T Aii = *pV;
    
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    T sumOfSquares = SumOfSquares_StreamUnaligned_Parallel(pV + 1, Mi - 1, numberOfThreads);
#else
    T sumOfSquares = SumOfSquares_Sequential(pV + 1, pV + Mi);
#endif

    T normAii = Sqrt(Pow<2>(Aii) + sumOfSquares);
    *pV -= normAii;

    T norm = Sqrt(Pow<2>(*pV) + sumOfSquares);
    if (norm < tolerance)
    {
      *pV += normAii;  // Restore.
      continue;
    }

    T div = T(1)/norm;
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    ScalarMultiplication_StreamUnaligned_Parallel(pV, div, Mi, numberOfThreads);
#else
    ScalarMultiplication_Sequential(pV, div, pV + Mi);
#endif

    T dotAii = (Aii * (Aii - normAii) + sumOfSquares) * div;
    for (I32 j = i+1; j < N; j++)
    {
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
      T dotA = DotProduct_StreamUnaligned_Parallel(pV, At + P[j]*rowSize + i, Mi, numberOfThreads);
      AdditionScaled_StreamUnaligned_Parallel(At + P[j]*rowSize + i, pV, T(-2) * dotA, Mi,
                                              numberOfThreads);
#else
      T dotA = DotProduct_Sequential(pV, At + P[j]*rowSize + i, pV + Mi);
      AdditionScaled_Sequential(At + P[j]*rowSize + i, pV, T(-2) * dotA, At + P[j]*rowSize + i + Mi);
#endif
    }
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
    T dotB = DotProduct_StreamUnaligned_Parallel(pV, x + i, Mi, numberOfThreads);
    AdditionScaled_StreamUnaligned_Parallel(x + i, pV, T(-2) * dotB, Mi, numberOfThreads);
#else
    T dotB = DotProduct_Sequential(pV, x + i, pV + Mi);
    AdditionScaled_Sequential(x + i, pV, T(-2) * dotB, x + i + Mi);
#endif

    *pV = Aii + *pV * T(-2) * dotAii;
  }

  I32 rank = 0;

  for (I32 j = N - 1; j >= 0; j--)
  {
    for (I32 k = j + 1; k < N; k++)
      x[j] -= At[P[k] * rowSize + j] * x[k];

    T d = At[P[j] * rowSize + j];

    if (Abs(d) < tolerance)
    {
      x[j] = 0;
    }
    else
    {
      x[j] /= d;
      rank++;
    }
  }

  memcpy(pQtb, &x[0], N*sizeof(x[0]));
  for (I32 i = 0; i < N; i++)
    x[P[i]] = pQtb[i];

  return rank;
}
template <class T>
static I32 SolveHouseholderQRTransposed(DynamicVector<T>& x, DynamicMatrix<T>& At,
                                        DynamicVector<I32>& permutation,
                                        DynamicVector<T>& permuteBuffer,
                                        const T& tolerance = Epsilon<T>())
{
  I32 M = At.Cols();
  I32 N = At.Rows();
  assert(M == (I32)x.Size());

  permutation.Resize(N);
  permuteBuffer.Resize(N);

  I32 rank = SolveHouseholderQRTransposed(x.Begin(), At[0], M, N, At.RowSize(),
                                          permutation.Begin(), permuteBuffer.Begin(), tolerance);

  x.Resize(N);

  return rank;
}
template <class T>
static I32 SolveHouseholderQRTransposed(DynamicVector<T>& x,
                                        DynamicMatrix<T>& At,
                                        const T& tolerance = Epsilon<T>())
{
  I32 N = At.Rows();

  DynamicVector<I32> permutation(N);
  DynamicVector<T> permuteBuffer(N);

  return SolveHouseholderQRTransposed(x, At, permutation, permuteBuffer, tolerance);
}

}  // namespace MTL

#endif // MTL_QR_H
