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


#ifndef MTL_GIVENS_H
#define MTL_GIVENS_H

#include "DynamicVector.h"

namespace MTL
{

template<class T>
MTL_INLINE static void GivensRotation(T& a, T& b, const T& c, const T& s)
{
  T t = a;
  a = c*t + s*b;
  b = c*b - s*t;
}

template<I32 M, I32 N, class T>
MTL_INLINE static void GivensRotation(Matrix<M,N,T>& A, I32 i, I32 j, const T& c, const T& s)
{
  for (I32 k = 0; k < M; k++)
    GivensRotation(A[k][i], A[k][j], c, s);
}

template<I32 M, I32 N, class T>
MTL_INLINE static void GivensRotation(Matrix<M,N,T>& A, I32 i, I32 j, const T& c, const T& s,
                                      T& a, T& b)
{
  a = 0;
  b = 0;

  for (I32 k = 0; k < M; k++)
  {
    T tx = A[k][i];
    T ty = A[k][j];
    A[k][i] = c*tx + s*ty;
    A[k][j] = c*ty - s*tx;

    a += Square(tx);
    b += Square(ty);
  }
}

template<class T>
MTL_INLINE static void GivensRotation_Sequential(T* x, T* y, const T& c, const T& s, const T* xEnd)
{
  for (; x < xEnd; x++, y++)
    GivensRotation(*x, *y, c, s);
}
template<class T>
MTL_INLINE static void GivensRotation_Sequential(T* x, T* y, const T& c, const T& s, const T* xEnd,
                                                 T& a, T& b)
{
  for (; x < xEnd; x++, y++)
  {
    a += Square(*x);
    b += Square(*y);
    GivensRotation(*x, *y, c, s);
  }
}

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
template<class T>
MTL_INLINE static void GivensRotation_StreamAligned_Sequential(T* x, T* y, const T& c, const T& s,
                                                               SizeType N)
{
  const T* xEnd = x + N;

  if (XX<T>::DoSSE(N))
  {
    XX<T> xC(c);
    XX<T> xS(s);

    FOR_STREAM2(x, y, N)
    {
      XX<T> xX, xY;
      xX.LoadPackedAligned(x);
      xY.LoadPackedAligned(y);
      XX<T> xRX = xC*xX + xS*xY;
      XX<T> xRY = xC*xY - xS*xX;
      xRX.StorePackedAligned(x);
      xRY.StorePackedAligned(y);
    }
  }
  GivensRotation_Sequential(x, y, c, s, xEnd);
}
template<class T>
MTL_INLINE static void GivensRotation_StreamAligned_Sequential(T* x, T* y, const T& c, const T& s,
                                                               SizeType N, T& a, T& b)
{
  const T* xEnd = x + N;

  if (XX<T>::DoSSE(N))
  {
    XX<T> xA = XX<T>::Zeros();
    XX<T> xB = XX<T>::Zeros();

    XX<T> xC(c);
    XX<T> xS(s);

    FOR_STREAM2(x, y, N)
    {
      XX<T> xX, xY;
      xX.LoadPackedAligned(x);
      xY.LoadPackedAligned(y);
      XX<T> xRX = xC*xX + xS*xY;
      XX<T> xRY = xC*xY - xS*xX;
      xRX.StorePackedAligned(x);
      xRY.StorePackedAligned(y);

      xA += Square(xX);
      xB += Square(xY);
    }

    a = Sum< XX<T>::Increment >(xA.pData());
    b = Sum< XX<T>::Increment >(xB.pData());
  }

  T aa, bb;
  GivensRotation_Sequential(x, y, c, s, xEnd, aa, bb);

  a += aa;
  b += bb;
}

template<class T>
MTL_INLINE static void GivensRotation_StreamAligned_Parallel(T* x, T* y, const T& c, const T& s,
                                                             SizeType N)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(N, numberOfThreads))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, N, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      GivensRotation_StreamAligned_Sequential(x + offsets[i], y + offsets[i], c, s, subSizes[i]);
  }
  else
#endif
    GivensRotation_StreamAligned_Sequential(x, y, c, s, N);
}
template<class T>
MTL_INLINE static void GivensRotation_StreamAligned_Parallel(T* x, T* y, const T& c, const T& s,
                                                             SizeType N, T& a, T& b)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(N, numberOfThreads))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, N, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      GivensRotation_StreamAligned_Sequential(x + offsets[i], y + offsets[i], c, s,
                                              subSizes[i], a, b);
  }
  else
#endif
    GivensRotation_StreamAligned_Sequential(x, y, c, s, N, a, b);
}
#endif  // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX

}  // namespace MTL

#endif // MTL_GIVENS_H
