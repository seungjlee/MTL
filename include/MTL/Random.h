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


#ifndef MTL_RANDOM_H
#define MTL_RANDOM_H

#include "DynamicMatrix.h"

namespace MTL
{

template<U64 Multiplier, U64 Increment>
class RandomLinearCongruential
{
public:
  MTL_INLINE RandomLinearCongruential(U64 seed = 0)
    : x_(seed) {}

  MTL_INLINE U32 GetNext()
  {
    x_ = x_ * Multiplier + Increment;

    return U32(x_);
  }

  template<class T>
  MTL_INLINE T GetNext(T start, T end)
  {
    return GetNextInRange(start, end - start);
  }

  template<class T>
  MTL_INLINE T GetNextInRange(T start, T range)
  {
    return (T(GetNext()) / T(U32(-1))) * range + start;
  }

  template<class T>
  MTL_INLINE MTL::DynamicMatrix<T> DynamicMatrix(I32 M, I32 N, T start, T end)
  {
    MTL::DynamicMatrix<T> m(M, N);
    for (I32 i = 0; i < m.Rows(); i++)
	  Array(m[i], N, start, end);

    return m;
  }

  template<class T>
  MTL_INLINE MTL::DynamicVector<T> DynamicVector(SizeType size, T start, T end)
  {
    MTL::DynamicVector<T> v(size);
    Array(&v[0], size, start, end);

    return v;
  }

  template<I32 M, I32 N, class T>
  MTL_INLINE MTL::Matrix<M,N,T> Matrix(T start, T end)
  {
    MTL::Matrix<M,N,T> m;
    Array(m[0], M*N, start, end);

    return m;
  }

  template<class T>
  MTL_INLINE void Array(T* pDst, SizeType size, T start, T end)
  {
    Array(pDst, pDst + size, start, end);
  }

  template<class T>
  MTL_INLINE void Array(T* pDst, const T* pDstEnd, T start, T end)
  {
    T range = end - start;

    for (; pDst < pDstEnd; pDst++)
      *pDst = GetNextInRange(start, range);
  }

  MTL_INLINE MTL::DynamicMatrix<I32> DynamicMatrixI32(I32 M, I32 N)
  {
    MTL::DynamicMatrix<I32> m(M, N);
    for (I32 i = 0; i < m.Rows(); i++)
      ArrayI32(m[i], N);

    return m;
  }
  MTL_INLINE MTL::DynamicVector<I32> DynamicVectorI32(SizeType size)
  {
    MTL::DynamicVector<T> v(size);
    ArrayI32(&v[0], size);

    return v;
  }
  MTL_INLINE void ArrayI32(I32* pDst, SizeType size)
  {
    Array(pDst, pDst + size, start, end);
  }
  MTL_INLINE void ArrayI32(I32* pDst, const I32* pDstEnd)
  {
    for (; pDst < pDstEnd; pDst++)
      *pDst = GetNext();
  }

private:
  U64 x_;
};

//
//  The multiplier and increment are from MMIX by Donald Knuth.
//  I wrote the Python code below to verify the randomness using these 2 numbers.
//
//    from time import clock
//    from numpy import arange
//    from numpy import uint32, uint64
//    from matplotlib.pyplot import show, hist
//
//    N = 100000000
//    X = uint32(arange(N))
//    X.ravel()
//
//    a = uint64(6364136223846793005)
//    c = uint64(1442695040888963407)
//
//    Seed = uint64(0)
//    x = Seed
//    start = clock()
//    for i in range(0, N-1):
//      x = x * a + c
//      X[i] = uint32(x)
//    elapsed = (clock() - start)
//    print('Random generation time: {0:.3f} secs'.format(elapsed))
//
//    start = clock()
//    hist(X,1000)
//    elapsed = (clock() - start)
//    print('Histogram time: {0:.3f} secs'.format(elapsed))
//    show()
//
typedef RandomLinearCongruential<6364136223846793005, 1442695040888963407> Random;

}  // namespace MTL

#endif // MTL_RANDOM_H
