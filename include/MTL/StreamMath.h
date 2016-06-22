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

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX

#ifndef MTL_STREAM_MATH_H
#define MTL_STREAM_MATH_H

#include <MTL/SSE.h>

namespace MTL
{

template <class T> class X128;

MTL_INLINE static F32 Abs(const F32& a)
{
  X128<F32> X;
  X.LoadSingle(&a);
  return Abs(X)[0];
}
MTL_INLINE static F64 Abs(const F64& a)
{
  X128<F64> X;
  X.LoadSingle(&a);
  return Abs(X)[0];
}

template <class T>
MTL_INLINE static T Sqrt(const T& a)
{
  X128<T> X;
  X.LoadSingle(&a);
  return X.SqrtSingle()[0];
}

}  // namespace MTL

#endif  // MTL_STREAM_MATH_H

#endif  // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX
