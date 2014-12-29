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


#ifndef MTL_STREAM_H
#define MTL_STREAM_H

#if !defined(MTL_ENABLE_SSE) && !defined(MTL_ENABLE_AVX)
#define MTL_ENABLE_SSE 1
#define MTL_ENABLE_AVX 0
#endif

#if defined(MTL_ENABLE_AVX) && MTL_ENABLE_AVX
  #ifndef MTL_ENABLE_SSE
    #define MTL_ENABLE_SSE 1
  #endif
#endif

#if MTL_ENABLE_AVX
  #define MTL_STREAM_BITS 256
#elif MTL_ENABLE_SSE
  #define MTL_STREAM_BITS 128
#else
  #define MTL_STREAM_BITS 8
#endif

#define MTL_STREAM_BYTES MTL_STREAM_BITS / 8

#define X__X(Bits)           X ## Bits
#define X__X_SetPacked(Bits) X ## Bits ## _SetPacked
#define X_X(Bits)            X__X(Bits)
#define X_X_SetPacked(Bits)  X__X_SetPacked(Bits)

#define XX           X_X(MTL_STREAM_BITS)
#define XX_SetPacked X_X_SetPacked(MTL_STREAM_BITS)


#define MTL_STREAM_EXTRA_INTEGER_OPERATORS(Bits)                                       \
MTL_INLINE X__X(Bits)& operator+=(const X__X(Bits)& y)  { return *this = *this + y; }  \
MTL_INLINE X__X(Bits)& operator-=(const X__X(Bits)& y)  { return *this = *this - y; }  \
MTL_INLINE X__X(Bits)& operator&=(const X__X(Bits)& y)  { return *this = *this & y; }  \
MTL_INLINE X__X(Bits)& operator|=(const X__X(Bits)& y)  { return *this = *this | y; }  \

#define MTL_STREAM_EXTRA_OPERATORS(Bits)                                               \
MTL_STREAM_EXTRA_INTEGER_OPERATORS(Bits)                                               \
MTL_INLINE X__X(Bits)& operator*=(const X__X(Bits)& y)  { return *this = *this * y; }  \
MTL_INLINE X__X(Bits)& operator/=(const X__X(Bits)& y)  { return *this = *this / y; }  \


#define FOR_STREAM(p, size)                                              \
  const T* pStreamEnd = p + MTL::XX<T>::StreamSize(size);                \
  for (; p < pStreamEnd; p += MTL::XX<T>::Increment)

#define FOR_STREAM_TYPE(p, size, Type)                                   \
  const Type* pStreamEnd = p + MTL::XX<Type>::StreamSize(size);          \
  for (; p < pStreamEnd; p += MTL::XX<Type>::Increment)

#define FOR_STREAM2(p1, p2, size)                                        \
  const T* pStreamEnd = p1 + MTL::XX<T>::StreamSize(size);               \
  for (; p1 < pStreamEnd; p1 += MTL::XX<T>::Increment,                   \
                          p2 += MTL::XX<T>::Increment)

#define FOR_STREAM2_TYPE(p1, p2, size, Type)                             \
  const Type* pStreamEnd = p1 + MTL::XX<Type>::StreamSize(size);         \
  for (; p1 < pStreamEnd; p1 += MTL::XX<Type>::Increment,                \
                          p2 += MTL::XX<Type>::Increment)

#define FOR_STREAM2_TYPE2(p1, p2, size, Type, Type2)                     \
  const Type* pStreamEnd = p1 + MTL::XX<Type2>::StreamSize(size);        \
  for (; p1 < pStreamEnd; p1 += MTL::XX<Type2>::Increment,               \
                          p2 += MTL::XX<Type2>::Increment)

#endif  // MTL_STREAM_H
