//
// Math Template Library
//
// Copyright (c) 2022: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_BYTE_SWAP_H
#define MTL_BYTE_SWAP_H

#include<MTL/Exception.h>

namespace MTL
{
#ifdef __GNUC__
template <typename T> MTL_INLINE T _bswap_x86(T bytes)
{
  asm("bswap %0" : "=r" (bytes));
  return bytes;
}
template <> MTL_INLINE uint16_t _bswap_x86(uint16_t bytes)
{
  asm("xchg %%al, %%ah" : "=a" (bytes));
  return bytes;
}
#elif defined(WIN32)
MTL_INLINE uint16_t _vs_byte_swap(uint16_t bytes) { return _byteswap_ushort(bytes); }
MTL_INLINE uint32_t _vs_byte_swap(uint32_t bytes) { return _byteswap_ulong (bytes); }
MTL_INLINE uint64_t _vs_byte_swap(uint64_t bytes) { return _byteswap_uint64(bytes); }
#endif

template <typename T> MTL_INLINE T ByteSwapOptimized(T bytes)
{
#ifdef __GNUC__
  return _bswap_x86(bytes);
#elif defined(WIN32)
  return _vs_byte_swap(bytes);
#endif
}

template <typename T> MTL_INLINE static T ByteSwap(T bytes)
{
  if (sizeof(T) == 2)
    return (T)ByteSwapOptimized((uint16_t&)bytes);
  else if (sizeof(T) == 4)
    return (T)ByteSwapOptimized((uint32_t&)bytes);
  else if (sizeof(T) == 8)
    return (T)ByteSwapOptimized((uint64_t&)bytes);
  else
    MTL_THROW("Unsupported number of bytes for input! " + std::to_string(sizeof(T)));
}

}
#endif  // MTL_BYTE_SWAP_H
