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
template <typename T> MTL_INLINE T _bwap_x86(T bytes)
{
  asm("bswap %0" : "=r" (bytes));
  return bytes;
}

template <typename T> MTL_INLINE static T ByteSwap(T bytes)
{
  if (sizeof(T) == 2)
    return _bwap_x86((uint16_t&)bytes);
  else if (sizeof(T) == 4)
    return _bwap_x86((uint32_t&)bytes);
  else if (sizeof(T) == 8)
    return _bwap_x86((uint64_t&)bytes);
  else
    MTL_THROW("Unsupported number of bytes for input! " + std::to_string(sizeof(T)));
}

}
#endif  // MTL_BYTE_SWAP_H
