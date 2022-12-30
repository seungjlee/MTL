//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_DEFINITIONS_H
#define MTL_DEFINITIONS_H

#include <cstdint>
#include <iostream>
#include <string>
#include <MTL/Macros.h>

namespace MTL
{

typedef int8_t     I8;
typedef int16_t   I16;
typedef int32_t   I32;
typedef int64_t   I64;
typedef uint8_t    U8;
typedef uint16_t  U16;
typedef uint32_t  U32;
typedef uint64_t  U64;

typedef float  F32;
typedef double F64;

typedef size_t SizeType;

#define ConsoleOut std::wcout
using String = std::wstring;
using OutputStream = std::wostream;
}  // namespace MTL

#ifndef MTL_INLINE
  #define MTL_INLINE inline
#endif

#define MTL__FILE__   TO_WCHAR(__FILE__)

// Default definitions for stream instructions.
#if !defined(MTL_ENABLE_SSE) && !defined(MTL_ENABLE_AVX) && !defined(MTL_ENABLE_AVX512)
#define MTL_ENABLE_SSE    1
#define MTL_ENABLE_AVX    0
#define MTL_ENABLE_AVX512 0
#endif

#if defined(MTL_ENABLE_AVX) && MTL_ENABLE_AVX
  #ifndef MTL_ENABLE_SSE
    #define MTL_ENABLE_SSE 1
  #endif
#endif

#if !defined(MTL_ENABLE_AVX)
#define MTL_ENABLE_AVX 0
#endif

#if !defined(MTL_ENABLE_AVX512)
#define MTL_ENABLE_AVX512 0
#endif

// Note that I have not really tested builds with OpenMP disabled.
#ifndef MTL_ENABLE_OPENMP
  #define MTL_ENABLE_OPENMP 1
#endif

#define MTL_THROW(MSG)                                                                        \
{                                                                                             \
  char _msg_[1024];                                                                           \
  std::snprintf(_msg_, sizeof(_msg_), "%s, %s, line %d: ", __FUNCTION__, __FILE__, __LINE__); \
  throw MTL::Exception(std::string(_msg_) + MSG);                                             \
}

#ifdef WIN32

#include <corecrt.h>
#pragma warning(push)
#pragma warning(disable: _UCRT_DISABLED_WARNINGS)
_UCRT_DISABLE_CLANG_WARNINGS
_CRT_BEGIN_C_HEADER
_ACRTIMP void __cdecl _wassert(
  _In_z_ wchar_t const* _Message,
  _In_z_ wchar_t const* _File,
  _In_   unsigned       _Line
);
_CRT_END_C_HEADER

#define MTL_ALWAYS_ASSERT(expression) (void)(                                                 \
            (!!(expression)) ||                                                               \
            (_wassert(_CRT_WIDE(#expression), _CRT_WIDE(__FILE__), (unsigned)(__LINE__)), 0))
_UCRT_RESTORE_CLANG_WARNINGS
#pragma warning(pop) // _UCRT_DISABLED_WARNINGS

#else

extern "C" void __assert_fail(const char *__assertion, const char *__file,
  unsigned int __line, const char *__function)
  __THROW __attribute__((__noreturn__));
#define MTL_ALWAYS_ASSERT(expr)							\
     (static_cast <bool> (expr)						  \
      ? void (0)							              \
      : __assert_fail (#expr, __FILE__, __LINE__, __extension__ __PRETTY_FUNCTION__))

#endif

#endif  // MTL_DEFINITIONS_H
