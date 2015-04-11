//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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

#include <iostream>
#include <string>

namespace MTL
{

typedef signed   char       I8;
typedef          short     I16;
typedef          long      I32;
typedef          long long I64;
typedef unsigned char       U8;
typedef unsigned short     U16;
typedef unsigned long      U32;
typedef unsigned long long U64;

typedef float  F32;
typedef double F64;

typedef size_t SizeType;

typedef std::wstring  String;
typedef std::wostream OutputStream;

#ifndef MTL_INLINE
  #define MTL_INLINE inline
#endif

#define TOWCHAR(x)   L ## x
#define MTL_FILE(x)  TOWCHAR(x)
#define MTL__FILE__  MTL_FILE(__FILE__)

}  // namespace MTL


#endif  // MTL_DEFINITIONS_H
