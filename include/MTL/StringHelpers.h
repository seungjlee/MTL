//
// Math Template Library
//
// Copyright (c) 2016-2019: Seung Jae Lee, https://github.com/seungjlee/MTL
//
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

#ifndef MTL_STRING_HELPERS_H
#define MTL_STRING_HELPERS_H

#include <MTL/Definitions.h>
#include <algorithm>
#include <wctype.h>
#include <codecvt>
#include <locale>

namespace MTL
{

// Call this to switch between std::cout (printf) and std::wcout (wprintf).
static MTL_INLINE void ResetOutputStream()
{
#ifdef WIN32
  FILE* f;
  errno_t n = freopen_s(&f, nullptr, "w", stdout);
#else
  FILE* f = freopen(nullptr, "w", stdout);
#endif
}

static MTL_INLINE void SetToLowerCase(std::wstring& str)
{
  std::transform(str.begin(), str.end(), str.begin(), towlower);
}
static MTL_INLINE void SetToLowerCase(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), tolower);
}
static MTL_INLINE void SetToUpperCase(std::wstring& str)
{
  std::transform(str.begin(), str.end(), str.begin(), towupper);
}
static MTL_INLINE void SetToUpperCase(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), toupper);
}

static MTL_INLINE std::wstring ToLowerCase(const std::wstring& input)
{
  std::wstring output = input;
  SetToLowerCase(output);
  return output;
}
static MTL_INLINE std::string ToLowerCase(const std::string& input)
{
  std::string output = input;
  SetToLowerCase(output);
  return output;
}
static MTL_INLINE std::wstring ToUpperCase(const std::wstring& input)
{
  std::wstring output = input;
  SetToUpperCase(output);
  return output;
}
static MTL_INLINE std::string ToUpperCase(const std::string& input)
{
  std::string output = input;
  SetToUpperCase(output);
  return output;
}

static MTL_INLINE std::string ToUTF8(const std::wstring& str)
{
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.to_bytes(str);
}
static MTL_INLINE std::wstring ToUTF16(const std::string& str)
{
#ifdef WIN32
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
#else  // g++ 5.4
  std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
#endif
  return converter.from_bytes(str);
}

}  // namespace MTL


#endif  // MTL_STRING_HELPERS_H
