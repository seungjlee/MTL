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
#include <codecvt>
#include <cstring>
#include <locale>
#include <wctype.h>

namespace MTL
{

// Call this to switch between std::cout (printf) and std::wcout (wprintf).
inline void ResetOutputStream()
{
#ifdef WIN32
  FILE* f;
  freopen_s(&f, nullptr, "w", stdout);
#else
  freopen(nullptr, "w", stdout);
#endif
}

inline char ToLower(char input)
{
  return (char)tolower(input);
}
inline char ToUpper(char input)
{
  return (char)toupper(input);
}
inline void SetToLowerCase(std::wstring& str)
{
  std::transform(str.begin(), str.end(), str.begin(), towlower);
}
inline void SetToLowerCase(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), ToLower);
}
inline void SetToUpperCase(std::wstring& str)
{
  std::transform(str.begin(), str.end(), str.begin(), towupper);
}
inline void SetToUpperCase(std::string& str)
{
  std::transform(str.begin(), str.end(), str.begin(), ToUpper);
}

inline std::wstring ToLowerCase(const std::wstring& input)
{
  std::wstring output = input;
  SetToLowerCase(output);
  return output;
}
inline std::string ToLowerCase(const std::string& input)
{
  std::string output = input;
  SetToLowerCase(output);
  return output;
}
inline std::wstring ToUpperCase(const std::wstring& input)
{
  std::wstring output = input;
  SetToUpperCase(output);
  return output;
}
inline std::string ToUpperCase(const std::string& input)
{
  std::string output = input;
  SetToUpperCase(output);
  return output;
}

inline std::string ToUTF8(const std::wstring& str)
{
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
  return converter.to_bytes(str);
}
inline std::wstring ToUTF16(const std::string& str)
{
#ifdef WIN32
  std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
#else  // g++ 5.4
  std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
#endif
  return converter.from_bytes(str);
}

template <class CharType>
inline std::basic_string<CharType> ToString(const char* ptr);
template <>
inline std::basic_string<char> ToString(const char* ptr)
{
  return std::basic_string<char>(ptr);
}
inline std::basic_string<wchar_t> ToString(const char* ptr)
{
  return ToUTF16(ToString<char>(ptr));
}

template <class ... Args>
static std::string StringPrintf(const char* format, Args ... args)
{
  size_t formatLength = strlen(format);

  if (formatLength == 0)
    return "";

  std::string buffer;  // Use it as container.

  const std::size_t N = sizeof...(Args);  // Number of arguments.
  buffer.resize(formatLength + N * 16 + 16);
  uint64_t bytes = snprintf(&buffer[0], buffer.size(), format, args ...);

  if (bytes >= buffer.size())
  {
    buffer.resize(bytes + 1);
    snprintf(&buffer[0], buffer.size(), format, args ...);
  }

  return buffer.data();  // Make it resize properly.
}
template <class ... Args>
inline std::wstring StringPrintf(const wchar_t* format, Args ... args)
{
  return ToUTF16(StringPrintf(ToUTF8(format).c_str(), args ...));
}

template <class CharType, class ... Args>
inline std::basic_string<CharType> StringPrint(const CharType* format, Args ... args)
{
  return StringPrintf(format, args ...);
}

}  // namespace MTL


#endif  // MTL_STRING_HELPERS_H
