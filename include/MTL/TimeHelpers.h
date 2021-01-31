//
// Math Template Library
//
// Copyright (c) 2021: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_TIME_HELPERS_H
#define MTL_TIME_HELPERS_H

#include <chrono>
#include <iomanip>
#include <sstream>
#include <MTL/StringHelpers.h>

namespace MTL
{

template<class CharT> inline std::basic_string<CharT> GetTimeWithUnit(double time, const char* unit);
template<> inline std::basic_string<char> GetTimeWithUnit(double seconds, const char* unit)
{
  return StringPrint<char>("%.2f %s", seconds, unit);
}
template<> inline std::basic_string<wchar_t> GetTimeWithUnit(double seconds, const char* unit)
{
  return StringPrint<wchar_t>(L"%.2f %s", seconds, unit);
}

template<class CharT>
inline std::basic_string<CharT> GetTimeString(double seconds)
{
  if (seconds < 1e-3)
    return GetTimeWithUnit<CharT>(seconds * 1e6, "microseconds");

  if (seconds < 1)
    return GetTimeWithUnit<CharT>(seconds * 1e3, "milliseconds");

  return GetTimeWithUnit<CharT>(seconds, "seconds");
}
inline std::wstring GetTime(double seconds)
{
  return GetTimeString<wchar_t>(seconds);
}
inline std::string GetTimeUTF8(double seconds)
{
  return GetTimeString<char>(seconds);
}

#ifndef MTL_DATE_TIME_FORMAT
#define MTL_DATE_TIME_FORMAT "%F %T"
#endif

template<class CharT> inline const CharT* DateTimeFormat();
template<> inline const wchar_t* DateTimeFormat()
{
  return TO_WCHAR(MTL_DATE_TIME_FORMAT);
}
template<> inline const char* DateTimeFormat()
{
  return MTL_DATE_TIME_FORMAT;
}

template<class CharT>
inline std::basic_string<CharT> CurrentDateTime()
{
  std::basic_stringstream<CharT> stream;

  auto now = std::chrono::system_clock::now();
  std::chrono::milliseconds msecs =
    std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

  time_t currentTime = std::chrono::system_clock::to_time_t(now);

#ifdef WIN32
  tm tbuffer;
  localtime_s(&tbuffer, &currentTime);
  stream << std::put_time(&tbuffer, DateTimeFormat<CharT>());
#else
  stream << std::put_time(localtime(&currentTime), DateTimeFormat<CharT>());
#endif
  stream << StringPrintf(".%03lld", msecs.count()).c_str();

  return stream.str();
}
inline std::wstring GetCurrentDateTime()
{
  return CurrentDateTime<wchar_t>();
}
inline std::string GetCurrentDateTimeUTF8()
{
  return CurrentDateTime<char>();
}

}  // namespace MTL


#endif  // MTL_STRING_HELPERS_H
