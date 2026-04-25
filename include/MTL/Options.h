//
// Math Template Library
//
// Copyright (c) 2026: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_OPTIONS_H
#define MTL_OPTIONS_H

#include <MTL/Definitions.h>
#include <MTL/StringHelpers.h>

#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#ifndef MTL_OPTIONS_CASE_INSENSITIVE
#define MTL_OPTIONS_CASE_INSENSITIVE 1
#endif

namespace MTL
{
namespace Options
{

// Convert a value to its display string. wstringstream handles all the common
// types; specialize bool and string types for nicer output.
template <class T>
inline String FormatOptionValue(const T& value)
{
  std::wstringstream ss;
  ss << value;
  return ss.str();
}
inline String FormatOptionValue(bool value)        { return value ? L"true" : L"false"; }
inline String FormatOptionValue(const String& s)   { return s; }
inline String FormatOptionValue(const std::string& s) { return ToUTF16(s); }

// Parse a value from its command-line string form. Returns false on failure.
template <class T>
inline bool ParseOptionValue(const String& text, T& out)
{
  std::wstringstream ss(text);
  ss >> out;
  return !ss.fail() && ss.eof();
}
inline bool ParseOptionValue(const String& text, bool& out)
{
  String lower = ToLowerCase(text);
  if (lower == L"true" || lower == L"1" || lower == L"yes" || lower == L"on" || lower.empty())
  {
    out = true;
    return true;
  }
  if (lower == L"false" || lower == L"0" || lower == L"no" || lower == L"off")
  {
    out = false;
    return true;
  }
  return false;
}
inline bool ParseOptionValue(const String& text, String& out)
{
  out = text;
  return true;
}
inline bool ParseOptionValue(const String& text, std::string& out)
{
  out = ToUTF8(text);
  return true;
}

class OptionBase
{
public:
  OptionBase(const String& name, const String& typeName, const String& help)
    : Name_(name), TypeName_(typeName), Help_(help) {}
  virtual ~OptionBase() = default;

  virtual bool Parse(const String& text) = 0;
  virtual bool IsFlag() const = 0;  // True for bool options that may appear without a value.
  virtual String DefaultText() const = 0;
  virtual String CurrentText() const = 0;

  const String& Name()     const { return Name_;     }
  const String& TypeName() const { return TypeName_; }
  const String& Help()     const { return Help_;     }

private:
  String Name_;
  String TypeName_;
  String Help_;
};

template <class T>
class Option : public OptionBase
{
public:
  Option(const String& name, const String& typeName, const T& defaultValue, const String& help)
    : OptionBase(name, typeName, help), Value_(defaultValue), Default_(defaultValue) {}

  bool Parse(const String& text)  override { return ParseOptionValue(text, Value_); }
  bool IsFlag() const             override { return std::is_same<T, bool>::value; }
  String DefaultText() const      override { return FormatOptionValue(Default_); }
  String CurrentText() const      override { return FormatOptionValue(Value_); }

  operator const T&() const { return Value_; }
  const T& Value() const    { return Value_; }
  T&       Value()          { return Value_; }

private:
  T Value_;
  T Default_;
};

inline String NormalizeName(const String& name)
{
#if MTL_OPTIONS_CASE_INSENSITIVE
  return ToLowerCase(name);
#else
  return name;
#endif
}

inline std::map<String, OptionBase*>& Registry()
{
  static std::map<String, OptionBase*> registry;
  return registry;
}

inline void Register(OptionBase* option)
{
  Registry()[NormalizeName(option->Name())] = option;
}

template <class T>
struct Registrar
{
  explicit Registrar(Option<T>* option) { Register(option); }
};

inline void PrintHelp(OutputStream& out = ConsoleOut)
{
  out << L"Options:" << std::endl;
  for (const auto& kv : Registry())
  {
    const OptionBase* option = kv.second;
    out << L"  --" << option->Name()
        << L"=<" << option->TypeName() << L">"
        << L"   [default: " << option->DefaultText() << L"]" << std::endl
        << L"      " << option->Help() << std::endl;
  }
}

namespace detail
{
inline String StripLeadingDashes(const String& token)
{
  size_t start = 0;
  while (start < token.size() && start < 2 && token[start] == L'-')
    start++;
  return token.substr(start);
}

inline void ReportAndExit(const String& message, int code)
{
  ConsoleOut << message << std::endl;
  PrintHelp();
  std::exit(code);
}
}  // namespace detail

// Parses options from `args` (which must NOT include the program name).
// Returns the positional (non-option) arguments. On --help, prints help and
// exits with code 0. On unknown or malformed options, prints help and exits
// with code 1 unless `allowUnknown` is true, in which case unknown tokens are
// passed through as positional.
inline std::vector<String> ParseArgs(const std::vector<String>& args, bool allowUnknown = false)
{
  std::vector<String> positional;
  for (size_t i = 0; i < args.size(); i++)
  {
    const String& token = args[i];
    if (token.size() < 2 || token[0] != L'-')
    {
      positional.push_back(token);
      continue;
    }

    String body = detail::StripLeadingDashes(token);
    String name, value;
    bool hasInlineValue = false;

    size_t eq = body.find(L'=');
    if (eq != String::npos)
    {
      name = body.substr(0, eq);
      value = body.substr(eq + 1);
      hasInlineValue = true;
    }
    else
    {
      name = body;
    }

    String normalized = NormalizeName(name);
    if (normalized == L"help" || normalized == L"h")
    {
      PrintHelp();
      std::exit(0);
    }

    auto found = Registry().find(normalized);
    if (found == Registry().end())
    {
      if (allowUnknown)
      {
        positional.push_back(token);
        continue;
      }
      detail::ReportAndExit(L"Unknown option: --" + name, 1);
    }

    OptionBase* option = found->second;
    if (!hasInlineValue)
    {
      if (option->IsFlag())
      {
        value = L"true";
      }
      else if (i + 1 < args.size())
      {
        value = args[++i];
      }
      else
      {
        detail::ReportAndExit(L"Missing value for --" + name, 1);
      }
    }

    if (!option->Parse(value))
      detail::ReportAndExit(L"Invalid value for --" + name + L": '" + value + L"'", 1);
  }
  return positional;
}

inline std::vector<String> Parse(int argc, char** argv, bool allowUnknown = false)
{
  std::vector<String> args;
  args.reserve(argc > 1 ? argc - 1 : 0);
  for (int i = 1; i < argc; i++)
    args.push_back(ToUTF16(argv[i]));
  return ParseArgs(args, allowUnknown);
}

inline std::vector<String> Parse(int argc, wchar_t** argv, bool allowUnknown = false)
{
  std::vector<String> args;
  args.reserve(argc > 1 ? argc - 1 : 0);
  for (int i = 1; i < argc; i++)
    args.push_back(argv[i]);
  return ParseArgs(args, allowUnknown);
}

// Convenience overload for any container with begin()/end() yielding String,
// e.g. MTL::Test::Arguments() (a DynamicVector<String> that already includes
// the program name as element 0).
template <class Container>
inline std::vector<String> ParseFromContainer(const Container& argv, bool skipProgramName = true,
                                              bool allowUnknown = false)
{
  std::vector<String> args;
  bool skipped = !skipProgramName;
  for (const String& a : argv)
  {
    if (!skipped) { skipped = true; continue; }
    args.push_back(a);
  }
  return ParseArgs(args, allowUnknown);
}

}  // namespace Options
}  // namespace MTL


// MTL_OPTION(Type, Default, Name, "Help text") declares a global option named
// FLAGS_<Name> of the given Type, registers it, and exposes it on the command
// line as --<Name>=<value>. Names are matched case-insensitively by default
// (override via -DMTL_OPTIONS_CASE_INSENSITIVE=0). For bool options, --<Name>
// alone is equivalent to --<Name>=true.
#define MTL_OPTION(Type, Default, Name, Help)                                                      \
  inline ::MTL::Options::Option<Type> FLAGS_##Name {                                               \
    L"" #Name, L"" #Type, (Default), ::MTL::ToUTF16(std::string(Help)) };                          \
  inline ::MTL::Options::Registrar<Type> MTL_OptionRegistrar_##Name { &FLAGS_##Name };

#endif  // MTL_OPTIONS_H
