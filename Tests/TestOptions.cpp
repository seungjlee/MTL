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

#include <MTL/Tools/Test.h>
#include <MTL/Options.h>

using namespace MTL;

// Declare options at file scope so the MTL_OPTION macro and Registrar are exercised.
MTL_OPTION(bool,        true,    OptVerbose,    "Verbose flag.")
MTL_OPTION(int,         42,      OptIters,      "Iteration count.")
MTL_OPTION(double,      1.5,     OptScale,      "Scale factor.")
MTL_OPTION(std::string, "alpha", OptName,       "Name string.")
MTL_OPTION(String,      L"beta", OptWideName,   "Wide name string.")

TEST(Test_FormatAndParseValues)
{
  MTL_VERIFY(Options::FormatOptionValue(true)  == L"true");
  MTL_VERIFY(Options::FormatOptionValue(false) == L"false");
  MTL_VERIFY(Options::FormatOptionValue(7)     == L"7");
  MTL_VERIFY(Options::FormatOptionValue(String(L"hi")) == L"hi");
  MTL_VERIFY(Options::FormatOptionValue(std::string("hi")) == L"hi");

  bool b = false;
  MTL_VERIFY( Options::ParseOptionValue(L"true", b)  && b == true);
  MTL_VERIFY( Options::ParseOptionValue(L"YES", b)   && b == true);
  MTL_VERIFY( Options::ParseOptionValue(L"On", b)    && b == true);
  MTL_VERIFY( Options::ParseOptionValue(L"1", b)     && b == true);
  MTL_VERIFY( Options::ParseOptionValue(L"", b)      && b == true);
  MTL_VERIFY( Options::ParseOptionValue(L"false", b) && b == false);
  MTL_VERIFY( Options::ParseOptionValue(L"NO", b)    && b == false);
  MTL_VERIFY( Options::ParseOptionValue(L"off", b)   && b == false);
  MTL_VERIFY( Options::ParseOptionValue(L"0", b)     && b == false);
  MTL_VERIFY(!Options::ParseOptionValue(L"maybe", b));

  int i = 0;
  MTL_VERIFY(Options::ParseOptionValue(L"123", i) && i == 123);
  MTL_VERIFY(!Options::ParseOptionValue(L"abc", i));

  double d = 0;
  MTL_VERIFY(Options::ParseOptionValue(L"3.25", d) && d == 3.25);

  String ws;
  MTL_VERIFY(Options::ParseOptionValue(L"hello", ws) && ws == L"hello");
  std::string ns;
  MTL_VERIFY(Options::ParseOptionValue(L"world", ns) && ns == "world");
}

TEST(Test_OptionAccessors)
{
  MTL_VERIFY(FLAGS_OptIters.Value() == 42);
  MTL_VERIFY(FLAGS_OptScale.DefaultText() == L"1.5");
  MTL_VERIFY(FLAGS_OptName.CurrentText() == L"alpha");
  MTL_VERIFY(FLAGS_OptVerbose.IsFlag() == true);
  MTL_VERIFY(FLAGS_OptIters.IsFlag() == false);

  // Implicit conversion via operator const T&().
  int copy = FLAGS_OptIters;
  MTL_VERIFY(copy == 42);

  // Mutable access.
  FLAGS_OptIters.Value() = 7;
  MTL_VERIFY(FLAGS_OptIters.Value() == 7);
  FLAGS_OptIters.Value() = 42;
}

TEST(Test_ParseArgs_HappyPaths)
{
  // Inline =value form, mixed-case names, and bare bool flag.
  std::vector<String> args = {
    L"--OptIters=10",
    L"--optscale", L"2.5",
    L"--OptVerbose",
    L"--OptName=hello",
    L"positional1",
    L"positional2",
  };
  std::vector<String> positional = Options::ParseArgs(args);
  MTL_EQUAL(positional.size(), size_t(2));
  MTL_VERIFY(positional[0] == L"positional1");
  MTL_VERIFY(positional[1] == L"positional2");
  MTL_EQUAL(FLAGS_OptIters.Value(),     10);
  MTL_EQUAL(FLAGS_OptScale.Value(),     2.5);
  MTL_EQUAL(FLAGS_OptVerbose.Value(),   true);
  MTL_VERIFY(FLAGS_OptName.Value() == "hello");

  // Bare-dash and double-dash both accepted; explicit bool=false.
  std::vector<String> args2 = { L"-OptVerbose=false", L"--OptIters=99" };
  Options::ParseArgs(args2);
  MTL_EQUAL(FLAGS_OptVerbose.Value(), false);
  MTL_EQUAL(FLAGS_OptIters.Value(),   99);

  // Wide string option round-trip.
  std::vector<String> args3 = { L"--OptWideName=gamma" };
  Options::ParseArgs(args3);
  MTL_VERIFY(FLAGS_OptWideName.Value() == L"gamma");
}

TEST(Test_ParseArgs_AllowUnknown)
{
  std::vector<String> args = {
    L"--unknown=value",
    L"--OptIters=5",
    L"-x",
    L"plain",
  };
  std::vector<String> positional = Options::ParseArgs(args, /*allowUnknown=*/true);
  MTL_EQUAL(FLAGS_OptIters.Value(), 5);
  MTL_EQUAL(positional.size(), size_t(3));
  MTL_VERIFY(positional[0] == L"--unknown=value");
  MTL_VERIFY(positional[1] == L"-x");
  MTL_VERIFY(positional[2] == L"plain");
}

TEST(Test_ParseArgs_NonOptionPositional)
{
  // A token shorter than 2 chars or not starting with '-' is positional.
  std::vector<String> args = { L"a", L"-", L"b" };
  std::vector<String> positional = Options::ParseArgs(args);
  MTL_EQUAL(positional.size(), size_t(3));
}

TEST(Test_ParseFromCharArgv)
{
  const char* argv[] = { "prog", "--OptIters=21", "leftover" };
  std::vector<String> positional = Options::Parse(3, const_cast<char**>(argv));
  MTL_EQUAL(FLAGS_OptIters.Value(), 21);
  MTL_EQUAL(positional.size(), size_t(1));
  MTL_VERIFY(positional[0] == L"leftover");
}

TEST(Test_ParseFromWcharArgv)
{
  const wchar_t* argv[] = { L"prog", L"--OptIters=33", L"end" };
  std::vector<String> positional = Options::Parse(3, const_cast<wchar_t**>(argv));
  MTL_EQUAL(FLAGS_OptIters.Value(), 33);
  MTL_EQUAL(positional.size(), size_t(1));
}

TEST(Test_ParseFromContainer)
{
  std::vector<String> argv = { L"progName", L"--OptIters=44", L"trailing" };
  std::vector<String> positional = Options::ParseFromContainer(argv);
  MTL_EQUAL(FLAGS_OptIters.Value(), 44);
  MTL_EQUAL(positional.size(), size_t(1));
  MTL_VERIFY(positional[0] == L"trailing");

  // Without skipProgramName: the first token is treated as positional too.
  std::vector<String> argv2 = { L"--OptIters=55", L"x" };
  std::vector<String> p2 = Options::ParseFromContainer(argv2, /*skipProgramName=*/false);
  MTL_EQUAL(FLAGS_OptIters.Value(), 55);
  MTL_EQUAL(p2.size(), size_t(1));
}

TEST(Test_PrintHelp)
{
  // PrintHelp writes to ConsoleOut; redirect by capturing into Out() instead.
  // Just verify it does not crash and walks the registry.
  Options::PrintHelp(Out());
}
