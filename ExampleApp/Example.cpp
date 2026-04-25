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

#include <MTL/Options.h>
#include <MTL/Tools/Test.h>

MTL_OPTION(bool,        true,           verbose,    "Print extra status output.")
MTL_OPTION(int,         10,             iterations, "Number of iterations to run.")
MTL_OPTION(double,      1.5,            scale,      "Scale factor applied to the result.")
MTL_OPTION(std::string, "world",        greet,      "Name to greet.")

MTL_APP()
{
  // Test::Arguments() already contains argv[0] (the program name) followed by
  // the rest, so skip the first element when parsing.
  MTL::Options::ParseFromContainer(MTL::Test::Arguments());

  if (FLAGS_verbose)
  {
    ConsoleOut << L"verbose    = " << MTL::Options::FormatOptionValue(FLAGS_verbose.Value())    << std::endl;
    ConsoleOut << L"iterations = " << FLAGS_iterations.Value() << std::endl;
    ConsoleOut << L"scale      = " << FLAGS_scale.Value()      << std::endl;
    ConsoleOut << L"greet      = " << MTL::ToUTF16(FLAGS_greet.Value()) << std::endl;
  }

  double accumulator = 0.0;
  for (int i = 0; i < FLAGS_iterations; i++)
    accumulator += FLAGS_scale * (i + 1);

  ConsoleOut << L"Hello, " << MTL::ToUTF16(FLAGS_greet.Value()) << L"!" << std::endl;
  ConsoleOut << L"Sum after " << FLAGS_iterations.Value()
             << L" iterations: " << accumulator << std::endl;
}
