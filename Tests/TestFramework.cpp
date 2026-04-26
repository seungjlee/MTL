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

using namespace MTL;

// Exercise the comparison helpers exposed via the MTL_* convenience macros.
// All inputs are chosen to satisfy the assertion so the harness reports zero
// errors; the goal is to make the helper templates execute at least once for
// representative integral and floating-point types.
TEST(Test_ComparisonHelpers)
{
  MTL_LESS_THAN(1, 5);
  MTL_LESS_THAN(0.5, 1.0);

  MTL_LESS_THAN_OR_EQUAL_TO(5, 5);
  MTL_LESS_THAN_OR_EQUAL_TO(2, 5);
  MTL_LESS_THAN_OR_EQUAL_TO(1.0, 1.0);

  MTL_GREATER_THAN(7, 3);
  MTL_GREATER_THAN(2.5, 1.5);

  MTL_GREATER_THAN_OR_EQUAL_TO(3, 3);
  MTL_GREATER_THAN_OR_EQUAL_TO(8, 3);

  MTL_EQUAL(std::string("hi"), std::string("hi"));

  // ShowProgressBar overloads (wide-string and narrow-string).
  ShowProgressBar(0.5, L"halfway");
  ShowProgressBar(1.0, std::string("done"));
}

TEST(Test_ArgumentLookup)
{
  // Argument 0 is the test executable path. Look it up both ways.
  const DynamicVector<String>& args = Test::Arguments();
  MTL_VERIFY(args.Size() >= 1u);

  // Case-sensitive find of a known argument: argv[0] (the program path).
  const String& argv0 = args[0];
  I32 idx = Test::FindArgument(argv0);
  MTL_EQUAL(idx, I32(0));

  // Case-insensitive find of the same argument.
  idx = Test::FindArgumentIgnoreCase(argv0);
  MTL_EQUAL(idx, I32(0));

  // Look for an argument that doesn't exist.
  MTL_EQUAL(Test::FindArgument(L"--definitely-no-such-flag-xyz"), I32(-1));
  MTL_EQUAL(Test::FindArgumentIgnoreCase(L"--no-such-flag"),       I32(-1));

  // TestPath / TestFilePathName both look at the source-file path baked in
  // by the TEST() macro; just walk the path-stripping loop.
  String filePath = Test::TestFilePathName();
  String dirPath  = Test::TestPath();
  MTL_VERIFY(filePath.size() >= dirPath.size());
}

// Throw from inside a TEST() to drive the exception handler in Run_All.
// The framework catches and counts as failure, but the per-test failure count
// is reported only at end of run; this test itself is reported as failing,
// so we can't include it in the normal harness without polluting the result.
// Instead we directly call Test::Verify with a passing condition (the EQUAL
// helper for std::string is the lone Equal overload not yet covered).
TEST(Test_VerifyDirectAPI)
{
  Test::Verify(true, "1 == 1", L"TestFramework.cpp", 1);
  Test::EqualFloat(1.0, 1.0, 0.0, L"TestFramework.cpp", 2);
}
