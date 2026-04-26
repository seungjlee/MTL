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

// Negative tests for the test framework itself: drive every comparison-helper
// failure branch and the exception handlers in Test::Run_All. We provide our
// own main so an "expected" failure count means success and the executable
// still exits with status 0.
#define MTL_TEST_NO_MAIN
#include <MTL/Tools/Test.h>
#include <MTL/Exception.h>

using namespace MTL;

// Each test below intentionally records exactly one expected failure. Update
// this constant whenever you add or remove an expected-failure assertion.
// Counts the negative TEST blocks below plus the Initialize/Shutdown throws
// and the catch(...) handler in Run_All.
static const U64 kExpectedFailures = 16;

TEST(Negative_VerifyFalse)
{
  Test::Verify(false, "intentional verify(false)", L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_EqualMismatchTemplate)
{
  // Mismatched ints exercise the templated Equal failure branch.
  Test::Equal(int(7), int(11), L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_EqualMismatchString)
{
  Test::Equal(std::string("aaa"), std::string("bbb"), L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_EqualFloatTooFar)
{
  Test::EqualFloat(1.0, 2.0, 0.01, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_EqualFloatNan)
{
  // NaN compares unequal to itself; this exercises the actual != actual branch.
  double nan = std::numeric_limits<double>::quiet_NaN();
  Test::EqualFloat(nan, 0.0, 1e-9, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_LessThanEqualLimit)
{
  // actual == limit: exercises the "equal to limit" branch of LessThan.
  Test::LessThan(5, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_LessThanGreaterThanLimit)
{
  // actual > limit: exercises the strictly-greater branch of LessThan.
  Test::LessThan(7, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_LessThanOrEqualToGreater)
{
  Test::LessThanOrEqualTo(7, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_GreaterThanEqualLimit)
{
  Test::GreaterThan(5, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_GreaterThanLessLimit)
{
  Test::GreaterThan(3, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_GreaterThanOrEqualToLess)
{
  Test::GreaterThanOrEqualTo(3, 5, L"TestFrameworkNegative.cpp", __LINE__);
}

TEST(Negative_ThrowsMTLException)
{
  // Exercises the catch (const Exception&) handler in Run_All.
  throw Exception(L"Intentional MTL::Exception from negative test.");
}

TEST(Negative_ThrowsStdException)
{
  // Exercises the catch (const std::exception&) handler in Run_All.
  throw std::runtime_error("Intentional std::exception from negative test.");
}

TEST(Negative_ThrowsUnknown)
{
  // Exercises the catch (...) handler in Run_All.
  throw 42;
}

// Initialize throws an MTL::Exception so the catch(Exception) branch in
// Test::Initialize executes. Shutdown throws std::runtime_error so the
// catch(std::exception) branch in Test::Shutdown executes. Each throw counts
// as one expected failure.
TEST_INITIALIZE()
{
  throw Exception(L"Intentional Initialize throw from negative test.");
}

TEST_SHUTDOWN()
{
  throw std::runtime_error("Intentional Shutdown throw from negative test.");
}

int main(int argc, char* argv[])
{
  // Replicate the argv-to-Arguments wiring that Test.h's default main does.
  DynamicVector<String>& arguments =
    const_cast<DynamicVector<String>&>(Test::Arguments());
  for (int i = 0; i < argc; i++)
    arguments.PushBack(ToUTF16(argv[i]));

  // Drives the -DisableColorRGB24 fallback branch in Run_All.
  arguments.PushBack(String(L"-DisableColorRGB24"));

  Test::RunAll();

  U64 actual = Test::TotalNumberOfFailures();
  if (actual != kExpectedFailures)
  {
    ConsoleOut << COLOR_RED << L"TestFrameworkNegative: expected " << kExpectedFailures
               << L" failures, observed " << actual << COLOR_RESET << std::endl;
    return 1;
  }

  ConsoleOut << COLOR_LGREEN << L"TestFrameworkNegative: all " << actual
             << L" expected failure paths exercised." << COLOR_RESET << std::endl;
  return 0;
}
