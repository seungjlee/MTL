//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee
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


#ifndef MTL_TEST_H
#define MTL_TEST_H

#include "CPU.h"
#include "Exception.h"
#include "Timer.h"
#include "DynamicVector.h"

//
// Macros.
//
#define TEST(TestName)                                                                             \
class Test_ ## TestName : public Test                                                              \
{                                                                                                  \
public:                                                                                            \
  Test_ ## TestName(const String& testName)                                                        \
    : Test(testName)  { FilePath_ = MTL__FILE__; }                                                 \
  virtual void Run();                                                                              \
};                                                                                                 \
static Test_ ## TestName                                                                           \
 Test_ ## TestName ## _Instance_(String(TOWCHAR(#TestName)));                                      \
void Test_ ## TestName::Run()

#define MTL_VERIFY(Expression)  MTL::Test::Verify(Expression, #Expression, MTL__FILE__, __LINE__)
#define MTL_EQUAL(Actual, Expected)  MTL::Test::Equal(Actual, Expected, MTL__FILE__, __LINE__)
#define MTL_EQUAL_FLOAT(Actual, Expected, Tolerance) \
MTL::Test::EqualFloat(double(Actual), double(Expected), double(Tolerance), MTL__FILE__, __LINE__)

namespace MTL
{

class Test
{
public:
  Test(const String& name)
    : Name_(name), NumberOfFailures_(0), TimeElapsed_(0)
  { List_.PushBack(this); }

  virtual ~Test()  {}

  virtual void Run() = 0;

  static void RunAll()
  {
    Out() << std::endl;
    for (U32 i = 0; i < List_.Size(); i++)
    {
      Out() << "[" << List_[i]->Name_ << "]" << " Begins..." << std::endl;

      try
      {
        Timer timer(true);
        List_[i]->Run();
        timer.Stop();
        List_[i]->TimeElapsed_ = timer.Seconds();
      }
      catch(const Exception& e)
      {
        List_[i]->NumberOfFailures_++;
        Out() << "ERROR: " << e.Message() << std::endl;
      }
      catch(const std::exception& e)
      {
        List_[i]->NumberOfFailures_++;
        Out() << "std::exception: " << e.what() << std::endl;
      }
      catch(...)
      {
        List_[i]->NumberOfFailures_++;
        Out() << "UNEXPECTED ERROR!" << std::endl;
      }

      Out() << "[" << List_[i]->Name_ << "]" << " Ends.";
      if (List_[i]->NumberOfFailures_ > 0)
        Out() << "  FAILED!";

      Out() << "  Number of errors: " << List_[i]->NumberOfFailures_
            << ", Time: " << float(List_[i]->TimeElapsed_) << " seconds."
            << std::endl << std::endl;

      TotalNumberOfFailures_ += List_[i]->NumberOfFailures_;
      TotalTimeElapsed_ += List_[i]->TimeElapsed_;
    }

    Out() << "Total number of errors: " << TotalNumberOfFailures() << std::endl;
    Out() << "Total time: " << TotalTimeElapsed_ << " seconds." << std::endl;
  }

  static OutputStream& Out()
  {
    return Output_;
  }

  static const DynamicVector<String>& Arguments()  { return Arguments_; }

  static U64 TotalNumberOfFailures()  { return TotalNumberOfFailures_; }

protected:
  U64 NumberOfFailures_;
  double TimeElapsed_;  // In seconds.

  static String FilePath_;

  void Verify(bool e, char* expression, wchar_t* file, U64 line)
  {
    if (!e)
    {
      NumberOfFailures_++;
      Out() << "[In file '" << file << "' - line " << line  << "]" << std::endl
            << "  '" << expression << "' failed!" << std::endl;
    }
  }

  template <class T1, class T2>
  void Equal(const T1& actual, const T2& expected, wchar_t* file, U64 line)
  {
    if (actual != expected)
    {
      NumberOfFailures_++;
      Out() << "[File '" << file << "' - line " << line << "]" << std::endl
            << "  Actual value is '" << actual << "' but '"
            << expected << "' is expected!'" << std::endl;
    }
  }

  void EqualFloat(double actual, double expected, double tolerance, wchar_t* file, U64 line)
  {
    double difference = actual - expected;
    if (MTL::Abs(difference) > tolerance || actual != actual)
    {
      NumberOfFailures_++;
      Out() << "[File '" << file << "' - line " << line  << "]" << std::endl
            << "  Actual value is " << actual << " but " << expected
            << " is expected; difference is: " << difference << std::endl;
    }
  }

private:
  String Name_;

  static OutputStream& Output_;
  static DynamicVector<Test*> List_;
  static DynamicVector<String> Arguments_;

  static U64 TotalNumberOfFailures_;
  static double TotalTimeElapsed_;  // In seconds.
};

String Test::FilePath_;
OutputStream& Test::Output_ = std::wcout;
DynamicVector<Test*> Test::List_;
DynamicVector<String> Test::Arguments_;
U64 Test::TotalNumberOfFailures_;
double Test::TotalTimeElapsed_;

}  // namespace MTL


#if defined(WIN32) || defined(WIN64)
  // For debugging memory leaks with Visual C++.
  #define _CRTDBG_MAP_ALLOC
  #include <crtdbg.h>

  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#endif  // WIN32

int wmain(int argc, wchar_t* argv[])
{
#if defined(WIN32) || defined(WIN64)

#ifndef _DEBUG
  // Avoid pop up windows when not debugging.
  DWORD dwMode = SetErrorMode(SEM_NOGPFAULTERRORBOX);
  SetErrorMode(dwMode | SEM_NOGPFAULTERRORBOX);
#endif

  // Report memory leaks with Visual C++.
  _CrtSetReportMode(_CRT_WARN,   _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG);

  _CrtSetReportMode
    (_CRT_ERROR,  _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG | _CRTDBG_MODE_WNDW);
  _CrtSetReportMode
    (_CRT_ASSERT, _CRTDBG_MODE_FILE | _CRTDBG_MODE_DEBUG | _CRTDBG_MODE_WNDW);

  _CrtSetReportFile(_CRT_WARN,   _CRTDBG_FILE_STDOUT);
  _CrtSetReportFile(_CRT_ERROR,  _CRTDBG_FILE_STDOUT);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);

  _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

  try
  {
    std::wcout << std::endl << "Number Of Actual Cores: "
               << MTL::CPU::Instance().NumberOfCores() << std::endl;
    std::wcout << "Initial Number Of OpenMP Threads: "
               << MTL::CPU::Instance().NumberOfThreads() << std::endl;

    MTL::DynamicVector<MTL::String>& arguments =
      const_cast<MTL::DynamicVector<MTL::String>&>(MTL::Test::Arguments());

    for (int i = 0; i < argc; i++)
      arguments.PushBack(argv[i]);

    MTL::Test::RunAll();
  }
  catch(const MTL::Exception& e)
  {
    std::wcout << "ERROR: " << e.Message() << std::endl;
  }
  catch(...)
  {
    std::wcout << "UNEXPECTED ERROR!" << std::endl;
  }

  return (int)MTL::Test::TotalNumberOfFailures();
}

#endif  // MTL_TEST_H
