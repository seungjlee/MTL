//
// Math Template Library
//
// Copyright (c) 2014-2019: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#include <MTL/CPU.h>
#include <MTL/Colors.h>
#include <MTL/DynamicVector.h>
#include <MTL/Exception.h>
#include <MTL/Timer.h>
#include <MTL/ProgressBar.h>

//
// Macros.
//
#ifndef COLOR_ERROR
#define COLOR_ERROR COLOR_FG(255,55,55)
#endif

#define TEST_MACRO(TestName, Enabled)                                                              \
class Test_ ## TestName : public MTL::Test                                                         \
{                                                                                                  \
public:                                                                                            \
  Test_ ## TestName(const MTL::String& testName)                                                   \
    : Test(testName, Enabled)  { FilePath_ = MTL__FILE__; }                                        \
  virtual void Run();                                                                              \
};                                                                                                 \
static Test_ ## TestName                                                                           \
 Test_ ## TestName ## _Instance_(MTL::String(TOWCHAR(#TestName)));                                 \
void Test_ ## TestName::Run()

#define TEST(TestName)           TEST_MACRO(TestName, true)
#define DISABLED_TEST(TestName)  TEST_MACRO(TestName, false)

#define TEST_INITIALIZE()                                                                          \
static void MTL_TestInitialize_();                                                                 \
static int MTL_Test_Initialze_ = MTL::Test::SetInitializeFunction(MTL_TestInitialize_);            \
static void MTL_TestInitialize_()

#define TEST_SHUTDOWN()                                                                            \
static void MTL_TestShutdown_();                                                                   \
static int MTL_Test_Shutdown_ = MTL::Test::SetShutdownFunction(MTL_TestShutdown_);                 \
static void MTL_TestShutdown_()

#define MTL_VERIFY(Expression)                                                                     \
MTL::Test::Verify(Expression, #Expression, MTL::String(MTL__FILE__), __LINE__)

#define MTL_EQUAL(Actual, Expected)                                                                \
MTL::Test::Equal(Actual, Expected, MTL::String(MTL__FILE__), __LINE__)

#define MTL_EQUAL_FLOAT(Actual, Expected, Tolerance)                                               \
MTL::Test::EqualFloat(double(Actual), double(Expected), double(Tolerance),                         \
                      MTL::String(MTL__FILE__), __LINE__)

#define MTL_LESS_THAN(Actual, Limit)                                                               \
MTL::Test::LessThan(Actual, Limit, MTL::String(MTL__FILE__), __LINE__)

#define MTL_LESS_THAN_OR_EQUAL_TO(Actual, Limit)                                                   \
MTL::Test::LessThanOrEqualTo(Actual, Limit, MTL::String(MTL__FILE__), __LINE__)

#define MTL_GREATER_THAN(Actual, Limit)                                                            \
MTL::Test::GreaterThan(Actual, Limit, MTL::String(MTL__FILE__), __LINE__)

#define MTL_GREATER_THAN_OR_EQUAL_TO(Actual, Limit)                                                \
MTL::Test::GreaterThanOrEqualTo(Actual, Limit, MTL::String(MTL__FILE__), __LINE__)

#define MTL_APP()                                                                                  \
static void MTL_App_();                                                                            \
static int MTL_App_Initialze_ = MTL::Test::SetAppFunction(MTL_App_);                               \
static void MTL_App_()

namespace MTL
{

static ProgressBar TestProgressBar;

static void ShowProgressBar(double percent, const String& message = L"", int barLength = 50,
                            const ColorRGB& barColor = ColorRGB(0, 255, 0),
                            const ColorRGB& textColor = ColorRGB(0, 255, 255), uint16_t indent = 2)
{
  TestProgressBar.Update(percent, message, barLength, barColor, textColor, indent);
}
static void ShowProgressBar(double percent, const std::string& message, int barLength = 50,
                            const ColorRGB& barColor = ColorRGB(0, 255, 0),
                            const ColorRGB& textColor = ColorRGB(0, 255, 255), uint16_t indent = 2)
{
  ShowProgressBar(percent, ToUTF16(message), barLength, barColor, textColor, indent);
}

static const wchar_t* ErrorColor = COLOR_ERROR;
static const wchar_t* InitializeShutdownColor     = COLOR_FG(60, 120, 240);
static const wchar_t* InitializeShutdownTimeColor = COLOR_FG(60, 200, 240);
static const wchar_t* TestBeginEndColor           = COLOR_FG(120, 100, 255);
static const wchar_t* TestEndTimeColor            = COLOR_FG(180, 180, 255);

class Test
{
public:
  Test(const String& name, bool enabled)
    : Name_(name), TimeElapsed_(0)
  {
    if (enabled)
      List_.PushBack(this);
  }

  virtual ~Test()  {}

  virtual void Run() = 0;

  // For tracking memory leaks with Visual Studio builds.
  static void SetAllocationBreak(long allocationNumber)
  {
#if defined(WIN32) || defined(WIN64)
    _CrtSetBreakAlloc(allocationNumber);
#endif
  }

  static void RunAll()
  {
#if defined(WIN32) || defined(WIN64)
    __try
    {
#endif
      Run_All();
#if defined(WIN32) || defined(WIN64)
    }
    __except(EXCEPTION_EXECUTE_HANDLER)  
    {
      Out() << ErrorColor << L"ERROR: SEH exception caught!" << COLOR_RESET << std::endl;
      TotalNumberOfFailures_++;
    }
#endif
  }

  static OutputStream& Out()
  {
    return Output_;
  }

  static const DynamicVector<String>& Arguments()  { return Arguments_; }

  static I32 FindArgument(const String& str)
  {
    for (U32 i = 0; i < Arguments_.Size(); i++)
    {
      if (str == Arguments_[i])
        return i;
    }

    return -1;
  }

  static I32 FindArgumentIgnoreCase(const String& str)
  {
    String lowerCaseStr = ToLowerCase(str);

    for (U32 i = 0; i < Arguments_.Size(); i++)
    {
      if (lowerCaseStr == ToLowerCase(Arguments_[i]))
        return i;
    }

    return -1;
  }

  static U64 TotalNumberOfFailures()  { return TotalNumberOfFailures_; }

  static const String& TestFilePathName() throw()
  { return FilePath_; }

  static String TestPath() throw()
  {
    String fullPath = FilePath_;
    for (U32 i = 0; i < fullPath.length(); i++)
      if (fullPath[i] == L'\\')
        fullPath[i] = L'/';

    I32 slashIndex = -1;
    for (U32 i = 0; i < fullPath.length(); i++)
      if (fullPath[i] == L'/')
        slashIndex = i;

    if (slashIndex != -1)
      fullPath.erase(slashIndex + 1, fullPath.length() - slashIndex - 1);

    return fullPath;
  }

  static void Verify(bool e, const char* expression, const String& file, U64 line)
  {
    if (!e)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[In file '" << file << L"' - line " << line  << L"]" << std::endl
            << L"  '" << expression << L"' failed!" <<  std::endl;
    }
  }

  template <class T>
  static void Equal(const T& actual, const T& expected, const String& file, U64 line)
  {
    if (actual != expected)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual << L"' but '"
            << expected << L"' is expected!'" << std::endl;
    }
  }

  static void EqualFloat(double actual, double expected, double tolerance, const String& file,
                         U64 line)
  {
    double difference = actual - expected;
    if (MTL::Abs(difference) > tolerance || actual != actual)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line  << L"]" << std::endl
            << L"  Actual value is " << actual << L" but " << expected
            << L" is expected; difference is: " << difference << L", tolerance is: "
            << tolerance << std::endl;
    }
  }

  template <class T>
  static void LessThan(const T& actual, const T& limit, const String& file, U64 line)
  {
    if (actual == limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual
            << L"' which is equal to the limit '"
            << limit << L"' but it is expected to be less than this limit." << std::endl;
    }
    else if (actual > limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual
            << L"' which is greater the limit '"
            << limit << L"'" << std::endl;
    }
  }
  template <class T>
  static void LessThanOrEqualTo(const T& actual, const T& limit, const String& file, U64 line)
  {
    if (actual > limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual << L"' which is greater than the limit '"
            << limit << L"'" << std::endl;
    }
  }

  template <class T>
  static void GreaterThan(const T& actual, const T& limit, const String& file, U64 line)
  {
    if (actual == limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual
            << L"' which is equal to the limit '"
            << limit << L"' but it is expected to be greater than this limit." << std::endl;
    }
    else if (actual < limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual
            << L"' which is less the limit '"
            << limit << L"'" << std::endl;
    }
  }
  template <class T>
  static void GreterThanOrEqualTo(const T& actual, const T& limit, const String& file, U64 line)
  {
    if (actual < limit)
    {
      TotalNumberOfFailures_++;
      ColorScope c(ErrorColor);
      Out() << L"[File '" << file << L"' - line " << line << L"]" << std::endl
            << L"  Actual value is '" << actual << L"' which is less than the limit '"
            << limit << L"'" << std::endl;
    }
  }

  static int SetInitializeFunction(void(*initf)())  { Initialize_ = initf;  return 0; }
  static int SetShutdownFunction(void(*sdf)())      { Shutdown_   = sdf;    return 0; }
  static int SetAppFunction(void(*app)())           { App_        = app;    return 0; }

protected:
  double TimeElapsed_;  // In seconds.

  static String FilePath_;

private:
  String Name_;

  static OutputStream& Output_;
  static DynamicVector<Test*> List_;
  static DynamicVector<String> Arguments_;

  static U64 TotalNumberOfFailures_;
  static double TotalTimeElapsed_;  // In seconds.

  static void (*Initialize_)();
  static void (*Shutdown_)();
  static void (*App_)();

  static void Initialize()
  {
    if (Initialize_)
    {
      try
      {
        Out() << InitializeShutdownColor << L"[Initialize_Test] Begins..." << COLOR_RESET << std::endl;
        Timer timer(true);
        Initialize_();
        timer.Stop();

        Out() << InitializeShutdownColor << L"[Initialize_Test] Ends.";
        Out() << InitializeShutdownTimeColor << L"  Time: " << timer.Seconds() << L" seconds." << COLOR_RESET << std::endl << std::endl;
      }
      catch (const Exception& e)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Initialize: ERROR: " << e.Message() << std::endl;
      }
      catch (const std::exception& e)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Initialize: std::exception: " << ToUTF16(e.what()) << std::endl;
      }
      catch (...)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Initialize: UNEXPECTED EXCEPTION!" << std::endl;
      }
    }
  }

  static void Shutdown()
  {
    if (Shutdown_)
    {
      try
      {
        Out() << InitializeShutdownColor << L"[Shutdown_Test] Begins..." << COLOR_RESET << std::endl;
        Timer timer(true);
        Shutdown_();
        timer.Stop();

        Out() << InitializeShutdownColor << L"[Shutdown_Test] Ends.";
        Out() << InitializeShutdownTimeColor << L"  Time: " << timer.Seconds() << L" seconds." << COLOR_RESET << std::endl << std::endl;
      }
      catch (const Exception& e)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Shutdown: ERROR: " << e.Message() << std::endl;
      }
      catch (const std::exception& e)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Shutdown: std::exception: " << ToUTF16(e.what()) << std::endl;
      }
      catch (...)
      {
        TotalNumberOfFailures_++;
        ColorScope c(ErrorColor);
        Out() << L"Test::Shutdown: UNEXPECTED EXCEPTION!" << std::endl;
      }
    }
  }

  static void Run_All()
  {
    if (List_.Size() > 0)
    {
      if (FindArgumentIgnoreCase(L"-DisableProgressBar") > 0 || FindArgumentIgnoreCase(L"-NoProgressBar") > 0)
        TestProgressBar.Disable();

      if (FindArgumentIgnoreCase(L"-DisableColorRGB24") > 0 || FindArgumentIgnoreCase(L"-NoColorRGB24") > 0)
      {
        // Some outputs don't support all ansi color escape codes.
        // Override with more compatible colors.
        ErrorColor = COLOR_RED;
        InitializeShutdownColor = COLOR_MAGENTA;
        InitializeShutdownTimeColor = COLOR_LBLUE;
        TestBeginEndColor = COLOR_LBLUE;
        TestEndTimeColor = COLOR_LCYAN;
      }

      Out() << COLOR_LCYAN << std::endl << L"Number Of Actual Cores: "
            << MTL::CPU::Instance().NumberOfCores() << std::endl;
      Out() << L"Initial Number Of OpenMP Threads: "
            << MTL::CPU::Instance().NumberOfThreads() << COLOR_RESET << std::endl;

      Out() << std::endl;
      Initialize();

      for (U32 i = 0; i < List_.Size(); i++)
      {
        Out() << TestBeginEndColor << L"[" << List_[i]->Name_ << L"]" << L" Begins..." << COLOR_RESET << std::endl;

        U64 numberOfFailuresBeforeRun = TotalNumberOfFailures_;
        try
        {
          ResetOutputStream();

          Timer timer(true);
          List_[i]->Run();
          timer.Stop();
          List_[i]->TimeElapsed_ = timer.Seconds();

          ResetOutputStream();
        }
        catch (const Exception& e)
        {
          ResetOutputStream();

          TotalNumberOfFailures_++;
          ColorScope c(ErrorColor);
          Out() << L"ERROR: " << e.Message() << std::endl;
        }
        catch (const std::exception& e)
        {
          ResetOutputStream();

          TotalNumberOfFailures_++;
          ColorScope c(ErrorColor);
          Out() << L"std::exception: " << ToUTF16(e.what()) << std::endl;
        }
        catch (...)
        {
          ResetOutputStream();

          TotalNumberOfFailures_++;
          ColorScope c(ErrorColor);
          Out() << L"UNEXPECTED EXCEPTION!" << std::endl;
        }

        Out() << TestBeginEndColor << L"[" << List_[i]->Name_ << L"]";
        if (TotalNumberOfFailures_ > numberOfFailuresBeforeRun)
        {
          Out() << COLOR_LRED << L" Ends with errors.";
        }
        else
        {
          Out() << L" Ends.";
        }
        Out() << TestEndTimeColor << L"  Time: " << float(List_[i]->TimeElapsed_) << L" seconds."
	      << COLOR_RESET << std::endl << std::endl;

        TotalTimeElapsed_ += List_[i]->TimeElapsed_;

        Out().flush();  // Make sure we at least flush the output at the end of every test (std::endl is supposed to flush for some streams).
      }
      Shutdown();

      if (TotalNumberOfFailures() == 0)
        Out() << COLOR_LGREEN;
      else
        Out() << COLOR_LRED;

      Out() << L"Total number of errors: " << TotalNumberOfFailures() << COLOR_RESET << std::endl;
      Out() << COLOR_LCYAN << L"Total time: " << TotalTimeElapsed_ << L" seconds." << COLOR_RESET << std::endl;
    }
    if (App_)
    {
      App_();
    }
  }
};

String Test::FilePath_;
OutputStream& Test::Output_ = ConsoleOut;
DynamicVector<Test*> Test::List_;
DynamicVector<String> Test::Arguments_;
U64 Test::TotalNumberOfFailures_;
double Test::TotalTimeElapsed_;
void(*Test::Initialize_)() = 0;
void(*Test::Shutdown_)() = 0;
void(*Test::App_)() = 0;

}  // namespace MTL


#ifndef MTL_TEST_NO_MAIN

#if defined(WIN32) || defined(WIN64)
  // For debugging memory leaks with Visual C++.
  #define _CRTDBG_MAP_ALLOC
  #include <crtdbg.h>

  #ifndef WIN32_LEAN_AND_MEAN
  #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>

  #ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
  #define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
  #endif

  static bool EnableVTMode()
  {
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hOut == INVALID_HANDLE_VALUE)
      return false;

    DWORD dwMode = 0;
    if (!GetConsoleMode(hOut, &dwMode))
      return false;

    if (!SetConsoleMode(hOut, dwMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING))
      return false;

    return true;
  }
#endif  // WIN32

#if defined(WIN32) || defined(WIN64)
int wmain(int argc, wchar_t* argv[])
#else
int main(int argc, char* argv[])
#endif
{
#if defined(WIN32) || defined(WIN64)
  EnableVTMode();

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

  int status = 0;

  try
  {
    MTL::DynamicVector<MTL::String>& arguments =
      const_cast<MTL::DynamicVector<MTL::String>&>(MTL::Test::Arguments());

    for (int i = 0; i < argc; i++)
#if defined(WIN32) || defined(WIN64)
    {
      arguments.PushBack(argv[i]);
    }
#else
    {
      arguments.PushBack(MTL::ToUTF16(argv[i]));
    }
#endif

    MTL::Test::RunAll();
  }
  catch(const MTL::Exception& e)
  {
    MTL::ColorScope c(MTL::ErrorColor);
    ConsoleOut << L"MAIN() - MTL::Exception: " << e.Message() << std::endl;

#if defined(WIN32) || defined(WIN64)
    status = 100000000;
#else
    status = -1;
#endif
  }
  catch(const std::exception& e)
  {
    MTL::ColorScope c(MTL::ErrorColor);
    ConsoleOut << L"MAIN() - std::exception: " << MTL::ToUTF16(e.what()) << std::endl;

#if defined(WIN32) || defined(WIN64)
    status = 200000000;
#else
    status = -2;
#endif
  }
  catch(...)
  {
    MTL::ColorScope c(MTL::ErrorColor);
    ConsoleOut << L"MAIN() - UNEXPECTED EXCEPTION!" << std::endl;

#if defined(WIN32) || defined(WIN64)
    status = 300000000;
#else
    status = -3;
#endif
  }

#if defined(WIN32) || defined(WIN64)
  return status + (int)MTL::Test::TotalNumberOfFailures();
#else
  return status ? status : MTL::Min(255, (int)MTL::Test::TotalNumberOfFailures());
#endif
}
#endif  // #ifndef MTL_TEST_NO_MAIN

#endif  // MTL_TEST_H
