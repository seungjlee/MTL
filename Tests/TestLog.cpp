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

#include <MTL/Log.h>
#include <MTL/Timer.h>
#include <MTL/Tools/Test.h>

#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace MTL;

static std::string ReadFileContents(const char* path)
{
  std::ifstream file(path, std::ios::binary);
  std::stringstream stream;
  stream << file.rdbuf();
  return stream.str();
}

static size_t CountOccurrences(const std::string& text, const std::string& substring)
{
  size_t count = 0;
  for (size_t pos = text.find(substring); pos != std::string::npos;
       pos = text.find(substring, pos + substring.size()))
    count++;
  return count;
}

// Basic usage: stream values with the level macros. Colors and millisecond
// timestamps are on by default.
TEST(TestLogBasics)
{
  Log& TheLog = Log::Instance();

  MTL_LOG(Info)    << "Info message with a number: " << 42;
  MTL_LOG(Debug)   << "Debug message; flag is " << true << " and pi is about " << 3.14159265;
  MTL_LOG(Warning) << "Warning message: " << 2.5 << " seconds remaining.";
  MTL_LOG(Error)   << "Error message with code " << -7 << '.';

  // Short forms and the supported operand types.
  MTL_LOG_I << "Strings: " << std::string("std::string") << ", "
            << std::wstring(L"std::wstring") << ", " << L"wide literal.";
  MTL_LOG_I << "Integers: " << I64(-1234567890123LL) << ' ' << U32(4000000000u)
            << ", pointer: " << &TheLog;

  // Makes all queued messages visible before the test ends.
  TheLog.Flush();
}

// Messages below the minimum level are discarded without even evaluating their arguments.
TEST(TestLogLevelFiltering)
{
  Log& TheLog = Log::Instance();
  TheLog.SetLevel(LogLevel::Warning);

  int evaluations = 0;
  MTL_LOG(Debug) << "Never logged: " << ++evaluations;
  MTL_LOG(Info)  << "Never logged: " << ++evaluations;
  MTL_EQUAL(evaluations, 0);

  MTL_LOG(Warning) << "Warnings still get through; evaluations = " << ++evaluations;
  MTL_EQUAL(evaluations, 1);

  MTL_VERIFY(!TheLog.IsEnabled(LogLevel::Info));
  MTL_VERIFY(TheLog.IsEnabled(LogLevel::Error));

  TheLog.SetLevel(LogLevel::Debug);  // Back to the default.
  TheLog.Flush();
}

// Verbose messages only show up when verbose mode is turned on.
// File output gets the same lines with the color escape codes stripped.
TEST(TestLogVerboseAndFileOutput)
{
  const char* kPath = "TestLog_Output.txt";

  Log& TheLog = Log::Instance();
  TheLog.OpenFile(kPath);

  MTL_LOG(Verbose) << "Hidden; verbose mode is off by default.";

  TheLog.SetVerbose(true);
  MTL_LOG(Verbose) << "Visible verbose message.";
  TheLog.SetVerbose(false);
  TheLog.SetLevel(LogLevel::Debug);

  MTL_LOG(Info) << "Info line for the file.";

  TheLog.Flush();
  TheLog.CloseFile();

  std::string contents = ReadFileContents(kPath);
  MTL_VERIFY(contents.find("Hidden") == std::string::npos);
  MTL_VERIFY(contents.find("V Visible verbose message.") != std::string::npos);
  MTL_VERIFY(contents.find("I Info line for the file.") != std::string::npos);
  MTL_VERIFY(contents.find('\033') == std::string::npos);  // No colors in the file.

  // Timestamps look like "2026-07-26 15:04:05.123 I ...".
  size_t position = contents.find(" I Info line");
  MTL_VERIFY(position != std::string::npos);
  MTL_EQUAL(contents[position - 4], '.');  // Millisecond separator.

  remove(kPath);
}

// Concurrent logging from multiple threads; every line arrives intact.
TEST(TestLogMultithreaded)
{
  const char* kPath = "TestLog_Threads.txt";
  const int kThreads = 4;
  const int kMessagesPerThread = 250;

  Log& TheLog = Log::Instance();
  TheLog.EnableConsole(false);
  TheLog.OpenFile(kPath);

  std::vector<std::thread> threads;
  for (int t = 0; t < kThreads; t++)
  {
    threads.emplace_back([t]()
    {
      for (int i = 0; i < kMessagesPerThread; i++)
        MTL_LOG(Info) << "Worker " << t << " message " << i << " end.";
    });
  }
  for (std::thread& thread : threads)
    thread.join();

  TheLog.Flush();
  TheLog.CloseFile();
  TheLog.EnableConsole(true);

  std::string contents = ReadFileContents(kPath);
  MTL_EQUAL(CountOccurrences(contents, "Worker "), size_t(kThreads * kMessagesPerThread));
  MTL_EQUAL(CountOccurrences(contents, " end.\n"), size_t(kThreads * kMessagesPerThread));

  remove(kPath);
}

// Logging calls only compose and queue; the flush thread does the writing.
TEST(TestLogThroughput)
{
  const char* kPath = "TestLog_Throughput.txt";
  const int kMessages = 20000;

  Log& TheLog = Log::Instance();
  TheLog.EnableConsole(false);  // Time the logging cost without terminal output.
  TheLog.OpenFile(kPath);

  Timer timer(true);
  for (int i = 0; i < kMessages; i++)
    MTL_LOG(Info) << "Throughput message " << i << " with a bit of payload text.";
  double queueSeconds = timer.Seconds();

  TheLog.Flush();
  double totalSeconds = timer.Seconds();

  TheLog.CloseFile();
  TheLog.EnableConsole(true);

  std::string contents = ReadFileContents(kPath);
  MTL_EQUAL(CountOccurrences(contents, "Throughput message "), size_t(kMessages));

  MTL_LOG_I << "Queued " << kMessages << " messages in " << GetTimeUTF8(queueSeconds)
            << " (" << queueSeconds / kMessages * 1e9 << " nanoseconds per message); "
            << "written and flushed in " << GetTimeUTF8(totalSeconds) << '.';
  TheLog.Flush();

  remove(kPath);
}
