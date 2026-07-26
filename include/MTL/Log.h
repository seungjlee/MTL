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

#ifndef MTL_LOG_H
#define MTL_LOG_H

//
// A simple high performance logger. Similar to glog but much simpler.
//
//   MTL_LOG(Info)    << "Value is " << 42;
//   MTL_LOG_W        << "Low disk space: " << percent << '%';
//   MTL_LOG(Verbose) << "Only shown when verbose mode is on.";
//
// - Messages are composed in reusable thread local buffers and queued in a reusable
//   double buffer; steady state logging does not allocate.
// - A background thread does all the writing and flushing (console and optional file)
//   so logging calls return quickly.
// - Colors are on by default: light cyan timestamps with per-level message colors.
// - Timestamps show milliseconds by default (see MTL_DATE_TIME_FORMAT in TimeHelpers.h
//   for the date-time portion).
//

#include <MTL/Colors.h>
#include <MTL/StringHelpers.h>
#include <MTL/Tools/Event.h>
#include <MTL/Tools/Lock.h>
#include <MTL/Tools/SpinMutex.h>

#include <atomic>
#include <chrono>
#include <ctime>
#include <mutex>
#include <stdio.h>
#include <string>
#include <thread>
#include <type_traits>

#ifndef MTL_LOG_DATE_TIME_FORMAT
#define MTL_LOG_DATE_TIME_FORMAT "%F %T"
#endif

namespace MTL
{

enum class LogLevel : int
{
  Verbose = 0,
  Debug   = 1,
  Info    = 2,
  Warning = 3,
  Error   = 4,
  None    = 5   // Use with SetLevel() to turn off all output.
};

class Log
{
public:
  static Log& Instance()
  {
    static Log TheLog;
    return TheLog;
  }

  // Minimum level that gets logged. Defaults to Debug so that verbose messages are off.
  void SetLevel(LogLevel level)     { MinLevel_ = int(level);                     }
  LogLevel Level() const            { return LogLevel(MinLevel_.load());          }

  void SetVerbose(bool on = true)
  {
    Verbose_ = on;
    if (on)
      MinLevel_ = int(LogLevel::Verbose);
  }
  bool Verbose() const              { return Verbose_;                            }

  bool IsEnabled(LogLevel level) const
  {
    return int(level) >= MinLevel_.load(std::memory_order_relaxed) &&
           (level != LogLevel::Verbose || Verbose_.load(std::memory_order_relaxed));
  }

  void EnableColors(bool on = true) { Colors_ = on;                               }
  bool ColorsEnabled() const        { return Colors_;                             }

  void EnableConsole(bool on = true){ Console_ = on;                              }

  // Optional file output. Color escape codes are stripped from the file output.
  void OpenFile(const std::string& path, bool append = false)
  {
    GenericLock<std::mutex> lock(WriteMutex_);
    if (File_)
      fclose(File_);
#ifdef WIN32
    fopen_s(&File_, path.c_str(), append ? "ab" : "wb");
#else
    File_ = fopen(path.c_str(), append ? "ab" : "wb");
#endif
  }
  void OpenFile(const std::wstring& path, bool append = false)
  {
    OpenFile(ToUTF8(path), append);
  }
  void CloseFile()
  {
    GenericLock<std::mutex> lock(WriteMutex_);
    if (File_)
    {
      fclose(File_);
      File_ = nullptr;
    }
  }

  // Queues a composed message; the flush thread picks it up. Fast: one spin lock and a copy.
  void Submit(const char* data, size_t size)
  {
    {
      GenericLock<SpinMutex<0>> lock(QueueMutex_);
      Front_.append(data, size);
    }
    DataReady_.Signal();
  }

  // Blocks until everything queued so far has been written out.
  void Flush()
  {
    WriteOut();
  }

private:
  Log()
  {
    Front_.reserve(kInitialBufferSize);
    Back_.reserve(kInitialBufferSize);
    Thread_ = std::thread(&Log::FlushThread, this);
  }
  ~Log()
  {
    Running_ = false;
    DataReady_.Signal();
    if (Thread_.joinable())
      Thread_.join();
    WriteOut();
    if (File_)
      fclose(File_);
  }
  Log(const Log&) = delete;
  Log& operator=(const Log&) = delete;

  void FlushThread()
  {
    while (Running_)
    {
      DataReady_.Wait(kFlushPeriodMilliseconds);
      WriteOut();
    }
  }

  void WriteOut()
  {
    GenericLock<std::mutex> writeLock(WriteMutex_);
    {
      GenericLock<SpinMutex<0>> lock(QueueMutex_);
      Front_.swap(Back_);
    }
    if (Back_.empty())
      return;

    if (Console_)
    {
      WidenTo(Back_, WideBuffer_);
      fputws(WideBuffer_.c_str(), stdout);
      fflush(stdout);
    }
    if (File_)
    {
      StripColors(Back_, FileBuffer_);
      fwrite(FileBuffer_.data(), 1, FileBuffer_.size(), File_);
      fflush(File_);
    }
    Back_.clear();  // Keeps capacity; no allocation in steady state.
  }

  // Widens into a reused buffer. Fast path for ASCII which covers typical log output;
  // falls back to a full UTF-8 conversion otherwise.
  static void WidenTo(const std::string& src, std::wstring& dst)
  {
    bool ascii = true;
    for (char c : src)
    {
      if (uint8_t(c) >= 0x80)
      {
        ascii = false;
        break;
      }
    }
    if (ascii)
    {
      dst.resize(src.size());
      for (size_t i = 0; i < src.size(); i++)
        dst[i] = wchar_t(uint8_t(src[i]));
    }
    else
    {
      dst = ToUTF16(src);
    }
  }

  // Removes ANSI escape sequences ("\033[...m") for clean file output.
  static void StripColors(const std::string& src, std::string& dst)
  {
    dst.clear();
    for (size_t i = 0; i < src.size(); i++)
    {
      if (src[i] == '\033' && i + 1 < src.size() && src[i + 1] == '[')
      {
        i += 2;
        while (i < src.size() && src[i] != 'm')
          i++;
      }
      else
      {
        dst += src[i];
      }
    }
  }

  static const size_t kInitialBufferSize = 1 << 16;
  static const int64_t kFlushPeriodMilliseconds = 100;

  std::atomic<int>  MinLevel_{int(LogLevel::Debug)};
  std::atomic<bool> Verbose_{false};
  std::atomic<bool> Colors_{true};
  std::atomic<bool> Console_{true};
  std::atomic<bool> Running_{true};

  SpinMutex<0> QueueMutex_;
  std::string Front_;        // Producers append here.
  std::string Back_;         // Flush thread writes from here.
  std::mutex WriteMutex_;
  std::wstring WideBuffer_;  // Reused console conversion buffer.
  std::string FileBuffer_;   // Reused color-stripped file buffer.
  Event DataReady_;
  FILE* File_ = nullptr;
  std::thread Thread_;
};

//
// Composes a single log line in a reusable thread local buffer and submits it on destruction.
//
class LogRecord
{
public:
  LogRecord(LogLevel level)
    : Buffer_(ComposeBuffer()), Colors_(Log::Instance().ColorsEnabled())
  {
    Buffer_.clear();
    if (Colors_)
      Buffer_ += COLOR_UTF8_LCYAN;
    AppendTimestamp(Buffer_);
    Buffer_ += ' ';
    if (Colors_)
      Buffer_ += LevelColor(level);
    Buffer_ += LevelChar(level);
    Buffer_ += ' ';
  }
  ~LogRecord()
  {
    if (Colors_)
      Buffer_ += COLOR_UTF8_RESET;
    Buffer_ += '\n';
    Log::Instance().Submit(Buffer_.data(), Buffer_.size());
  }

  LogRecord& operator<<(const char* text)          { Buffer_ += text;               return *this; }
  LogRecord& operator<<(const std::string& text)   { Buffer_ += text;               return *this; }
  LogRecord& operator<<(const wchar_t* text)       { Buffer_ += ToUTF8(text);       return *this; }
  LogRecord& operator<<(const std::wstring& text)  { Buffer_ += ToUTF8(text);       return *this; }
  LogRecord& operator<<(char c)                    { Buffer_ += c;                  return *this; }
  LogRecord& operator<<(bool b)                    { Buffer_ += b ? "true" : "false"; return *this; }

  template<class T, typename std::enable_if<std::is_integral<T>::value &&
                                            std::is_signed<T>::value, int>::type = 0>
  LogRecord& operator<<(T value)
  {
    char text[24];
    Buffer_.append(text, snprintf(text, sizeof(text), "%lld", (long long)value));
    return *this;
  }
  template<class T, typename std::enable_if<std::is_integral<T>::value &&
                                            std::is_unsigned<T>::value, int>::type = 0>
  LogRecord& operator<<(T value)
  {
    char text[24];
    Buffer_.append(text, snprintf(text, sizeof(text), "%llu", (unsigned long long)value));
    return *this;
  }
  template<class T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
  LogRecord& operator<<(T value)
  {
    char text[32];
    Buffer_.append(text, snprintf(text, sizeof(text), "%.9g", double(value)));
    return *this;
  }
  LogRecord& operator<<(const void* pointer)
  {
    char text[24];
    Buffer_.append(text, snprintf(text, sizeof(text), "%p", pointer));
    return *this;
  }

private:
  static std::string& ComposeBuffer()
  {
    static thread_local std::string Buffer;
    return Buffer;
  }

  // Caches the formatted date-time per second so most calls only format milliseconds.
  static void AppendTimestamp(std::string& buffer)
  {
    int64_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now().time_since_epoch()).count();
    time_t seconds = time_t(milliseconds / 1000);

    static thread_local time_t CachedSeconds = -1;
    static thread_local char CachedText[64];
    if (seconds != CachedSeconds)
    {
      tm timeInfo;
#ifdef WIN32
      localtime_s(&timeInfo, &seconds);
#else
      localtime_r(&seconds, &timeInfo);
#endif
      strftime(CachedText, sizeof(CachedText), MTL_LOG_DATE_TIME_FORMAT, &timeInfo);
      CachedSeconds = seconds;
    }
    buffer += CachedText;

    char msText[8];
    buffer.append(msText, snprintf(msText, sizeof(msText), ".%03d", int(milliseconds % 1000)));
  }

  static char LevelChar(LogLevel level)
  {
    static const char kChars[] = { 'V', 'D', 'I', 'W', 'E' };
    return kChars[int(level)];
  }
  static const char* LevelColor(LogLevel level)
  {
    static const char* const kColors[] =
    {
      COLOR_UTF8_LBLACK,   // Verbose
      COLOR_UTF8_LGREEN,   // Debug
      COLOR_UTF8_LWHITE,   // Info
      COLOR_UTF8_LYELLOW,  // Warning
      COLOR_UTF8_LRED      // Error
    };
    return kColors[int(level)];
  }

  std::string& Buffer_;
  bool Colors_;
};

// Helper for the MTL_LOG macro; makes the expression void so it composes with if/else.
struct LogVoidify
{
  void operator&(const LogRecord&) {}
};

}  // namespace MTL

// Usage: MTL_LOG(Info) << "message"; with levels Verbose, Debug, Info, Warning, Error.
// Arguments are not evaluated when the level is disabled.
#define MTL_LOG(Level)                                                                             \
  !MTL::Log::Instance().IsEnabled(MTL::LogLevel::Level) ? (void)0 :                                \
  MTL::LogVoidify() & MTL::LogRecord(MTL::LogLevel::Level)

#define MTL_LOG_V  MTL_LOG(Verbose)
#define MTL_LOG_D  MTL_LOG(Debug)
#define MTL_LOG_I  MTL_LOG(Info)
#define MTL_LOG_W  MTL_LOG(Warning)
#define MTL_LOG_E  MTL_LOG(Error)

#endif  // MTL_LOG_H
