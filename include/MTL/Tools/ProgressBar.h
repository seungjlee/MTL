//
// Math Template Library
//
// Copyright (c) 2016-2019: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_UTILITIES_H
#define MTL_UTILITIES_H

#include <MTL/Tools/WorkerThread.h>

#include <clocale>

#ifndef MAX_PROGRESS_BUFFER_SIZE
#define MAX_PROGRESS_BUFFER_SIZE 512
#endif

namespace MTL
{

// Switch the C locale to UTF-8 at static-initialization time so std::wcout can
// encode wide block glyphs to UTF-8 bytes the terminal can render. Without
// this the default "C" locale drops non-ASCII wide chars and prints them as
// '?'. (Note: std::wcout.imbue(std::locale("en_US.UTF-8")) is not sufficient
// in libstdc++; the wide-stream codecvt follows the C locale.)
namespace detail
{
inline bool SetUtf8CLocale()
{
  static const char* candidates[] =
  {
    "en_US.UTF-8", "en_US.utf8", "C.UTF-8", "C.utf8", ""
  };
  for (const char* name : candidates)
  {
    if (std::setlocale(LC_ALL, name) != nullptr)
      return true;
  }
  return false;
}
static const bool kCLocaleIsUtf8 = SetUtf8CLocale();
}  // namespace detail

struct ProgressData
{
  ProgressData() {}
  ProgressData(double percent, const String& message, bool extraPrecision, int barLength,
               const ColorRGB& barColor, const ColorRGB& textColor, uint16_t indent)
    : Percent(percent), Message(message), ExtraPrecision(extraPrecision), BarLength(barLength),
      BarColor(barColor), TextColor(textColor), Indent(indent)
  {
  }

  double Percent;
  String Message;
  bool ExtraPrecision;
  int BarLength;
  ColorRGB BarColor;
  ColorRGB TextColor;
  uint16_t Indent;
};

class ProgressBarWorker : public WorkerThread<ProgressData>
{
public:
  ProgressBarWorker(bool finalUpdateIsSynchronous = true)
    : WorkerThread<ProgressData>(L"ProgressBarWorker"), Enabled_(true),
      FinalUpdateIsSynchronous_(finalUpdateIsSynchronous), LastIntegerPercentage_(0)
  {
    MaxWorkQueueSize(1);
    Start();
  }

  ~ProgressBarWorker()
  {
    Shutdown();
  }

  virtual void QueueWork(const ProgressData& dataIn)
  {
    if (Enabled_)
    {
      ProgressData data = dataIn;
      data.Percent = Limit(data.Percent, 0.0, 1.0);

      const int resolution = data.ExtraPrecision ? 10000 : 1000;

      if (data.Percent == 0.0 || LastIntegerPercentage_ != int(resolution * data.Percent))
      {
        LastIntegerPercentage_ = int(resolution * data.Percent);
        WorkerThread<ProgressData>::QueueWork(data);

        if (FinalUpdateIsSynchronous_ && data.Percent == 1.0)
          Finish_.Wait();
      }
    }
  }

  void Disable()  { Enabled_ = false; }
  void Enable()   { Enabled_ = true;  }

protected:
  virtual void ProcessWork(const std::vector<ProgressData>& data)
  {
    if (Enabled_)
    {
      const ProgressData& d = data.back();
      ShowProgressBar(d.Percent, d.Message, d.ExtraPrecision, d.BarLength, d.BarColor, d.TextColor, d.Indent);
    }
  }

private:
  bool Enabled_;
  bool FinalUpdateIsSynchronous_;
  int LastIntegerPercentage_;
  Event Finish_;

  void ShowProgressBar(double percent, const String& message, bool extraPrecision, int barLength,
                       const ColorRGB& barColor, const ColorRGB& textColor, uint16_t indent)
  {
    char buf[MAX_PROGRESS_BUFFER_SIZE];

    std::wcout << L"\r";
    for (uint16_t i = 0; i < indent; i++)
      std::wcout << L" ";

    if (barLength > 0)
    {
      // Unicode left-block characters give 8 sub-cell steps per character cell,
      // i.e. 8x the resolution of a plain space-or-block bar.
      static const wchar_t* kPartialBlock[8] =
      {
        L" ",       // 0/8: empty (background only)
        L"\u258F",  // 1/8
        L"\u258E",  // 2/8
        L"\u258D",  // 3/8
        L"\u258C",  // 4/8
        L"\u258B",  // 5/8
        L"\u258A",  // 6/8
        L"\u2589"   // 7/8
      };
      static const wchar_t* kFullBlock = L"\u2588";  // 8/8
      static const ColorRGB kInitialBackground(50, 50, 50);

      const int subCells = barLength * 8;
      const int filledSubCells = int(subCells * percent + 0.5);
      const int fullCells = filledSubCells / 8;
      const int partialIndex = filledSubCells % 8;

      std::wcout << ColorRGB::BackgroundColor(kInitialBackground)
                 << ColorRGB::ForegroundColor(barColor);

      int i = 0;
      for (; i < fullCells; i++)
        std::wcout << kFullBlock;
      if (partialIndex > 0 && i < barLength)
      {
        std::wcout << kPartialBlock[partialIndex];
        i++;
      }
      for (; i < barLength; i++)
        std::wcout << L" ";

      std::wcout << COLOR_RESET << L" ";
    }

    {
      ColorScope cs(textColor);
      const char* percentFormat = extraPrecision ? "%.2f%%" : "%.1f%%";

      int bytes = snprintf(&buf[0], sizeof(buf), percentFormat, float(100.0 * percent));
      buf[bytes] = 0;

      std::wcout << buf;
    }

    std::wcout << " " << message;
    std::wcout.flush();

    if (FinalUpdateIsSynchronous_ && percent >= 1.0)
    {
      Finish_.Signal();
      return;
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
};

class ProgressBar
{
public:
  // If true and percent == 1.0, the update call is synchronous.
  // Note that this should only be changed during set up since it could cause undesired waits without any signals.
  bool FinalUpdateIsSynchronous = true;

  ProgressBar(bool FinalUpdateIsSynchronous = true)
    : Worker_(FinalUpdateIsSynchronous)
  {
  }

  void Update(double percent, const String& message = L"", bool extraPrecision = false, int barLength = 50,
              const ColorRGB& barColor = ColorRGB(0, 255, 0), const ColorRGB& textColor = ColorRGB(0, 255, 255), uint16_t indent = 2)
  {
    Worker_.QueueWork(ProgressData(percent, message, extraPrecision, barLength, barColor, textColor, indent));
  }
  
  void Disable()  { Worker_.Disable(); }
  void Enable()   { Worker_.Enable();  }
private:
  ProgressBarWorker Worker_;
};

}  // namespace MTL


#endif  // MTL_UTILITIES_H
