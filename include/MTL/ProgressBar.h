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

#include <MTL/WorkerThread.h>

#ifndef MAX_PROGRESS_BUFFER_SIZE
#define MAX_PROGRESS_BUFFER_SIZE 512
#endif

namespace MTL
{

struct ProgressData
{
  ProgressData() {}
  ProgressData(double percent, const String& message, int barLength, const ColorRGB& barColor, const ColorRGB& textColor, uint16_t indent)
    : Percent(percent), Message(message), BarLength(barLength), BarColor(barColor), TextColor(textColor), Indent(indent)
  {
  }

  double Percent;
  String Message;
  int BarLength;
  ColorRGB BarColor;
  ColorRGB TextColor;
  uint16_t Indent;
};

class ProgressBarWorker : public WorkerThread<ProgressData>
{
public:
  ProgressBarWorker(bool finalUpdateIsSynchronous = true)
    : WorkerThread<ProgressData>(L"ProgressBarWorker"), LastIntegerPercentage_(0),
      FinalUpdateIsSynchronous_(finalUpdateIsSynchronous), Enabled_(true)
  {
    MaxWorkQueueSize(1);
  }

  virtual void QueueWork(const ProgressData& dataIn)
  {
    if (Enabled_)
    {
      ProgressData data = dataIn;
      data.Percent = Limit(data.Percent, 0.0, 1.0);

      const int resolution = 1000;

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
  virtual void ProcessWork(const MTL::DynamicVector<ProgressData>& data)
  {
    for (const ProgressData& d : data)
    {
      ShowProgressBar(d.Percent, d.Message, d.BarLength, d.BarColor, d.TextColor, d.Indent);
    }
  }

private:
  bool Enabled_;
  bool FinalUpdateIsSynchronous_;
  int LastIntegerPercentage_;
  Event Finish_;

  void ShowProgressBar(double percent, const String& message, int barLength, const ColorRGB& barColor, const ColorRGB& textColor, uint16_t indent)
  {
    int numberOfBlocksToPrint = int(barLength * percent);
    char buf[MAX_PROGRESS_BUFFER_SIZE];

#ifdef WIN32
    // This code does not work on Ubuntu with default settings.
    {
      ColorScope cs(barColor);

      int index = 0;
      buf[index++] = '\r';
      for (uint16_t i = 0; i < indent; i++)
        buf[index++] = ' ';
      buf[index++] = '[';

      for (int i = 0; i < barLength; i++)
      {
        if (i < numberOfBlocksToPrint)
          buf[index++] = '\xFE';
        else
          buf[index++] = ' ';
      }
      buf[index++] = ']';
      buf[index] = 0;
      std::wcout << buf;
    }
    {
      ColorScope cs(textColor);

      int bytes = snprintf(&buf[0], sizeof(buf), " %.1f%%", float(100.0 * percent));
      buf[bytes] = 0;

      std::wcout << buf;
    }
#else
    String fgColor = ColorRGB::ForegroundColor(textColor);
    String bgColor = ColorRGB::BackgroundColor(barColor);
    String endColor = ColorRGB::BackgroundColor(barColor * 0.2);

    std::wcout << "\r";
    for (uint16_t i = 0; i < indent; i++)
      std::wcout << " ";
    std::wcout << endColor;
    std::wcout << " ";
    std::wcout << bgColor;
    int i = 0;
    for (; i < numberOfBlocksToPrint; i++)
      std::wcout << " ";

    std::wcout << COLOR_BG(0,0,0);
    for (; i < barLength; i++)
      std::wcout << " ";

    std::wcout << endColor;
    std::wcout << " " << COLOR_RESET;
    std::wcout << fgColor;
    int bytes = snprintf(&buf[0], MAX_PROGRESS_BUFFER_SIZE, " %.1f%%", float(100.0 * percent));
    std::wcout << buf << COLOR_RESET;
#endif

    std::wcout << " " << message;
    std::wcout.flush();

    if (FinalUpdateIsSynchronous_ && percent >= 1.0)
      Finish_.Signal();
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

  void Update(double percent, const String& message = L"", int barLength = 50,
              const ColorRGB& barColor = ColorRGB(0, 255, 0), const ColorRGB& textColor = ColorRGB(0, 255, 255), uint16_t indent = 2)
  {
    Worker_.QueueWork(ProgressData(percent, message, barLength, barColor, textColor, indent));
  }
  
  void Disable()  { Worker_.Disable(); }
  void Enable()   { Worker_.Enable();  }
private:
  ProgressBarWorker Worker_;
};

}  // namespace MTL


#endif  // MTL_UTILITIES_H
