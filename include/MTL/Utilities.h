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

namespace MTL
{

struct ProgressData
{
  ProgressData() {}
  ProgressData(double percent, bool showFractions, int barLength, const wchar_t* color, int indent)
    : Percent(percent), ShowFractions(showFractions), BarLength(barLength), Color(color), Indent(indent)
  {
  }

  double Percent;
  bool ShowFractions;
  int BarLength;
  const wchar_t* Color;
  int Indent;
};

#ifndef MAX_PROGRESS_BUFFER_SIZE
#define MAX_PROGRESS_BUFFER_SIZE 512
#endif

class ProgressBarWorker : public WorkerThread<ProgressData>
{
public:
  ProgressBarWorker() : WorkerThread<ProgressData>(L"ProgressBarWorker"), LastIntegerPercentage_(0)
  {
    MaxWorkQueueSize(1);
  }

  virtual void QueueWork(const ProgressData& data)
  {
    int resolution = data.ShowFractions ? 1000 : 100;

    if (data.Percent == 0.0 || LastIntegerPercentage_ != int(resolution * data.Percent))
    {
      LastIntegerPercentage_ = int(resolution * data.Percent);
      WorkerThread<ProgressData>::QueueWork(data);
    }
  }

protected:
  virtual void ProcessWork(const MTL::DynamicVector<ProgressData>& data)
  {
    for (const ProgressData& d : data)
    {
      ShowProgressBar(d.Percent, d.ShowFractions, d.BarLength, d.Color, d.Indent);
    }
  }

private:
  int LastIntegerPercentage_;
  void ShowProgressBar(double percent, bool showFractions, int barLength, const wchar_t* color, int indent)
  {
    int numberOfBlocksToPrint = int(barLength * percent);

    char buf[MAX_PROGRESS_BUFFER_SIZE];

    int index = 0;
    buf[index++] = '\r';
    for (int i = 0; i < indent; i++)
      buf[index++] = ' ';
    buf[index++] = '[';

    for (int i = 0; i < barLength; i++)
    {
      if (i < numberOfBlocksToPrint)
        buf[index++] = '\xFE';
      else
        buf[index++] = ' ';
    }

    const char* format = showFractions ? "] %.1f%%" : "] %.0f%%";
    int bytes = snprintf(&buf[index], MAX_PROGRESS_BUFFER_SIZE - index, format, float(100.0 * percent));
    assert(bytes > 0);
    index += bytes;
    buf[index] = 0;

    ColorScope cs(color);
    std::wcout << buf;
    std::wcout.flush();
  }
};

static void ShowProgressBar(double percent, bool showFractions = false, int barLength = 50,
                            const wchar_t* color = COLOR_FG(0, 255, 0), int indent = 2)
{
  static ProgressBarWorker worker;
  worker.QueueWork(ProgressData(percent, showFractions, barLength, color, indent));
}

}  // namespace MTL


#endif  // MTL_UTILITIES_H
