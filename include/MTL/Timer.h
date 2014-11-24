//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_TIMER_H
#define MTL_TIMER_H

#include "Definitions.h"
#include <sys/timeb.h>

#if defined(WIN32) || defined(WIN64)
#define WINDOWS_LEAN_AND_MEAN 1
#include <windows.h>

static double InitClockFactorToMilliseconds()
{
  LARGE_INTEGER frequency;
  QueryPerformanceFrequency(&frequency);
  return 1000. / frequency.QuadPart;
}

static const double kClockFactorToMilliseconds = InitClockFactorToMilliseconds();

#endif

namespace MTL
{

class Timer
{
public:
  Timer(bool start = false)
    : Running_(false), ElapsedTimeMilliseconds_(0)
  {
    if (start)
      Start();
  }

  void Start()
  {
    if (!Running_)
    {
      StartTimeMilliseconds_ = getCurrentTimeInMilliseconds();
      Running_ = true;
    }
  }

  void Stop()
  {
    if (Running_)
    {
      ElapsedTimeMilliseconds_ += getElapsedTimeInMilliseconds();

      Running_ = false;
    }
  }

  void Reset()
  {
    ElapsedTimeMilliseconds_ = 0;
    if (Running_)
      StartTimeMilliseconds_ = getCurrentTimeInMilliseconds();
  }

  void ResetAndStart()  { Reset(); Start(); }

  double Seconds() const  { return Milliseconds() * 1e-3; }
  double Milliseconds() const
  {
    if (Running_)
      return
        ElapsedTimeMilliseconds_ + getElapsedTimeInMilliseconds();
    else
      return ElapsedTimeMilliseconds_;
  }

private:
  bool Running_;
  double ElapsedTimeMilliseconds_;
  double StartTimeMilliseconds_;

  double getCurrentTimeInMilliseconds() const
  {
#if defined(WIN32) || defined(WIN64)
    LARGE_INTEGER count;
    QueryPerformanceCounter(&count);
    return double(count.QuadPart) * kClockFactorToMilliseconds;
#else
    timeb t;
    ftime(&t);
    return double(t.time * 1000 + t.millitm);
#endif
  }

  double getElapsedTimeInMilliseconds() const
  {
    return getCurrentTimeInMilliseconds() - StartTimeMilliseconds_;
  }
};

}  // namespace MTL

#endif  // MTL_TIMER_H
