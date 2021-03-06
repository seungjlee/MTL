//
// Math Template Library
//
// Copyright (c) 2019: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_EVENT_H
#define MTL_EVENT_H

#include <condition_variable>

namespace MTL
{

class Event
{
public:
  Event() : Ready_(false)
  {
  }

  void Wait()
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    while (!Ready_)
      ConditionVariable_.wait(Mutex_);
    Ready_ = false;
  }

  // Returns false if it times out (TimeOut in milliseconds).
  bool Wait(int64_t TimeOut)
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    while (!Ready_)
    {
      std::cv_status status = ConditionVariable_.wait_for(Mutex_, std::chrono::milliseconds(TimeOut));
      if (status == std::cv_status::timeout)
        return false;
    }

    Ready_ = false;
    return true;
  }

  void Signal()
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    Ready_ = true;
    ConditionVariable_.notify_all();
  }

private:
  std::mutex Mutex_;
  std::condition_variable_any ConditionVariable_;
  bool Ready_;
};

}  // namespace MTL

#endif  // MTL_EVENT_H
