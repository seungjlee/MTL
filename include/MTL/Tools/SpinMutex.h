//
// Math Template Library
//
// Copyright (c) 2020: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_SPIN_MUTEX_H
#define MTL_SPIN_MUTEX_H

#include <atomic>
#include <thread>

namespace MTL
{

template <uint64_t NANO_SECONDS, bool YIELD = true>
class SpinMutex
{
public:
  SpinMutex() : Locked(false)
  {
  }

  void lock()
  {
    while (Locked.exchange(true))
    {
      if (YIELD)
      {
        if (NANO_SECONDS == 0)
        {
          std::this_thread::yield();
        }
        else
        {
          // Note that the minimum sleep time will depend on std implementation
          // and OS.
          std::this_thread::sleep_for(std::chrono::nanoseconds(NANO_SECONDS));
        }
      }
    }
  }
  void unlock()
  {
    Locked.store(false);
  }

private:
  std::atomic<bool> Locked;
};

}  // namespace MTL

#endif  // MTL_SPIN_MUTEX_H
