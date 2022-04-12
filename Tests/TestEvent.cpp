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

#include <MTL/Tools/Test.h>
#include <MTL/Tools/Event.h>

using namespace MTL;

static const int kInitialID = 7;
static const int kThreads = 32;
static const int kWork = 50000;

static std::array<MTL::Event, kThreads> events;

static int Count = 0;

static void ProcessWork(int ID, int nextID)
{
  for (int i = 0; i < kWork; i++)
  {
    events[ID].Wait();
    Count++;
    events[nextID].Signal();
  }
}

TEST(TestEvents)
{
  std::vector<std::thread> threads;

  for (int i = 0; i < kThreads; i++)
    threads.emplace_back(std::thread(ProcessWork, i, (i + 1) % kThreads));

  events[kInitialID].Signal();

  for (std::thread& t : threads)
    if (t.joinable())
      t.join();

  MTL_EQUAL(Count, kThreads * kWork);
}

TEST(TestEventWaitMilliseconds)
{
  Event event;
  bool signaled = event.Wait(50);
  MTL_EQUAL(signaled, false);
}
