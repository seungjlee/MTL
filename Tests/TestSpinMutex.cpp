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

#include <MTL/Tools/SpinMutex.h>
#include <MTL/Tools/Test.h>

using namespace MTL;

static const int Work = 4;
static const int MAX_THREADS = 256;
static int NumberOfThreads = 6;
static int NumberOfIterations = 40000;
static int64_t Counter = 0;
static std::array<double, MAX_THREADS> LockTimes;

template <class MutexClass>
static void ThreadProcess(uint32_t ID, MutexClass* mutex) {
  Timer timer(false);
  for (int i = 0; i < NumberOfIterations; i++) {
    timer.Start();
    std::lock_guard<MutexClass> lock(*mutex);
    timer.Stop();

    // Do some work.
    for (int w = 0; w < Work; w++)
      std::this_thread::yield();

    Counter++;
  }

  LockTimes[ID] = timer.Seconds();
}

template <class MutexClass>
static void TestMutex(const std::string& mutexName, int iterations)
{
  std::vector<std::thread> threads(NumberOfThreads);

  MutexClass mutex;
  Counter = 0;

  for (uint32_t tid = 0; tid < threads.size(); tid++)
    threads[tid] = std::thread(ThreadProcess<MutexClass>, tid, &mutex);

  for (auto& t : threads)
    t.join();

  double sum = 0;
  for (uint32_t t = 0; t < threads.size(); t++)
    sum += LockTimes[t];

  int64_t totalCount = threads.size() * iterations;
  printf("  [%-16s] -- Average lock latency: %6.2f usecs.\n", mutexName.c_str(), 1e6 * sum / totalCount);

  MTL_EQUAL(Counter, totalCount);
}

#define TEST_MUTEX(Iterations, MutexClass) \
  TestMutex<MutexClass>(#MutexClass, Iterations);

TEST(Test_SpinMutex)
{
  using SpinMutexNoYield = SpinMutex<0,false>;
  TEST_MUTEX(NumberOfIterations, std::mutex);
  TEST_MUTEX(NumberOfIterations, SpinMutex<0>);
  TEST_MUTEX(NumberOfIterations, SpinMutexNoYield);
  TEST_MUTEX(NumberOfIterations, SpinMutex<1>);
  TEST_MUTEX(NumberOfIterations, SpinMutex<10000>);
}