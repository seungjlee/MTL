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
#include <MTL/Tools/WorkerThread.h>

using namespace MTL;

static const int MaxCount = 100000;

class TestWorker : public WorkerThread<int>
{
public:
  int Count;
  int ID;
  TestWorker* pNext;

  TestWorker(int id, TestWorker* next = nullptr)
    : WorkerThread<int>(L"TestWorker" + std::to_wstring(id)),
      ID(id), Count(0), pNext(next)
  {
    MaxWorkQueueSize(MaxCount);
  }

  virtual void CleanupThread()
  {
    std::lock_guard<std::recursive_mutex> lock(QueueMutex_);
    ProcessWork(QueueData_);
    QueueData_.Clear();
  }

  virtual void ProcessWork(const MTL::DynamicVector<int>& data)
  {
    for (int val : data)
    {
      if (val != Count)
      {
        ColorScope cs(COLOR_RED);
        wprintf(L"%d: Invalid sequence %d, expected %d!\n", ID, val, Count);
        MTL_VERIFY(false);
      }

      if (pNext)
        pNext->QueueWork(val);

      Count++;
    }
  }
};

TEST(TestWorkerPipelines)
{
  const int kRepeats = 5;
  int Pipelines = 4;

  std::vector<std::shared_ptr<TestWorker>> master(Pipelines);
  std::vector<std::shared_ptr<TestWorker>> slave1(Pipelines);
  std::vector<std::shared_ptr<TestWorker>> slave2(Pipelines);

  for (int iteration = 0; iteration < kRepeats; iteration++)
  {
    for (int i = 0; i < Pipelines; i++)
    {
      slave2[i].reset(new TestWorker(2000 + i));
      slave1[i].reset(new TestWorker(2000 + i, slave2[i].get()));
      master[i].reset(new TestWorker(1000 + i, slave1[i].get()));
    }

    for (int k = 0; k < MaxCount; k++)
    {
      for (int i = 0; i < Pipelines; i++)
        master[i]->QueueWork(k);
    }

    for (int i = 0; i < Pipelines; i++)
    {
      master[i]->Shutdown();
      slave1[i]->Shutdown();
      slave2[i]->Shutdown();

      MTL_EQUAL(master[i]->Count, MaxCount);
      MTL_EQUAL(slave1[i]->Count, MaxCount);
      MTL_EQUAL(slave2[i]->Count, MaxCount);
    }
  }
}
