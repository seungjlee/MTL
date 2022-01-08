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


#ifndef MTL_WORKER_THREAD_H
#define MTL_WORKER_THREAD_H

#include <MTL/Exception.h>
#include <MTL/Colors.h>
#include <MTL/Tools/Event.h>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#if defined(WIN32) && !defined(MTL_NO_THREAD_NAME)
static bool SetThreadName(unsigned long dwThreadID, char* threadName);
#endif

namespace MTL
{

template<class DataType, class VectorClass = std::vector<DataType>>
class WorkerThread
{
public:
  WorkerThread(const std::wstring& name, bool start = true)
    : Running_(true), Name_(name), MaxWorkQueueSize_(10000)
  {
    if (start)
      Start();
  }
  WorkerThread(const std::string& name, bool start = true)
    : WorkerThread(ToUTF16(name), start)
  {
  }
  ~WorkerThread()
  {
    if (Thread_.joinable())
    {
      ResetOutputStream();
      std::wcerr << COLOR_YELLOW;
      std::wcerr << "WorkerThread::~WorkerThread: Shutdown() needs to be called before object destruction to ensure proper thread clean up.";
      std::wcerr << COLOR_RESET << std::endl;
    }
  }

  uint32_t MaxWorkQueueSize() const     { return MaxWorkQueueSize_;  }
  void MaxWorkQueueSize(uint32_t size)  { MaxWorkQueueSize_ = size;  }

  void Start()
  {
    Running_ = true;
    Thread_ = std::thread(&WorkerThread::ProcessThread, this);

#if defined(WIN32) && !defined(MTL_NO_THREAD_NAME)
    std::string name = ToUTF8(Name_);
    SetThreadName(GetThreadId(Thread_.native_handle()), &name[0]);
#endif
  }

  void Pause()
  {
    Running_ = false;
  }
  void Continue()
  {
    Running_ = true;
  }

  void ClearQueue()
  {
    std::lock_guard<std::recursive_mutex> lock(QueueMutex_);
    QueueData_.clear();
  }

  void Shutdown()
  {
    Running_ = false;
    if (Thread_.joinable())
    {
      ProcessData_.Signal();
      Thread_.join();
    }
  }

  virtual void QueueWork(const DataType& data)
  {
    std::lock_guard<std::recursive_mutex> lock(QueueMutex_);
    QueueData_.push_back(data);
    ProcessData_.Signal();
  }

  bool IsRunning() const
  {
    return Running_;
  }

protected:
  // Overridables to add code for in-thread initialization and cleanup.
  virtual void InitializeThread() {}
  virtual void CleanupThread() {}

  // Executes before waiting for work.
  virtual void ProcessBeforeWait() {}

  virtual void ProcessWork(const VectorClass& data) = 0;

  virtual void HandleException(const MTL::Exception& ex, const MTL::String& functionName)
  {
    ColorScope cs(COLOR_RED);
    wprintf(L"[%s] %s: Fatal error! ", Name_.c_str(), functionName.c_str());
    wprintf(L"%s\n", ex.Message().c_str());

    throw ex;
  }
  virtual void HandleException(const std::exception& ex, const MTL::String& functionName)
  {
    ColorScope cs(COLOR_RED);
    MTL::String msg = ToUTF16(ex.what());
    wprintf(L"[%s] %s: Fatal error! %s\n", Name_.c_str(), functionName.c_str(), msg.c_str());

    throw ex;
  }
  virtual void HandleException(const MTL::String& functionName)
  {
    ColorScope cs(COLOR_RED);
    wprintf(L"[%s] %s: Fatal error!\n", Name_.c_str(), functionName.c_str());

    MTL_THROW("Unexpected exception!");
  }

protected:
  bool Running_;
  MTL::String Name_;
  MTL::Event ProcessData_;
  std::thread Thread_;
  std::recursive_mutex QueueMutex_;
  VectorClass QueueData_;
  VectorClass ThreadWorkData_;
  uint32_t MaxWorkQueueSize_;

  void ProcessThread()
  {
    try
    {
      InitializeThread();
      while (true)
      {
        ProcessBeforeWait();
        ProcessData_.Wait();
        if (!Running_)
          break;

        if (QueueData_.size() > 0)
        {
          {
            std::lock_guard<std::recursive_mutex> lock(QueueMutex_);
            if (QueueData_.size() > MaxWorkQueueSize_)
            {
              uint32_t startOffset = (uint32_t)QueueData_.size() - MaxWorkQueueSize_;
              ThreadWorkData_ = std::vector<DataType>(QueueData_.begin() + startOffset, QueueData_.end());
            }
            else
            {
              ThreadWorkData_ = QueueData_;
            }
            QueueData_.clear();
          }
          ProcessWork(ThreadWorkData_);
        }
      }
      CleanupThread();
    }
    catch (const MTL::Exception& ex)
    {
      HandleException(ex, ToUTF16(__FUNCTION__));
    }
    catch (const std::exception& ex)
    {
      HandleException(ex, ToUTF16(__FUNCTION__));
    }
    catch (...)
    {
      HandleException(ToUTF16(__FUNCTION__));
    }
  }
};

}  // namespace MTL

#if defined(WIN32) && !defined(MTL_NO_THREAD_NAME)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <processthreadsapi.h>

static const DWORD MS_VC_EXCEPTION = 0x406D1388;
#pragma pack(push,8)
typedef struct tagTHREADNAME_INFO
{
  DWORD dwType;      // Must be 0x1000.
  LPCSTR szName;     // Pointer to name (in user addr space).
  DWORD dwThreadID;  // Thread ID (-1=caller thread).
  DWORD dwFlags;     // Reserved for future use, must be zero.
} THREADNAME_INFO;
#pragma pack(pop)

// Adapted from Microsoft's sample code.
static bool SetThreadName(unsigned long dwThreadID, char* threadName)
{
  THREADNAME_INFO info;
  info.dwType = 0x1000;
  info.szName = threadName;
  info.dwThreadID = dwThreadID;
  info.dwFlags = 0;

  __try
  {
    RaiseException(MS_VC_EXCEPTION, 0, sizeof(info) / sizeof(ULONG_PTR), (ULONG_PTR*)&info);
  }
  __except (EXCEPTION_EXECUTE_HANDLER)
  {
    return false;
  }

  return true;
}
#endif

#endif  // MTL_WORKER_THREAD_H
