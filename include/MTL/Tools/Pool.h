//
// Math Template Library
//
// Copyright (c) 2021: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#ifndef MTL_POOL_H
#define MTL_POOL_H

#include <cassert>
#include <mutex>
#include <vector>

namespace MTL
{

// Simple thread-safe pool.
template <class T>
class Pool
{
 public:
  Pool() : StackTop_(-1) {}
  ~Pool() {
    assert(StackTop_ == (int64_t)Elements_.size() - 1);
    for (auto ptr : Elements_)
      delete ptr;
  }

  T* CheckOut()
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    T* ptr;
    assert(StackTop_ >= -1);

    if (StackTop_ == -1)
    {
      ptr = new T();
      Elements_.resize(Elements_.size() + 1);
    }
    else
    {
      ptr = Elements_[StackTop_--];
    }
    return ptr;
  }
  void CheckIn(T* ptr)
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    StackTop_++;
    assert(StackTop_ < (int64_t)Elements_.size());
    Elements_[StackTop_] = ptr;
  }

  uint64_t Size() const
  {
    return Elements_.size();
  }

 protected:
  std::mutex Mutex_;
  std::vector<T*> Elements_;
  int64_t StackTop_;
};


}  // namespace MTL

#endif  // MTL_POOL_H
