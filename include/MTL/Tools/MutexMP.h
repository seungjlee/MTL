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


#ifndef MTL_MUTEX_MP_H
#define MTL_MUTEX_MP_H

#include <omp.h>

namespace MTL
{

// An alternative to std::mutex and std::recursive_mutex since g++ 7.5 still uses the inefficient pthread library
// which wastes CPU cycles while waiting for the mutex lock.
class MutexMP
{
public:
  MutexMP()
  {
    omp_init_lock(&Mutex);
  }
  ~MutexMP()
  {
    omp_destroy_lock(&Mutex);
  }
  void lock()
  {
    omp_set_lock(&Mutex);
  }
  void unlock()
  {
    omp_unset_lock(&Mutex);
  }

private:
  omp_lock_t Mutex;
};
class RecursiveMutexMP
{
public:
  RecursiveMutexMP()
  {
    omp_init_nest_lock(&Mutex);
  }
  ~RecursiveMutexMP()
  {
    omp_destroy_nest_lock(&Mutex);
  }
  void lock()
  {
    omp_set_nest_lock(&Mutex);
  }
  void unlock()
  {
    omp_unset_nest_lock(&Mutex);
  }

private:
  omp_nest_lock_t Mutex;
};

class LockMP
{
public:
  LockMP(MutexMP& m)
    : Mutex(m)
  {
    Mutex.lock();
  }
  ~LockMP()
  {
    Mutex.unlock();
  }

private:
  MutexMP& Mutex;
};
class RecursiveLockMP
{
public:
  RecursiveLockMP(RecursiveMutexMP& m)
    : Mutex(m)
  {
    Mutex.lock();
  }
  ~RecursiveLockMP()
  {
    Mutex.unlock();
  }

private:
  RecursiveMutexMP& Mutex;
};

}  // namespace MTL

#endif  // MTL_MUTEX_MP_H
