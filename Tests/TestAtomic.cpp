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

#include <MTL/Tools/Test.h>
#include <MTL/Math/Point2D.h>
#include <atomic>

using namespace MTL;

TEST(Test_Atomic)
{
  #define CHECK_ATOMIC_LOCK_FREE(T)                                                \
  printf("std::atomic<\033[96m" #T "\033[0m> is always lock free: %s\n\033[0m",            \
         (std::atomic<T>::is_always_lock_free ? "\033[92mTrue" : "\033[95mFalse"))

  CHECK_ATOMIC_LOCK_FREE(int);
  CHECK_ATOMIC_LOCK_FREE(int64_t);
  CHECK_ATOMIC_LOCK_FREE(double);
  CHECK_ATOMIC_LOCK_FREE(Point2D<float>);
  CHECK_ATOMIC_LOCK_FREE(Point2D<int32_t>);
  CHECK_ATOMIC_LOCK_FREE(Point2D<int64_t>);
  CHECK_ATOMIC_LOCK_FREE(Point2D<double>);
}