//
// Math Template Library
//
// Copyright (c) 2026: Seung Jae Lee, https://github.com/seungjlee/MTL
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

// Branch-coverage chip-aways for DynamicVector: growth from PushBack, Insert
// in the middle (drives CopyBackwards), self-assignment short-circuit, PopBack,
// and Reserve-then-Resize.

#include <MTL/Math/DynamicVector.h>
#include <MTL/Tools/Test.h>

using namespace MTL;

TEST(DynamicVectorBranches_PushBackGrowth)
{
  // Many calls to PushBack starting from a default-constructed vector force
  // the Size_ == BufferSize_ branch to trigger reallocation repeatedly.
  DynamicVector<I32> v;
  for (I32 i = 0; i < 1024; i++)
    v.PushBack(i);
  MTL_EQUAL(v.Size(), SizeType(1024));
  MTL_EQUAL(v[0], 0);
  MTL_EQUAL(v[1023], 1023);

  while (v.Size() > 0)
    v.PopBack();
  MTL_EQUAL(v.Size(), SizeType(0));
}

TEST(DynamicVectorBranches_InsertCopiesBackwards)
{
  // Insert at a position before the end triggers CopyBackwards to shift
  // existing elements without overwriting them.
  DynamicVector<I32> v;
  for (I32 i = 0; i < 8; i++)
    v.PushBack(i * 10);

  // Insert two elements in the middle.
  I32 inserted[2] = { 99, 98 };
  v.Insert(v.Begin() + 3, inserted, 2);

  MTL_EQUAL(v.Size(), SizeType(10));
  MTL_EQUAL(v[2], 20);   // unchanged before insertion point
  MTL_EQUAL(v[3], 99);
  MTL_EQUAL(v[4], 98);
  MTL_EQUAL(v[5], 30);   // shifted right
  MTL_EQUAL(v[9], 70);

  // Also drive the single-element overload.
  v.Insert(v.Begin() + 1, I32(-1));
  MTL_EQUAL(v[1], -1);
  MTL_EQUAL(v[2], 10);
}

TEST(DynamicVectorBranches_SelfAssignmentShortCircuit)
{
  // operator= has an `if (this != &rhs)` guard; assigning to itself drives
  // the false side of that branch.
  DynamicVector<F64> v;
  for (I32 i = 0; i < 5; i++)
    v.PushBack(F64(i) + 0.5);

  DynamicVector<F64>* pSelf = &v;
  v = *pSelf;  // self-assignment

  MTL_EQUAL(v.Size(), SizeType(5));
  for (I32 i = 0; i < 5; i++)
    MTL_EQUAL_FLOAT(v[i], F64(i) + 0.5, 0.0);
}

TEST(DynamicVectorBranches_ReserveThenResize)
{
  // Reserve grows the buffer without changing Size_; a subsequent Resize
  // smaller than the reserved capacity drives the in-place Size_=newSize
  // branch in Resize.
  DynamicVector<F32> v;
  v.Reserve(2048);
  v.Resize(64);
  MTL_EQUAL(v.Size(), SizeType(64));
  MTL_VERIFY(v.Capacity() >= SizeType(2048));
}
