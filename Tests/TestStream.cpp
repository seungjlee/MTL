//
// Math Template Library
//
// Copyright (c) 2017: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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

#include <MTL/Test.h>
#include <MTL/SSE.h>

using namespace MTL;

TEST(TestShuffle)
{
#if MTL_ENABLE_SSE
  X128<F32> x1234(1.f,2.f,3.f,4.f);
  X128<F32> x5678(5.f,6.f,7.f,8.f);

  X128<F32> x2468 = Shuffle<1|3<<2|1<<4|3<<6>(x1234,x5678);
  MTL_EQUAL(x2468[0],2.f);
  MTL_EQUAL(x2468[1],4.f);
  MTL_EQUAL(x2468[2],6.f);
  MTL_EQUAL(x2468[3],8.f);
#endif

#if MTL_ENABLE_AVX
#endif
}

TEST(TestXOR)
{
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
  XX<I32> xx = XX<I32>(1) ^ XX<I32>(-1);

  for (int k = 0; k < XX<I32>::Increment; k++)
    MTL_EQUAL(xx[k], -2);
#endif
}