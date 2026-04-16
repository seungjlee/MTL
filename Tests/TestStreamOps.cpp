//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/Stream/StreamArray.h>

//
// I8 operator+ and operator-
//
TEST(TestI8_AddSub)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I8> a((MTL::I8)10);
    MTL::X128<MTL::I8> b((MTL::I8)3);
    MTL::X128<MTL::I8> sum = a + b;
    MTL::X128<MTL::I8> diff = a - b;
    for (int i = 0; i < MTL::X128<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::I8)13);
      MTL_EQUAL(diff[i], (MTL::I8)7);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I8> a((MTL::I8)10);
    MTL::X256<MTL::I8> b((MTL::I8)3);
    MTL::X256<MTL::I8> sum = a + b;
    MTL::X256<MTL::I8> diff = a - b;
    for (int i = 0; i < MTL::X256<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::I8)13);
      MTL_EQUAL(diff[i], (MTL::I8)7);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I8> a((MTL::I8)10);
    MTL::X512<MTL::I8> b((MTL::I8)3);
    MTL::X512<MTL::I8> sum = a + b;
    MTL::X512<MTL::I8> diff = a - b;
    for (int i = 0; i < MTL::X512<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::I8)13);
      MTL_EQUAL(diff[i], (MTL::I8)7);
    }
  }
#endif
}

//
// U8 operator+ and operator-
//
TEST(TestU8_AddSub)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::U8> a((MTL::U8)200);
    MTL::X128<MTL::U8> b((MTL::U8)50);
    MTL::X128<MTL::U8> sum = a + b;
    MTL::X128<MTL::U8> diff = a - b;
    for (int i = 0; i < MTL::X128<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::U8)250);
      MTL_EQUAL(diff[i], (MTL::U8)150);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::U8> a((MTL::U8)200);
    MTL::X256<MTL::U8> b((MTL::U8)50);
    MTL::X256<MTL::U8> sum = a + b;
    MTL::X256<MTL::U8> diff = a - b;
    for (int i = 0; i < MTL::X256<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::U8)250);
      MTL_EQUAL(diff[i], (MTL::U8)150);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::U8> a((MTL::U8)200);
    MTL::X512<MTL::U8> b((MTL::U8)50);
    MTL::X512<MTL::U8> sum = a + b;
    MTL::X512<MTL::U8> diff = a - b;
    for (int i = 0; i < MTL::X512<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::U8)250);
      MTL_EQUAL(diff[i], (MTL::U8)150);
    }
  }
#endif
}

//
// I8 / U8 bitwise operators (&, |, ^)
//
TEST(TestI8_Bitwise)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I8> a((MTL::I8)0x0F);
    MTL::X128<MTL::I8> b((MTL::I8)0x33);
    MTL::X128<MTL::I8> andResult = a & b;
    MTL::X128<MTL::I8> orResult  = a | b;
    MTL::X128<MTL::I8> xorResult = a ^ b;
    for (int i = 0; i < MTL::X128<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(andResult[i], (MTL::I8)(0x0F & 0x33));
      MTL_EQUAL(orResult[i],  (MTL::I8)(0x0F | 0x33));
      MTL_EQUAL(xorResult[i], (MTL::I8)(0x0F ^ 0x33));
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I8> a((MTL::I8)0x0F);
    MTL::X256<MTL::I8> b((MTL::I8)0x33);
    MTL::X256<MTL::I8> andResult = a & b;
    MTL::X256<MTL::I8> orResult  = a | b;
    MTL::X256<MTL::I8> xorResult = a ^ b;
    for (int i = 0; i < MTL::X256<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(andResult[i], (MTL::I8)(0x0F & 0x33));
      MTL_EQUAL(orResult[i],  (MTL::I8)(0x0F | 0x33));
      MTL_EQUAL(xorResult[i], (MTL::I8)(0x0F ^ 0x33));
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I8> a((MTL::I8)0x0F);
    MTL::X512<MTL::I8> b((MTL::I8)0x33);
    MTL::X512<MTL::I8> andResult = a & b;
    MTL::X512<MTL::I8> orResult  = a | b;
    MTL::X512<MTL::I8> xorResult = a ^ b;
    for (int i = 0; i < MTL::X512<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(andResult[i], (MTL::I8)(0x0F & 0x33));
      MTL_EQUAL(orResult[i],  (MTL::I8)(0x0F | 0x33));
      MTL_EQUAL(xorResult[i], (MTL::I8)(0x0F ^ 0x33));
    }
  }
#endif
}

//
// I32 arithmetic and comparison operators
//
TEST(TestI32_Ops)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I32> a((MTL::I32)100);
    MTL::X128<MTL::I32> b((MTL::I32)42);
    MTL::X128<MTL::I32> sum  = a + b;
    MTL::X128<MTL::I32> diff = a - b;
    MTL::X128<MTL::I32> neg  = -a;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(sum[i], 142);
      MTL_EQUAL(diff[i], 58);
      MTL_EQUAL(neg[i], -100);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I32> a((MTL::I32)100);
    MTL::X256<MTL::I32> b((MTL::I32)42);
    MTL::X256<MTL::I32> sum  = a + b;
    MTL::X256<MTL::I32> diff = a - b;
    MTL::X256<MTL::I32> neg  = -a;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(sum[i], 142);
      MTL_EQUAL(diff[i], 58);
      MTL_EQUAL(neg[i], -100);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I32> a((MTL::I32)100);
    MTL::X512<MTL::I32> b((MTL::I32)42);
    MTL::X512<MTL::I32> sum  = a + b;
    MTL::X512<MTL::I32> diff = a - b;
    MTL::X512<MTL::I32> neg  = -a;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(sum[i], 142);
      MTL_EQUAL(diff[i], 58);
      MTL_EQUAL(neg[i], -100);
    }
  }
#endif
}

//
// I32 shift operators
//
TEST(TestI32_Shift)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I32> a((MTL::I32)0x100);
    MTL::X128<MTL::I32> left  = a << 4;
    MTL::X128<MTL::I32> right = a >> 4;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(left[i],  (MTL::I32)0x1000);
      MTL_EQUAL(right[i], (MTL::I32)0x10);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I32> a((MTL::I32)0x100);
    MTL::X256<MTL::I32> left  = a << 4;
    MTL::X256<MTL::I32> right = a >> 4;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(left[i],  (MTL::I32)0x1000);
      MTL_EQUAL(right[i], (MTL::I32)0x10);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I32> a((MTL::I32)0x100);
    MTL::X512<MTL::I32> left  = a << 4;
    MTL::X512<MTL::I32> right = a >> 4;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(left[i],  (MTL::I32)0x1000);
      MTL_EQUAL(right[i], (MTL::I32)0x10);
    }
  }
#endif
}

//
// I32 compound assignment operators (+=, -=, &=, |=, ^=, <<=, >>=)
//
TEST(TestI32_CompoundAssign)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I32> a((MTL::I32)10);
    MTL::X128<MTL::I32> b((MTL::I32)3);
    a += b;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 13);
    a -= b;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
    a <<= 2;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 40);
    a >>= 2;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I32> a((MTL::I32)10);
    MTL::X256<MTL::I32> b((MTL::I32)3);
    a += b;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 13);
    a -= b;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
    a <<= 2;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 40);
    a >>= 2;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I32> a((MTL::I32)10);
    MTL::X512<MTL::I32> b((MTL::I32)3);
    a += b;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 13);
    a -= b;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
    a <<= 2;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 40);
    a >>= 2;
    for (int i = 0; i < MTL::X512<MTL::I32>::Increment; i++)
      MTL_EQUAL(a[i], 10);
  }
#endif
}

//
// I16 multiply and MultiplyHi
//
TEST(TestI16_Multiply)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I16> a((MTL::I16)100);
    MTL::X128<MTL::I16> b((MTL::I16)200);
    MTL::X128<MTL::I16> prod = a * b;
    // 100 * 200 = 20000. Low 16 bits of 20000 = 20000 (fits in I16).
    for (int i = 0; i < MTL::X128<MTL::I16>::Increment; i++)
      MTL_EQUAL(prod[i], (MTL::I16)20000);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I16> a((MTL::I16)100);
    MTL::X256<MTL::I16> b((MTL::I16)200);
    MTL::X256<MTL::I16> prod = a * b;
    for (int i = 0; i < MTL::X256<MTL::I16>::Increment; i++)
      MTL_EQUAL(prod[i], (MTL::I16)20000);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::I16> a((MTL::I16)100);
    MTL::X512<MTL::I16> b((MTL::I16)200);
    MTL::X512<MTL::I16> prod = a * b;
    for (int i = 0; i < MTL::X512<MTL::I16>::Increment; i++)
      MTL_EQUAL(prod[i], (MTL::I16)20000);
  }
#endif
}

//
// F32 arithmetic
//
TEST(TestF32_Arithmetic)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(3.0f);
    MTL::X128<MTL::F32> b(2.0f);
    MTL::X128<MTL::F32> sum  = a + b;
    MTL::X128<MTL::F32> diff = a - b;
    MTL::X128<MTL::F32> prod = a * b;
    MTL::X128<MTL::F32> quot = a / b;
    MTL::X128<MTL::F32> neg  = -a;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0f);
      MTL_EQUAL(diff[i], 1.0f);
      MTL_EQUAL(prod[i], 6.0f);
      MTL_EQUAL(quot[i], 1.5f);
      MTL_EQUAL(neg[i], -3.0f);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(3.0f);
    MTL::X256<MTL::F32> b(2.0f);
    MTL::X256<MTL::F32> sum  = a + b;
    MTL::X256<MTL::F32> diff = a - b;
    MTL::X256<MTL::F32> prod = a * b;
    MTL::X256<MTL::F32> quot = a / b;
    MTL::X256<MTL::F32> neg  = -a;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0f);
      MTL_EQUAL(diff[i], 1.0f);
      MTL_EQUAL(prod[i], 6.0f);
      MTL_EQUAL(quot[i], 1.5f);
      MTL_EQUAL(neg[i], -3.0f);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(3.0f);
    MTL::X512<MTL::F32> b(2.0f);
    MTL::X512<MTL::F32> sum  = a + b;
    MTL::X512<MTL::F32> diff = a - b;
    MTL::X512<MTL::F32> prod = a * b;
    MTL::X512<MTL::F32> quot = a / b;
    MTL::X512<MTL::F32> neg  = -a;
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0f);
      MTL_EQUAL(diff[i], 1.0f);
      MTL_EQUAL(prod[i], 6.0f);
      MTL_EQUAL(quot[i], 1.5f);
      MTL_EQUAL(neg[i], -3.0f);
    }
  }
#endif
}

//
// F32 compound assignment operators (+=, -=, *=, /=)
//
TEST(TestF32_CompoundAssign)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(10.0f);
    MTL::X128<MTL::F32> b(3.0f);
    a += b;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 13.0f);
    a -= b;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
    a *= b;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 30.0f);
    a /= b;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(10.0f);
    MTL::X256<MTL::F32> b(3.0f);
    a += b;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 13.0f);
    a -= b;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
    a *= b;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 30.0f);
    a /= b;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(10.0f);
    MTL::X512<MTL::F32> b(3.0f);
    a += b;
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 13.0f);
    a -= b;
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
    a *= b;
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 30.0f);
    a /= b;
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(a[i], 10.0f);
  }
#endif
}

//
// F64 arithmetic
//
TEST(TestF64_Arithmetic)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F64> a(3.0);
    MTL::X128<MTL::F64> b(2.0);
    MTL::X128<MTL::F64> sum  = a + b;
    MTL::X128<MTL::F64> diff = a - b;
    MTL::X128<MTL::F64> prod = a * b;
    MTL::X128<MTL::F64> quot = a / b;
    MTL::X128<MTL::F64> neg  = -a;
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0);
      MTL_EQUAL(diff[i], 1.0);
      MTL_EQUAL(prod[i], 6.0);
      MTL_EQUAL(quot[i], 1.5);
      MTL_EQUAL(neg[i], -3.0);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F64> a(3.0);
    MTL::X256<MTL::F64> b(2.0);
    MTL::X256<MTL::F64> sum  = a + b;
    MTL::X256<MTL::F64> diff = a - b;
    MTL::X256<MTL::F64> prod = a * b;
    MTL::X256<MTL::F64> quot = a / b;
    MTL::X256<MTL::F64> neg  = -a;
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0);
      MTL_EQUAL(diff[i], 1.0);
      MTL_EQUAL(prod[i], 6.0);
      MTL_EQUAL(quot[i], 1.5);
      MTL_EQUAL(neg[i], -3.0);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F64> a(3.0);
    MTL::X512<MTL::F64> b(2.0);
    MTL::X512<MTL::F64> sum  = a + b;
    MTL::X512<MTL::F64> diff = a - b;
    MTL::X512<MTL::F64> prod = a * b;
    MTL::X512<MTL::F64> quot = a / b;
    MTL::X512<MTL::F64> neg  = -a;
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  5.0);
      MTL_EQUAL(diff[i], 1.0);
      MTL_EQUAL(prod[i], 6.0);
      MTL_EQUAL(quot[i], 1.5);
      MTL_EQUAL(neg[i], -3.0);
    }
  }
#endif
}

//
// F32 SquareRoot
//
TEST(TestF32_SquareRoot)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(16.0f);
    MTL::X128<MTL::F32> s = a.SquareRoot();
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(s[i], 4.0f);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(16.0f);
    MTL::X256<MTL::F32> s = a.SquareRoot();
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(s[i], 4.0f);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(16.0f);
    MTL::X512<MTL::F32> s = a.SquareRoot();
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(s[i], 4.0f);
  }
#endif
}

//
// F64 SquareRoot
//
TEST(TestF64_SquareRoot)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F64> a(25.0);
    MTL::X128<MTL::F64> s = a.SquareRoot();
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
      MTL_EQUAL(s[i], 5.0);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F64> a(25.0);
    MTL::X256<MTL::F64> s = a.SquareRoot();
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
      MTL_EQUAL(s[i], 5.0);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F64> a(25.0);
    MTL::X512<MTL::F64> s = a.SquareRoot();
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
      MTL_EQUAL(s[i], 5.0);
  }
#endif
}

//
// F32 comparisons
//
TEST(TestF32_Comparisons)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(3.0f);
    MTL::X128<MTL::F32> b(5.0f);
    MTL::X128<MTL::F32> c(3.0f);
    MTL::X128<MTL::F32> lt = a < b;
    MTL::X128<MTL::F32> gt = b > a;
    MTL::X128<MTL::F32> eq = a == c;
    MTL::X128<MTL::F32> le = a <= c;
    MTL::X128<MTL::F32> ge = a >= c;
    // True comparison results have all bits set (NaN pattern).
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_VERIFY(lt[i] != 0.0f);
      MTL_VERIFY(gt[i] != 0.0f);
      MTL_VERIFY(eq[i] != 0.0f);
      MTL_VERIFY(le[i] != 0.0f);
      MTL_VERIFY(ge[i] != 0.0f);
    }
    MTL::X128<MTL::F32> notLt = b < a;
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(notLt[i], 0.0f);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(3.0f);
    MTL::X256<MTL::F32> b(5.0f);
    MTL::X256<MTL::F32> c(3.0f);
    MTL::X256<MTL::F32> lt = a < b;
    MTL::X256<MTL::F32> gt = b > a;
    MTL::X256<MTL::F32> eq = a == c;
    MTL::X256<MTL::F32> le = a <= c;
    MTL::X256<MTL::F32> ge = a >= c;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_VERIFY(lt[i] != 0.0f);
      MTL_VERIFY(gt[i] != 0.0f);
      MTL_VERIFY(eq[i] != 0.0f);
      MTL_VERIFY(le[i] != 0.0f);
      MTL_VERIFY(ge[i] != 0.0f);
    }
    MTL::X256<MTL::F32> notLt = b < a;
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(notLt[i], 0.0f);
  }
#endif
}

//
// F64 comparisons
//
TEST(TestF64_Comparisons)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F64> a(3.0);
    MTL::X128<MTL::F64> b(5.0);
    MTL::X128<MTL::F64> eq = a == a;
    MTL::X128<MTL::F64> lt = a < b;
    MTL::X128<MTL::F64> gt = b > a;
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
    {
      MTL_VERIFY(eq[i] != 0.0);
      MTL_VERIFY(lt[i] != 0.0);
      MTL_VERIFY(gt[i] != 0.0);
    }
    MTL::X128<MTL::F64> notLt = b < a;
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
      MTL_EQUAL(notLt[i], 0.0);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F64> a(3.0);
    MTL::X256<MTL::F64> b(5.0);
    MTL::X256<MTL::F64> eq = a == a;
    MTL::X256<MTL::F64> lt = a < b;
    MTL::X256<MTL::F64> gt = b > a;
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
    {
      MTL_VERIFY(eq[i] != 0.0);
      MTL_VERIFY(lt[i] != 0.0);
      MTL_VERIFY(gt[i] != 0.0);
    }
    MTL::X256<MTL::F64> notLt = b < a;
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
      MTL_EQUAL(notLt[i], 0.0);
  }
#endif
}

//
// Min / Max
//
TEST(TestMinMax)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(3.0f);
    MTL::X128<MTL::F32> b(5.0f);
    MTL::X128<MTL::F32> mn = MTL::Min(a, b);
    MTL::X128<MTL::F32> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0f);
      MTL_EQUAL(mx[i], 5.0f);
    }
  }
  {
    MTL::X128<MTL::F64> a(3.0);
    MTL::X128<MTL::F64> b(5.0);
    MTL::X128<MTL::F64> mn = MTL::Min(a, b);
    MTL::X128<MTL::F64> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0);
      MTL_EQUAL(mx[i], 5.0);
    }
  }
  {
    MTL::X128<MTL::U8> a((MTL::U8)10);
    MTL::X128<MTL::U8> b((MTL::U8)20);
    MTL::X128<MTL::U8> mn = MTL::Min(a, b);
    MTL::X128<MTL::U8> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X128<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::U8)10);
      MTL_EQUAL(mx[i], (MTL::U8)20);
    }
  }
  {
    MTL::X128<MTL::I16> a((MTL::I16)-5);
    MTL::X128<MTL::I16> b((MTL::I16)5);
    MTL::X128<MTL::I16> mn = MTL::Min(a, b);
    MTL::X128<MTL::I16> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X128<MTL::I16>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::I16)-5);
      MTL_EQUAL(mx[i], (MTL::I16)5);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(3.0f);
    MTL::X256<MTL::F32> b(5.0f);
    MTL::X256<MTL::F32> mn = MTL::Min(a, b);
    MTL::X256<MTL::F32> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0f);
      MTL_EQUAL(mx[i], 5.0f);
    }
  }
  {
    MTL::X256<MTL::F64> a(3.0);
    MTL::X256<MTL::F64> b(5.0);
    MTL::X256<MTL::F64> mn = MTL::Min(a, b);
    MTL::X256<MTL::F64> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0);
      MTL_EQUAL(mx[i], 5.0);
    }
  }
  {
    MTL::X256<MTL::U8> a((MTL::U8)10);
    MTL::X256<MTL::U8> b((MTL::U8)20);
    MTL::X256<MTL::U8> mn = MTL::Min(a, b);
    MTL::X256<MTL::U8> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X256<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::U8)10);
      MTL_EQUAL(mx[i], (MTL::U8)20);
    }
  }
  {
    MTL::X256<MTL::I16> a((MTL::I16)-5);
    MTL::X256<MTL::I16> b((MTL::I16)5);
    MTL::X256<MTL::I16> mn = MTL::Min(a, b);
    MTL::X256<MTL::I16> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X256<MTL::I16>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::I16)-5);
      MTL_EQUAL(mx[i], (MTL::I16)5);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(3.0f);
    MTL::X512<MTL::F32> b(5.0f);
    MTL::X512<MTL::F32> mn = MTL::Min(a, b);
    MTL::X512<MTL::F32> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0f);
      MTL_EQUAL(mx[i], 5.0f);
    }
  }
  {
    MTL::X512<MTL::F64> a(3.0);
    MTL::X512<MTL::F64> b(5.0);
    MTL::X512<MTL::F64> mn = MTL::Min(a, b);
    MTL::X512<MTL::F64> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(mn[i], 3.0);
      MTL_EQUAL(mx[i], 5.0);
    }
  }
  {
    MTL::X512<MTL::U8> a((MTL::U8)10);
    MTL::X512<MTL::U8> b((MTL::U8)20);
    MTL::X512<MTL::U8> mn = MTL::Min(a, b);
    MTL::X512<MTL::U8> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X512<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::U8)10);
      MTL_EQUAL(mx[i], (MTL::U8)20);
    }
  }
  {
    MTL::X512<MTL::I16> a((MTL::I16)-5);
    MTL::X512<MTL::I16> b((MTL::I16)5);
    MTL::X512<MTL::I16> mn = MTL::Min(a, b);
    MTL::X512<MTL::I16> mx = MTL::Max(a, b);
    for (int i = 0; i < MTL::X512<MTL::I16>::Increment; i++)
    {
      MTL_EQUAL(mn[i], (MTL::I16)-5);
      MTL_EQUAL(mx[i], (MTL::I16)5);
    }
  }
#endif
}

//
// Abs
//
TEST(TestAbs)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(-7.5f);
    MTL::X128<MTL::F32> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 7.5f);
  }
  {
    MTL::X128<MTL::F64> a(-7.5);
    MTL::X128<MTL::F64> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 7.5);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(-7.5f);
    MTL::X256<MTL::F32> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 7.5f);
  }
  {
    MTL::X256<MTL::F64> a(-7.5);
    MTL::X256<MTL::F64> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 7.5);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(-7.5f);
    MTL::X512<MTL::F32> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 7.5f);
  }
  {
    MTL::X512<MTL::F64> a(-7.5);
    MTL::X512<MTL::F64> r = MTL::Abs(a);
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 7.5);
  }
#endif
}

//
// MultiplyAndAdd (FMA)
//
TEST(TestMultiplyAndAdd)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(2.0f);
    MTL::X128<MTL::F32> b(3.0f);
    MTL::X128<MTL::F32> c(10.0f);
    MTL::X128<MTL::F32> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 16.0f);
  }
  {
    MTL::X128<MTL::F64> a(2.0);
    MTL::X128<MTL::F64> b(3.0);
    MTL::X128<MTL::F64> c(10.0);
    MTL::X128<MTL::F64> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 16.0);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(2.0f);
    MTL::X256<MTL::F32> b(3.0f);
    MTL::X256<MTL::F32> c(10.0f);
    MTL::X256<MTL::F32> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 16.0f);
  }
  {
    MTL::X256<MTL::F64> a(2.0);
    MTL::X256<MTL::F64> b(3.0);
    MTL::X256<MTL::F64> c(10.0);
    MTL::X256<MTL::F64> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 16.0);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(2.0f);
    MTL::X512<MTL::F32> b(3.0f);
    MTL::X512<MTL::F32> c(10.0f);
    MTL::X512<MTL::F32> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
      MTL_EQUAL(r[i], 16.0f);
  }
  {
    MTL::X512<MTL::F64> a(2.0);
    MTL::X512<MTL::F64> b(3.0);
    MTL::X512<MTL::F64> c(10.0);
    MTL::X512<MTL::F64> r = MTL::MultiplyAndAdd(a, b, c);
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
      MTL_EQUAL(r[i], 16.0);
  }
#endif
}

//
// Conditional
//
TEST(TestConditional)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(1.0f);
    MTL::X128<MTL::F32> b(2.0f);
    MTL::X128<MTL::F32> condTrue  = a == a;  // All bits set.
    MTL::X128<MTL::F32> condFalse = a < a;   // All bits zero.
    MTL::X128<MTL::F32> r1 = MTL::Conditional(condTrue,  a, b);
    MTL::X128<MTL::F32> r2 = MTL::Conditional(condFalse, a, b);
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(r1[i], 1.0f);
      MTL_EQUAL(r2[i], 2.0f);
    }
  }
  {
    MTL::X128<MTL::F64> a(1.0);
    MTL::X128<MTL::F64> b(2.0);
    MTL::X128<MTL::F64> condTrue  = a == a;
    MTL::X128<MTL::F64> condFalse = a < a;
    MTL::X128<MTL::F64> r1 = MTL::Conditional(condTrue,  a, b);
    MTL::X128<MTL::F64> r2 = MTL::Conditional(condFalse, a, b);
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(r1[i], 1.0);
      MTL_EQUAL(r2[i], 2.0);
    }
  }
  {
    MTL::X128<MTL::I32> a((MTL::I32)10);
    MTL::X128<MTL::I32> b((MTL::I32)20);
    MTL::X128<MTL::I32> condTrue = a == a;
    MTL::X128<MTL::I32> r = MTL::Conditional(condTrue, a, b);
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(r[i], 10);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(1.0f);
    MTL::X256<MTL::F32> b(2.0f);
    MTL::X256<MTL::F32> condTrue  = a == a;
    MTL::X256<MTL::F32> condFalse = a < a;
    MTL::X256<MTL::F32> r1 = MTL::Conditional(condTrue,  a, b);
    MTL::X256<MTL::F32> r2 = MTL::Conditional(condFalse, a, b);
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(r1[i], 1.0f);
      MTL_EQUAL(r2[i], 2.0f);
    }
  }
  {
    MTL::X256<MTL::F64> a(1.0);
    MTL::X256<MTL::F64> b(2.0);
    MTL::X256<MTL::F64> condTrue  = a == a;
    MTL::X256<MTL::F64> condFalse = a < a;
    MTL::X256<MTL::F64> r1 = MTL::Conditional(condTrue,  a, b);
    MTL::X256<MTL::F64> r2 = MTL::Conditional(condFalse, a, b);
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(r1[i], 1.0);
      MTL_EQUAL(r2[i], 2.0);
    }
  }
  {
    MTL::X256<MTL::I32> a((MTL::I32)10);
    MTL::X256<MTL::I32> b((MTL::I32)20);
    MTL::X256<MTL::I32> condTrue = a == a;
    MTL::X256<MTL::I32> r = MTL::Conditional(condTrue, a, b);
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(r[i], 10);
  }
#endif
}

//
// Load / Store (unaligned)
//
TEST(TestLoadStore)
{
#if MTL_ENABLE_SSE
  {
    MTL::F32 src[4] = { 1.0f, 2.0f, 3.0f, 4.0f };
    MTL::F32 dst[4] = {};
    MTL::X128<MTL::F32> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 4; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
  {
    MTL::F64 src[2] = { 1.0, 2.0 };
    MTL::F64 dst[2] = {};
    MTL::X128<MTL::F64> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 2; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::F32 src[8] = { 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f };
    MTL::F32 dst[8] = {};
    MTL::X256<MTL::F32> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 8; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
  {
    MTL::F64 src[4] = { 1.0, 2.0, 3.0, 4.0 };
    MTL::F64 dst[4] = {};
    MTL::X256<MTL::F64> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 4; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
  {
    MTL::I8 src[32];
    MTL::I8 dst[32] = {};
    for (int i = 0; i < 32; i++) src[i] = (MTL::I8)(i + 1);
    MTL::X256<MTL::I8> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 32; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::F32 src[16];
    MTL::F32 dst[16] = {};
    for (int i = 0; i < 16; i++) src[i] = (MTL::F32)(i + 1);
    MTL::X512<MTL::F32> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 16; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
  {
    MTL::F64 src[8];
    MTL::F64 dst[8] = {};
    for (int i = 0; i < 8; i++) src[i] = (MTL::F64)(i + 1);
    MTL::X512<MTL::F64> x;
    x.Load(src);
    x.Store(dst);
    for (int i = 0; i < 8; i++)
      MTL_EQUAL(dst[i], src[i]);
  }
#endif
}

//
// Zeros() and Ones() static methods
//
TEST(TestZerosOnes)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> z = MTL::X128<MTL::F32>::Zeros();
    MTL::X128<MTL::F32> o = MTL::X128<MTL::F32>::Ones();
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0f);
      MTL_EQUAL(o[i], 1.0f);
    }
  }
  {
    MTL::X128<MTL::F64> z = MTL::X128<MTL::F64>::Zeros();
    MTL::X128<MTL::F64> o = MTL::X128<MTL::F64>::Ones();
    for (int i = 0; i < MTL::X128<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0);
      MTL_EQUAL(o[i], 1.0);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> z = MTL::X256<MTL::F32>::Zeros();
    MTL::X256<MTL::F32> o = MTL::X256<MTL::F32>::Ones();
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0f);
      MTL_EQUAL(o[i], 1.0f);
    }
  }
  {
    MTL::X256<MTL::F64> z = MTL::X256<MTL::F64>::Zeros();
    MTL::X256<MTL::F64> o = MTL::X256<MTL::F64>::Ones();
    for (int i = 0; i < MTL::X256<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0);
      MTL_EQUAL(o[i], 1.0);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> z = MTL::X512<MTL::F32>::Zeros();
    MTL::X512<MTL::F32> o = MTL::X512<MTL::F32>::Ones();
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0f);
      MTL_EQUAL(o[i], 1.0f);
    }
  }
  {
    MTL::X512<MTL::F64> z = MTL::X512<MTL::F64>::Zeros();
    MTL::X512<MTL::F64> o = MTL::X512<MTL::F64>::Ones();
    for (int i = 0; i < MTL::X512<MTL::F64>::Increment; i++)
    {
      MTL_EQUAL(z[i], 0.0);
      MTL_EQUAL(o[i], 1.0);
    }
  }
#endif
}

//
// F32 Reciprocal and ReciprocalPrecise (SSE and AVX only)
//
TEST(TestF32_Reciprocal)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(4.0f);
    MTL::X128<MTL::F32> r = a.Reciprocal();
    MTL::X128<MTL::F32> rp = a.ReciprocalPrecise();
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL_FLOAT(r[i], 0.25f, 0.002f);
      MTL_EQUAL_FLOAT(rp[i], 0.25f, 0.0001f);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(4.0f);
    MTL::X256<MTL::F32> r = a.Reciprocal();
    MTL::X256<MTL::F32> rp = a.ReciprocalPrecise();
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL_FLOAT(r[i], 0.25f, 0.002f);
      MTL_EQUAL_FLOAT(rp[i], 0.25f, 0.0001f);
    }
  }
#endif
}

//
// F32 ReciprocalSquareRoot and ReciprocalSquareRootPrecise (SSE and AVX only)
//
TEST(TestF32_ReciprocalSquareRoot)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(4.0f);
    MTL::X128<MTL::F32> r = a.ReciprocalSquareRoot();
    MTL::X128<MTL::F32> rp = a.ReciprocalSquareRootPrecise();
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL_FLOAT(r[i], 0.5f, 0.002f);
      MTL_EQUAL_FLOAT(rp[i], 0.5f, 0.0001f);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(4.0f);
    MTL::X256<MTL::F32> r = a.ReciprocalSquareRoot();
    MTL::X256<MTL::F32> rp = a.ReciprocalSquareRootPrecise();
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL_FLOAT(r[i], 0.5f, 0.002f);
      MTL_EQUAL_FLOAT(rp[i], 0.5f, 0.0001f);
    }
  }
#endif
}

//
// RoundedIntegers / TruncatedIntegers
//
TEST(TestF32_IntegerConversion)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> a(3.7f);
    MTL::X128<MTL::I32> rounded = a.RoundedIntegers();
    MTL::X128<MTL::I32> truncated = a.TruncatedIntegers();
    for (int i = 0; i < MTL::X128<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(rounded[i], 4);
      MTL_EQUAL(truncated[i], 3);
    }
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> a(3.7f);
    MTL::X256<MTL::I32> rounded = a.RoundedIntegers();
    MTL::X256<MTL::I32> truncated = a.TruncatedIntegers();
    for (int i = 0; i < MTL::X256<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(rounded[i], 4);
      MTL_EQUAL(truncated[i], 3);
    }
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> a(3.7f);
    MTL::X512<MTL::I32> rounded = a.RoundedIntegers();
    MTL::X512<MTL::I32> truncated = a.TruncatedIntegers();
    for (int i = 0; i < MTL::X512<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(rounded[i], 4);
      MTL_EQUAL(truncated[i], 3);
    }
  }
#endif
}

//
// I32 comparison operators (SSE/AVX only - AVX512 returns mask types)
//
TEST(TestI32_Comparisons)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::I32> a((MTL::I32)3);
    MTL::X128<MTL::I32> b((MTL::I32)5);
    MTL::X128<MTL::I32> eq = a == a;
    MTL::X128<MTL::I32> lt = a < b;
    MTL::X128<MTL::I32> gt = b > a;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(eq[i], -1);  // All bits set.
      MTL_EQUAL(lt[i], -1);
      MTL_EQUAL(gt[i], -1);
    }
    MTL::X128<MTL::I32> notEq = a == b;
    for (int i = 0; i < MTL::X128<MTL::I32>::Increment; i++)
      MTL_EQUAL(notEq[i], 0);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::I32> a((MTL::I32)3);
    MTL::X256<MTL::I32> b((MTL::I32)5);
    MTL::X256<MTL::I32> eq = a == a;
    MTL::X256<MTL::I32> lt = a < b;
    MTL::X256<MTL::I32> gt = b > a;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
    {
      MTL_EQUAL(eq[i], -1);
      MTL_EQUAL(lt[i], -1);
      MTL_EQUAL(gt[i], -1);
    }
    MTL::X256<MTL::I32> notEq = a == b;
    for (int i = 0; i < MTL::X256<MTL::I32>::Increment; i++)
      MTL_EQUAL(notEq[i], 0);
  }
#endif
}

//
// Multi-value Set constructors
//
TEST(TestMultiValueConstructors)
{
#if MTL_ENABLE_SSE
  {
    MTL::X128<MTL::F32> x(1.0f, 2.0f, 3.0f, 4.0f);
    MTL_EQUAL(x[0], 1.0f);
    MTL_EQUAL(x[1], 2.0f);
    MTL_EQUAL(x[2], 3.0f);
    MTL_EQUAL(x[3], 4.0f);
  }
  {
    MTL::X128<MTL::F64> x(1.0, 2.0);
    MTL_EQUAL(x[0], 1.0);
    MTL_EQUAL(x[1], 2.0);
  }
#endif
#if MTL_ENABLE_AVX
  {
    MTL::X256<MTL::F32> x(1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f);
    for (int i = 0; i < 8; i++)
      MTL_EQUAL(x[i], (float)(i + 1));
  }
  {
    MTL::X256<MTL::F64> x(1.0, 2.0, 3.0, 4.0);
    for (int i = 0; i < 4; i++)
      MTL_EQUAL(x[i], (double)(i + 1));
  }
#endif
#if MTL_ENABLE_AVX512
  {
    MTL::X512<MTL::F32> x(1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f,
                           9.f, 10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f);
    for (int i = 0; i < 16; i++)
      MTL_EQUAL(x[i], (float)(i + 1));
  }
  {
    MTL::X512<MTL::F64> x(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);
    for (int i = 0; i < 8; i++)
      MTL_EQUAL(x[i], (double)(i + 1));
  }
#endif
}

//
// XX<T> alias tests (uses the widest available SIMD)
//
TEST(TestXX_Alias)
{
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  {
    MTL::XX<MTL::F32> a(5.0f);
    MTL::XX<MTL::F32> b(3.0f);
    MTL::XX<MTL::F32> sum = a + b;
    MTL::XX<MTL::F32> diff = a - b;
    MTL::XX<MTL::F32> prod = a * b;
    for (int i = 0; i < MTL::XX<MTL::F32>::Increment; i++)
    {
      MTL_EQUAL(sum[i],  8.0f);
      MTL_EQUAL(diff[i], 2.0f);
      MTL_EQUAL(prod[i], 15.0f);
    }
  }
  {
    MTL::XX<MTL::F64> a(5.0);
    MTL::XX<MTL::F64> b(3.0);
    MTL::XX<MTL::F64> sum = a + b;
    for (int i = 0; i < MTL::XX<MTL::F64>::Increment; i++)
      MTL_EQUAL(sum[i], 8.0);
  }
  {
    MTL::XX<MTL::I8> a((MTL::I8)10);
    MTL::XX<MTL::I8> b((MTL::I8)3);
    MTL::XX<MTL::I8> sum = a + b;
    MTL::XX<MTL::I8> diff = a - b;
    for (int i = 0; i < MTL::XX<MTL::I8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::I8)13);
      MTL_EQUAL(diff[i], (MTL::I8)7);
    }
  }
  {
    MTL::XX<MTL::U8> a((MTL::U8)200);
    MTL::XX<MTL::U8> b((MTL::U8)50);
    MTL::XX<MTL::U8> sum = a + b;
    MTL::XX<MTL::U8> diff = a - b;
    for (int i = 0; i < MTL::XX<MTL::U8>::Increment; i++)
    {
      MTL_EQUAL(sum[i], (MTL::U8)250);
      MTL_EQUAL(diff[i], (MTL::U8)150);
    }
  }
#endif
}
