//
// Math Template Library
//
// Copyright (c) 2017: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#include <cmath>
#include <limits>

TEST(TestShuffle)
{
#if MTL_ENABLE_SSE
  MTL::X128<MTL::F32> x1234(1.f,2.f,3.f,4.f);
  MTL::X128<MTL::F32> x5678(5.f,6.f,7.f,8.f);

  MTL::X128<MTL::F32> x2468 = MTL::Shuffle<1|3<<2|1<<4|3<<6>(x1234,x5678);
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
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  MTL::XX<MTL::I32> xx = MTL::XX<MTL::I32>(1) ^ MTL::XX<MTL::I32>(-1);

  for (int k = 0; k < MTL::XX<MTL::I32>::Increment; k++)
    MTL_EQUAL(xx[k], -2);
#endif
}

TEST(TestToFloat)
{
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  const int size = 1111111;

  std::vector<MTL::I32>  u(size);
  std::vector<MTL::F32> v32(size);
  std::vector<MTL::F64> v64(size);

  for (int i = 0; i < size; i++)
    u[i] = i;

  MTL::Convert_StreamUnaligned_Parallel(&v32[0], &u[0], size, 0);
  MTL::Convert_StreamUnaligned_Parallel(&v64[0], &u[0], size, 0);

  for (int i = 0; i < size; i++)
  {
    MTL_EQUAL(v32[i], (MTL::F32)u[i]);
    MTL_EQUAL(v64[i], (MTL::F64)u[i]);
  }
#endif
}

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX || MTL_ENABLE_AVX512
#define MTL_TEST_FLOOR_CEIL 1
#else
#define MTL_TEST_FLOOR_CEIL 0
#endif

#if MTL_TEST_FLOOR_CEIL
// Values that cover fractions of both signs, exact integers, signed zeros, magnitudes at and past
// the point where every representable value is already an integer, and the non-finite inputs. The
// count is a multiple of 16 so it tiles evenly for every register width and element size.
template <class T> static std::vector<T> FloorCeilTestValues()
{
  const T kInfinity = std::numeric_limits<T>::infinity();
  const T kNaN      = std::numeric_limits<T>::quiet_NaN();

  return std::vector<T>
  {
    T( 0.0),        T(-0.0),        T( 1.0),        T(-1.0),
    T( 0.25),       T(-0.25),       T( 0.5),        T(-0.5),
    T( 0.75),       T(-0.75),       T( 1.5),        T(-1.5),
    T( 2.5),        T(-2.5),        T( 3.75),       T(-3.75),
    T( 100.5),      T(-100.5),      T( 1000000.5),  T(-1000000.5),
    T( 8388607.5),  T(-8388607.5),  T( 8388608.0),  T(-8388608.0),
    T( 4.0),        T(-4.0),        T( 1e30),       T(-1e30),
    kInfinity,      -kInfinity,     kNaN,           -kNaN
  };
}

template <class T> static void CheckFloorCeil(T actualFloor, T actualCeil, T input)
{
  if (std::isnan(input))
  {
    MTL_VERIFY(std::isnan(actualFloor));
    MTL_VERIFY(std::isnan(actualCeil));
    return;
  }

  const T expectedFloor = std::floor(input);
  const T expectedCeil  = std::ceil(input);

  MTL_EQUAL(actualFloor, expectedFloor);
  MTL_EQUAL(actualCeil, expectedCeil);

  // Zeros compare equal regardless of sign so the sign bit needs its own check. Rounding never
  // creates a zero out of a non-zero input, so this also confirms floor(-0.5) is -0 and not +0.
  MTL_EQUAL(std::signbit(actualFloor), std::signbit(expectedFloor));
  MTL_EQUAL(std::signbit(actualCeil), std::signbit(expectedCeil));
}

template <class XT, class T> static void TestFloorCeilRegister(const std::vector<T>& values)
{
  const MTL::SizeType increment = XT::Increment;

  for (MTL::SizeType offset = 0; offset + increment <= values.size(); offset += increment)
  {
    XT x(&values[offset]);
    XT xFloor = MTL::Floor(x);
    XT xCeil  = MTL::Ceil(x);

    for (MTL::SizeType k = 0; k < increment; k++)
      CheckFloorCeil<T>(xFloor[k], xCeil[k], values[offset + k]);
  }
}
#endif  // MTL_TEST_FLOOR_CEIL

TEST(TestFloorAndCeil)
{
#if MTL_TEST_FLOOR_CEIL
  std::vector<MTL::F32> values32 = FloorCeilTestValues<MTL::F32>();
  std::vector<MTL::F64> values64 = FloorCeilTestValues<MTL::F64>();

  TestFloorCeilRegister< MTL::X128<MTL::F32> >(values32);
  TestFloorCeilRegister< MTL::X128<MTL::F64> >(values64);

#if MTL_ENABLE_AVX || MTL_ENABLE_AVX512
  TestFloorCeilRegister< MTL::X256<MTL::F32> >(values32);
  TestFloorCeilRegister< MTL::X256<MTL::F64> >(values64);
#endif

#if MTL_ENABLE_AVX512
  TestFloorCeilRegister< MTL::X512<MTL::F32> >(values32);
  TestFloorCeilRegister< MTL::X512<MTL::F64> >(values64);
#endif

  // The XX alias resolves to the widest enabled register, which is what generic stream code uses.
  TestFloorCeilRegister< MTL::XX<MTL::F32> >(values32);
  TestFloorCeilRegister< MTL::XX<MTL::F64> >(values64);
#endif
}

TEST(TestFloorAndCeilLanes)
{
#if MTL_TEST_FLOOR_CEIL
  // Each lane is rounded independently: check a hand-written set of four values in a known order.
  MTL::X128<MTL::F32> x(-1.5f, -0.5f, 0.5f, 1.5f);

  MTL::X128<MTL::F32> xFloor = MTL::Floor(x);
  MTL_EQUAL(xFloor[0], -2.f);
  MTL_EQUAL(xFloor[1], -1.f);
  MTL_EQUAL(xFloor[2],  0.f);
  MTL_EQUAL(xFloor[3],  1.f);

  MTL::X128<MTL::F32> xCeil = MTL::Ceil(x);
  MTL_EQUAL(xCeil[0], -1.f);
  MTL_EQUAL(xCeil[1], -0.f);
  MTL_EQUAL(xCeil[2],  1.f);
  MTL_EQUAL(xCeil[3],  2.f);

  MTL::X128<MTL::F64> y(-2.25, 2.25);
  MTL_EQUAL(MTL::Floor(y)[0], -3.0);
  MTL_EQUAL(MTL::Floor(y)[1],  2.0);
  MTL_EQUAL(MTL::Ceil(y)[0],  -2.0);
  MTL_EQUAL(MTL::Ceil(y)[1],   3.0);
#endif
}
