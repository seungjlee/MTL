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

#include <MTL/Test.h>
#include <MTL/DynamicVectorOperators.h>
#include <MTL/Random.h>
#include <MTL/Utilities.h>

using namespace MTL;

static const double kTolF32 = 1e-6;
static const double kTolF64 = 1e-14;

TEST(TestClassNoDefaultConstructor)
{
  class X
  {
  public:
    X(int i) : test_(i) {}
    int test_;
  };

  DynamicVector<X> v(5, X(7));

  FOR_EACH_INDEX(v)
    MTL_EQUAL(v[vIndex].test_, 7);
}

TEST(TestSetAll)
{
  DynamicVector<double> v(100);
  v.SetAll(42);

  FOR_EACH_INDEX(v)
    MTL_EQUAL_FLOAT(v[vIndex], 42, kTolF64);
}

TEST(TestReductionsSmallVectorsF64)
{
  DynamicVector<F64> a1(1, 3);
  DynamicVector<F64> a2(1, -7);

  MTL_EQUAL_FLOAT(Sum(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 9, kTolF64);
  MTL_EQUAL_FLOAT(Min(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(Max(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 3, kTolF64);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 3, kTolF64);

  MTL_EQUAL_FLOAT(Sum(a2), -7, kTolF64);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 7, kTolF64);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 49, kTolF64);
  MTL_EQUAL_FLOAT(Min(a2), -7, kTolF64);
  MTL_EQUAL_FLOAT(Max(a2), -7, kTolF64);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 7, kTolF64);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 7, kTolF64);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 7, kTolF64);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -21, kTolF64);

  double aa1[] = { 3.0, 1.2,  4.5};
  double aa2[] = {-2.0, 5.0, -1.0};
  a1.Assign(aa1, 3);
  a2.Assign(aa2, 3);

  MTL_EQUAL_FLOAT(Sum(a1), 8.7, kTolF64);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 8.7, kTolF64);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 30.69, kTolF64);
  MTL_EQUAL_FLOAT(Min(a1), 1.2, kTolF64);
  MTL_EQUAL_FLOAT(Max(a1), 4.5, kTolF64);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 1.2, kTolF64);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 4.5, kTolF64);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 4.5, kTolF64);

  MTL_EQUAL_FLOAT(Sum(a2), 2, kTolF64);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 8, kTolF64);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 30, kTolF64);
  MTL_EQUAL_FLOAT(Min(a2), -2, kTolF64);
  MTL_EQUAL_FLOAT(Max(a2),  5, kTolF64);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 1, kTolF64);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 5, kTolF64);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 5, kTolF64);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -4.5, kTolF64);

  MTL_EQUAL_FLOAT(Mean(a1), 2.9, kTolF64);
  MTL_EQUAL_FLOAT(Variance(a1, Mean(a1)), 2.73, kTolF64);
  MTL_EQUAL_FLOAT(RMS(a1), Sqrt(30.69 /3), kTolF64);
  MTL_EQUAL_FLOAT(FrobeniusNorm(a1), Sqrt(30.69), kTolF64);
}

TEST(TestReductionsSmallVectorsF32)
{
  DynamicVector<F32> a1(1, 3);
  DynamicVector<F32> a2(1, -7);

  MTL_EQUAL_FLOAT(Sum(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 9, kTolF32);
  MTL_EQUAL_FLOAT(Min(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(Max(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 3, kTolF32);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 3, kTolF32);

  MTL_EQUAL_FLOAT(Sum(a2), -7, kTolF32);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 7, kTolF32);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 49, kTolF32);
  MTL_EQUAL_FLOAT(Min(a2), -7, kTolF32);
  MTL_EQUAL_FLOAT(Max(a2), -7, kTolF32);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 7, kTolF32);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 7, kTolF32);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 7, kTolF32);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -21, kTolF32);

  float aa1[] = { 3.0f, 1.2f,  4.5f};
  float aa2[] = {-2.0f, 5.0f, -1.0f};
  a1.Assign(aa1, 3);
  a2.Assign(aa2, 3);

  MTL_EQUAL_FLOAT(Sum(a1), 8.7, kTolF32);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 8.7, kTolF32);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 30.69, kTolF32);
  MTL_EQUAL_FLOAT(Min(a1), 1.2, kTolF32);
  MTL_EQUAL_FLOAT(Max(a1), 4.5, kTolF32);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 1.2, kTolF32);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 4.5, kTolF32);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 4.5, kTolF32);

  MTL_EQUAL_FLOAT(Sum(a2), 2, kTolF32);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 8, kTolF32);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 30, kTolF32);
  MTL_EQUAL_FLOAT(Min(a2), -2, kTolF32);
  MTL_EQUAL_FLOAT(Max(a2),  5, kTolF32);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 1, kTolF32);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 5, kTolF32);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 5, kTolF32);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -4.5, kTolF32);

  MTL_EQUAL_FLOAT(Mean(a1), 2.9, kTolF32);
  MTL_EQUAL_FLOAT(Variance(a1, Mean(a1)), 2.73, kTolF32);
  MTL_EQUAL_FLOAT(RMS(a1), Sqrt(30.69 / 3), kTolF32);
  MTL_EQUAL_FLOAT(FrobeniusNorm(a1), Sqrt(30.69), kTolF32);
}

TEST(TestHouseholderQR_Speed)
{
  static const double kTol = 1e-14;

  enum
  {
    N = 128*1024,
    kRepeats = 100
  };

  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    DynamicVector<F64> v = random.DynamicVector<F64>(N, -10, 10);

    DynamicVector<F64> v2 = v;
    SquareAll(v2);

    DynamicVector<F64> v1 = v2;
    SquareRootAll(v1);

    FOR_EACH_INDEX(v)
    {
      MTL_EQUAL_FLOAT(v2[vIndex], Square(v[vIndex]), kTol);
      MTL_EQUAL_FLOAT(v1[vIndex], Abs(v[vIndex]), kTol);
    }

    ShowProgressBar(double(i + 1) / kRepeats);
  }
  Out() << std::endl;
}

TEST(TestCastConstructor)
{
  enum
  {
    N = 64*1024*1024
  };

  Timer t;

  DynamicVector<F64> vF64(N);
  for (I32 i = 0; i < N; i++)
    vF64[i] = i+1;

  t.ResetAndStart();
  DynamicVector<I32> vI32 = DynamicVector<I32>(vF64);
  t.Stop();
  wprintf(L"  Conversion time: %.3f msecs.\n", t.Milliseconds());

  for (I32 i = 0; i < N; i++)
    MTL_EQUAL_FLOAT(vI32[i], vF64[i], NumericalEpsilon<F64>());
}

TEST(Test_AddBack_Insert)
{
  enum
  {
    kSize = 1000,
    kStep = 10
  };

  DynamicVector<I32> v;

  U32 oldSize = 0;
  U32 newSize = kStep;

  while (v.Size() <= kSize)
  {
    v.Resize(newSize);
    for (U32 i = oldSize; i < newSize; i++)
      v[i] = i;

    oldSize = newSize;
    newSize += kStep;
  }

  FOR_EACH_INDEX(v)
    MTL_EQUAL(v[vIndex], (int)vIndex);

  v.Insert(&v[11], 777);

  FOR_EACH_INDEX(v)
  {
    if (vIndex < 11)
      MTL_EQUAL(v[vIndex], (int)vIndex);
    else if(vIndex > 11)
      MTL_EQUAL(v[vIndex], (int)vIndex-1);
    else
      MTL_EQUAL(v[vIndex], 777);
  }
}

TEST(Test_SumOfSquaredDifferences)
{
  static const double kTol = 1e-14;

  enum
  {
    N = 128*1024,
    kRepeats = 777
  };

  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    DynamicVector<F64> v1 = random.DynamicVector<F64>(N, -1.0, 1.0);
    DynamicVector<F64> v2 = random.DynamicVector<F64>(N, -1.0, 1.0);

    F64 s1 = SumOfSquares(v1 - v2);
    F64 s2 = SumOfSquaredDifferences(v1, v2);

    MTL_EQUAL_FLOAT(s1, s2, kTol);

    ShowProgressBar(double(i + 1) / kRepeats);
  }
  Out() << std::endl;
}
