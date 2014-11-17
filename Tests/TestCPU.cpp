//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee
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
#include <MTL/CPU.h>
#include <MTL/DynamicVectorOperators.h>

using namespace MTL;

static const double kTol = 1e-13;

#define CHECK_FEATURE(FEATURE) \
  printf("  %-8s %s\n", #FEATURE, CPU::Instance().FEATURE()  ? "Yes" : "No");

TEST(TestCPU)
{
  CHECK_FEATURE(SSE);
  CHECK_FEATURE(SSE2);
  CHECK_FEATURE(AVX);
  CHECK_FEATURE(AVX2);
  CHECK_FEATURE(FMA);
}

TEST(TestMemoryBandwitdh)
{
  static const long kVectorSize = 4*1024*1024;
  static const long kTries = 5;

  Timer t;
  double bestTime;

  U64 maxNumberOfThreads = CPU::Instance().NumberOfThreads();

  // Used for overwriting the cache since we are trying to just test memory bandwidth.
  DynamicVector<double> tempV(kVectorSize);
  tempV.Zeros();

  for (long numberOfThreads = 1; numberOfThreads <= maxNumberOfThreads; numberOfThreads++)
  {
    DynamicVector<double> testV1(kVectorSize, 11);
    DynamicVector<double> testV2(kVectorSize, 22);
    double* p1 = testV1.Begin();
    double* p2 = testV2.Begin();

    printf("Number of Threads: %d\n", numberOfThreads);
    CPU::Instance().NumberOfThreads(numberOfThreads);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2.Zeros();
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Zero:                  %8.3f GB/s.\n",
           sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 0, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2.SetAll(101);
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Assign:                %8.3f GB/s.\n",
           sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 101, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 = testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Copy:                  %8.3f GB/s.\n",
           2 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 11, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 += 17;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Add Scalar:            %8.3f GB/s.\n",
           2 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 28, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 *= 7;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Multiply Scalar:       %8.3f GB/s.\n",
           2 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 77, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 /= 2;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Divide Scalar:         %8.3f GB/s.\n",
           2 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 5.5, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV1.SetAll(-123.456);
      testV2.SetAll(9999.777);
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 += testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Add Vectors:           %8.3f GB/s.\n",
           3 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 9999.777-123.456, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV1.SetAll(-123.456);
      testV2.SetAll(9999.777);
      tempV += numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 *= testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    printf("  Multiply Vectors:      %8.3f GB/s.\n",
           3 * sizeof(double) * kVectorSize * 1e-9 / bestTime);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], -123.456*9999.777, kTol);
  }
}
