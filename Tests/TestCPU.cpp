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

#include <MTL/Tools/Test.h>
#include <MTL/CPU.h>
#include <MTL/Math/DynamicVectorOperators.h>

using namespace MTL;

static const double kTol = 1e-13;

constexpr int MAX_TEST_THREADS = VALUE_DEBUG_RELEASE(2, 24);

#define CHECK_FEATURE(FEATURE) \
  wprintf(L"  %-8hs %hs\n", #FEATURE, CPU::Instance().FEATURE()  ? "Yes" : "No");

TEST(TestCPU)
{
  CHECK_FEATURE(SSE);
  CHECK_FEATURE(SSE2);
  CHECK_FEATURE(AVX);
  CHECK_FEATURE(AVX2);
  CHECK_FEATURE(FMA);
  CHECK_FEATURE(AVX512F);
}

TEST(TestMemoryBandwitdh)
{
  static const long kVectorSize = VALUE_DEBUG_RELEASE(2*1024*1024, 8*1024*1024);
  static const long kTries = 8;

  Timer t;
  double bestTime;

  int maxNumberOfThreads = Min(MAX_TEST_THREADS, (int)CPU::Instance().NumberOfThreads() * 2);

  // Used for overwriting the cache since we are trying to just test memory bandwidth.
  DynamicVector<double> tempV(kVectorSize);
  tempV.Zeros();

  double maxBandwidth = 0;
  int bestThreads = 0;

  for (int numberOfThreads = 2; numberOfThreads <= maxNumberOfThreads; numberOfThreads += 2)
  {
    DynamicVector<double> testV1(kVectorSize, 11);
    DynamicVector<double> testV2(kVectorSize, 22);

    wprintf(L"Number of Threads: %ld\n", numberOfThreads);
    CPU::Instance().NumberOfThreads(numberOfThreads);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2.Zeros();
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    double bandwidth = sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Zero:                  %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 0, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2.SetAll(101);
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Assign:                %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 101, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 += 17;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Add Scalar:            %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 28, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 *= 7;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Multiply Scalar:       %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 77, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV2 = testV1;
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 /= 2;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Divide Scalar:         %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 5.5, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 = testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = 2 * sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Copy:                  %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 11, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV1.SetAll(-123.456);
      testV2.SetAll(9999.777);
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 += testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = 2 * sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Add Vectors:           %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], 9999.777-123.456, kTol);

    bestTime = kINF;
    for (long i = 0; i < kTries; i++)
    {
      testV1.SetAll(-123.456);
      testV2.SetAll(9999.777);
      tempV += (double)numberOfThreads;  // Dirty the cache memory.

      t.ResetAndStart();
      testV2 *= testV1;
      t.Stop();
      bestTime = Min(bestTime, t.Seconds());
    }
    bandwidth = 2 * sizeof(double) * kVectorSize * 1e-9 / bestTime;
    if (bandwidth > maxBandwidth)
    {
      maxBandwidth = bandwidth;
      bestThreads = numberOfThreads;
    }
    wprintf(L"  Multiply Vectors:      %8.3f GB/s.\n", bandwidth);
    FOR_EACH_INDEX(testV2)
      MTL_EQUAL_FLOAT(testV2[testV2Index], -123.456*9999.777, kTol);
  }

  wprintf(COLOR_BOLD L"\nMax. Bandwidth: %.3f GB/s, %d threads.\n\n" COLOR_RESET, maxBandwidth, bestThreads);
}

TEST(TestStreamPerformance)
{
  static const long kVectorSize = 32*1024;
  static const long kIterations = 1000;
  static const long kTries = 20;

  Timer t;
  double bestTime;

  int maxNumberOfThreads = MAX_TEST_THREADS; // Min(MAX_TEST_THREADS, (int)CPU::Instance().NumberOfThreads());

  for (int numberOfThreads = 2; numberOfThreads <= maxNumberOfThreads; numberOfThreads *= 2)
  {
    DynamicVector<double> testV1(kVectorSize, -1);
    DynamicVector<double> testV2(kVectorSize, -2);

    wprintf(L"Number of Threads: %ld\n", numberOfThreads);
    CPU::Instance().NumberOfThreads(numberOfThreads);

    double dotProduct;
    bestTime = kINF;
    for (long k = 0; k < kTries; k++)
    {
      t.ResetAndStart();
      for (long i = 0; i < kIterations; i++)
        dotProduct = DotProduct(testV1, testV2);
      t.Stop();
      bestTime = Min(bestTime, t.Seconds() / kIterations);
    }
    wprintf(L"  Dot product:    %8.3f GFLOPS, %.6f msecs.\n",
            2 * kVectorSize * 1e-9 / bestTime, bestTime * 1e3);
    MTL_EQUAL_FLOAT(dotProduct, 2*kVectorSize, kTol);

    double sumOfSquares;
    bestTime = kINF;
    for (long k = 0; k < kTries; k++)
    {
      t.ResetAndStart();
      for (long i = 0; i < kIterations; i++)
        sumOfSquares = SumOfSquares(testV1);
      t.Stop();
      bestTime = Min(bestTime, t.Seconds() / kIterations);
    }
    wprintf(L"  Sum of Squares: %8.3f GFLOPS, %.6f msecs.\n",
            2 * kVectorSize * 1e-9 / bestTime, bestTime * 1e3);
    MTL_EQUAL_FLOAT(sumOfSquares, kVectorSize, kTol);
  }
}
