//
// Math Template Library
//
// Copyright (c) 2018: Seung Jae Lee, https://github.com/seungjlee/MTL
//
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

#ifndef MTL_HISTOGRAM_H
#define MTL_HISTOGRAM_H

#include "Math.h"
#include <MTL/CPU.h>
#include <MTL/OpenMP.h>
#include <vector>

#ifndef MTL_HISTOGRAM_PARALLEL_THRESHOLD
#define MTL_HISTOGRAM_PARALLEL_THRESHOLD 65536
#endif

namespace MTL
{
// Only designed to be used by 16-bit and 8-bit unsigned integers.
template <class T>
class Histogram
{
  static_assert(std::is_same<T,U8>::value || std::is_same<T,U16>::value,
                "Only 16-bit and 8-bit unsigned integer types supported!");

  static const int BINS = std::is_same<T,U8>::value ? 256 : 65536;

public:
  Histogram(const T* pData, U32 size)
  {
    NumberOfElements_ = size;
    memset(Bins_, 0, sizeof(Bins_));

#if MTL_ENABLE_OPENMP
    int numberOfThreads = (int)MTL::CPU::Instance().NumberOfThreads();
    if (numberOfThreads > 1 && size >= MTL_HISTOGRAM_PARALLEL_THRESHOLD)
    {
      // Per-thread partial histograms reduced into Bins_ at the end. Local
      // bins are heap-allocated because BINS=65536 (256 KiB for U16) is too
      // large for a typical thread stack.
      #pragma omp parallel num_threads(numberOfThreads)
      {
        std::vector<U32> localBins(BINS, 0);
        #pragma omp for nowait
        for (I32 i = 0; i < (I32)size; i++)
          localBins[pData[i]]++;
        #pragma omp critical
        for (int b = 0; b < BINS; b++)
          Bins_[b] += localBins[b];
      }
      return;
    }
#endif

    const T* pDataEnd = pData + size;
    for (; pData < pDataEnd; pData++)
      Bins_[*pData]++;
  }

  template <class MaskT>
  Histogram(const T* pData, const MaskT* pMask, U32 size)
  {
    NumberOfElements_ = 0;
    memset(Bins_, 0, sizeof(Bins_));

#if MTL_ENABLE_OPENMP
    int numberOfThreads = (int)MTL::CPU::Instance().NumberOfThreads();
    if (numberOfThreads > 1 && size >= MTL_HISTOGRAM_PARALLEL_THRESHOLD)
    {
      U32 localCounts = 0;
      #pragma omp parallel num_threads(numberOfThreads) reduction(+:localCounts)
      {
        std::vector<U32> localBins(BINS, 0);
        #pragma omp for nowait
        for (I32 i = 0; i < (I32)size; i++)
        {
          if (pMask[i])
          {
            localBins[pData[i]]++;
            localCounts++;
          }
        }
        #pragma omp critical
        for (int b = 0; b < BINS; b++)
          Bins_[b] += localBins[b];
      }
      NumberOfElements_ = localCounts;
      return;
    }
#endif

    const T* pDataEnd = pData + size;
    for (; pData < pDataEnd; pData++, pMask++)
    {
      if (*pMask)
      {
        Bins_[*pData]++;
        NumberOfElements_++;
      }
    }
  }

  const U32& operator[](U32 index) const  { return Bins_[index]; }
        U32& operator[](U32 index)        { return Bins_[index]; }

  T Peak() const
  {
    T peak = 0;
    U32 peakCount = Bins_[0];

    for (U32 i = 1; i < BINS; i++)
    {
      if (Bins_[i] > peakCount)
      {
        peak = i;
        peakCount = Bins_[i];
      }
    }

    return peak;
  }
  
  T LevelAtPercentage100(F64 percentage0to100) const
  {
    return LevelAtPercentage(percentage0to100 * 0.01);
  }
  T LevelAtPercentage(F64 percentage0to1) const
  {
    if (NumberOfElements_ > 0)
    {
      // 1-based position of the percentile.
      U32 position = U32((NumberOfElements_ + 1) * percentage0to1);
      position = Max(position, 1U);
      position = Min(position, NumberOfElements_);

      for (U32 i = 0, count = 0; i < BINS; i++)
      {
        count += Bins_[i];
        if (count >= position)
          return i;
      }
      return BINS-1;
    }

    return 0;
  }

  F64 Mean() const
  {
    F64 sum = 0;
    for (U32 i = 1; i < BINS; i++)
      sum += i * Bins_[i];

    return sum / NumberOfElements_;
  }

  U32 NumberOfValuesGreaterThan(U32 threshold) const
  {
    U32 sum = 0;
    for (U32 i = threshold+1; i < BINS; i++)
      sum += Bins_[i];

    return sum;
  }

private:
  U32 Bins_[BINS];
  U32 NumberOfElements_;
};

typedef Histogram<U8> Histogram8;
typedef Histogram<U16> Histogram16;

}  // namespace MTL


#endif  // MTL_HISTOGRAM_H
