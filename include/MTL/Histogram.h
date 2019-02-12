//
// Math Template Library
//
// Copyright (c) 2018: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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

    const T* pDataEnd = pData + size;
    for (; pData < pDataEnd; pData++)
      Bins_[*pData]++;
  }

  template <class MaskT>
  Histogram(const T* pData, const MaskT* pMask, U32 size)
  {
    NumberOfElements_ = 0;
    memset(Bins_, 0, sizeof(Bins_));

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
