//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_LARGE_INTEGER_H
#define MTL_LARGE_INTEGER_H

#include "DynamicVector.h"
#include "Exception.h"

namespace MTL
{

class LargeInteger
{
public:
  LargeInteger(U64 number)
  {
    if (number < kLimit)
    {
      Data_.Resize(1);
      Data_[0] = U32(number);
    }
    else if (number < Square(kLimit))
    {
      Data_.Resize(2);
      Data_[1] = U32(number / kLimit);
      Data_[0] = U32(number - kLimit * Data_[1]);
    }
    else
    {
      throw Exception(L"Currently unsupported!");
    }
  }

  LargeInteger operator+(const LargeInteger& rhs) const
  {
    if (*this < rhs)
      return rhs + *this;

    LargeInteger sum = *this;

    for (U32 i = 0; i < rhs.Data_.Size(); i++)
    {
      U64 sum64 = U64(Data_[i]) + U64(rhs.Data_[i]);

      if (sum64 > kLimit)
      {
        sum64 -= kLimit;

        if (i + 1 == sum.Data_.Size())
        {
          sum.Data_.PushBack(1);
        }
        else
        {
          sum.Data_[i+1]++;
        }
      }

      sum.Data_[i] = U32(sum64);
    }

    return sum;
  }

  LargeInteger operator*(const LargeInteger& rhs) const
  {
    if (*this < rhs)
      return rhs * *this;

    LargeInteger product(0);
    product.Data_.Resize(Data_.Size());
    product.Data_.Zeros();

    for (U32 i = 0; i < rhs.Data_.Size(); i++)
    {
      U64 product64 = U64(Data_[i]) * U64(rhs.Data_[i]);

    }

    return product;
  }

  bool operator<(const LargeInteger& rhs) const
  {
    if (Data_.Size() == rhs.Data_.Size())
    {
      for (int i = (int)Data_.Size() - 1; i >= 0; i--)
      {
        if (Data_[i] < rhs.Data_[i])
          return true;

        if (Data_[i] > rhs.Data_[i])
          return false;
      }

      return false;
    }

    return Data_.Size() < rhs.Data_.Size();
  }

  String GetString() const
  {
    String str;
    wchar_t buffer[16];
    wsprintf(buffer, L"%d", Data_.Back());
    str = buffer;

    for (I32 i = I32(Data_.Size()) - 2; i >= 0; i--)
    {
      wsprintf(buffer, L"%09d", Data_[i]);
      str += buffer;
    }

    return str;
  }

private:
  static const U64 kLimit = 1000000000LL;

  DynamicVector<U32> Data_;
};

}  // namespace MTL

#endif // MTL_LARGE_INTEGER_H
