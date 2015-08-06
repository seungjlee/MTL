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


#ifndef MTL_TRIGONOMETRY_H
#define MTL_TRIGONOMETRY_H

#include "Math.h"
#include <cmath>

namespace MTL
{

// Large offset to simplify computation by only dealing with positive numbers when we compute the
// correct sign.
static const double kTruncateOffset = 1e15;

template <class T> MTL_INLINE static void ComputeCosineSine(T& c, T& s, const T& angle)
{
  // Figure out if we should compute sine from cosine or viceversa. Better pick the best way to
  // avoid losing precision.
  U64 test = U64((angle+kPiOverFour)*kTwoOverPi+kTruncateOffset);
  if (test & 1)
  {
    c = cos(angle);
    s = Sqrt(1-Pow<2>(c));

    assert(angle*kOneOverPi >= -T(kTruncateOffset));
    s = U64(angle*kOneOverPi+T(kTruncateOffset)) & 1 ? -s : s;
  }
  else
  {
    s = sin(angle);
    c = Sqrt(1-Pow<2>(s));

    assert((angle+kPiOverTwo)*kOneOverPi >= -T(kTruncateOffset));
    c = U64((angle+kPiOverTwo)*kOneOverPi+T(kTruncateOffset)) & 1 ? -c : c;
  }
}

}  // namespace MTL

#endif // MTL_TRIGONOMETRY_H
