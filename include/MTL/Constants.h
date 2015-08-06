//
// Math Template Library
//
// Copyright (c) 2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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

#ifndef MTL_CONSTANTS_H
#define MTL_CONSTANTS_H

#include "Definitions.h"

namespace MTL
{

// Commonly used constants.
static const double kPi          = 3.141592653589793238460;
static const double kPiOverTwo   = 1.570796326794896619230;
static const double kPiOverFour  = 0.785398163397448309616;
static const double kTwoPi       = 2.0 * kPi;
static const double kOneOverPi   = 1./kPi;
static const double kTwoOverPi   = 2./kPi;
static const double kHalf        = 0.500000000000000000000;
static const double kOneThird    = 0.333333333333333333333;
static const double kTwoThirds   = 0.666666666666666666667;

// Angle conversion constants.
static const double kDegreesToRadians = kPi / 180.0;
static const double kRadiansToDegrees = 180.0 / kPi;

// Some helpers for floating point constants.
static const long kSign32  []   = { 0x80000000 };
static const long kNoSign32[]   = { 0x7FFFFFFF };
static const long kSign64  []   = { 0x00000000, 0x80000000 };
static const long kNoSign64[]   = { 0xFFFFFFFF, 0x7FFFFFFF };
static const long kInfinity64[] = { 0x00000000, 0x7FF00000 };

// Special floating point values.
static const double kNAN = (double&)*kNoSign64;
static const double kINF = (double&)*kInfinity64;

}  // namespace MTL


#endif // MTL_CONSTANTS_H
