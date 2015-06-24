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

#ifndef MTL_MATH_H
#define MTL_MATH_H

#include "Definitions.h"
#include <float.h>

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


template <class T> MTL_INLINE static T Conditional(bool condition, const T& a, const T& b)
{
  return condition ? a : b;
}

template <class T> MTL_INLINE static T Max(const T& a, const T& b)
{
  return Conditional(a > b, a, b);
}
template <class T> MTL_INLINE static T Max(const T& a, const T& b, const T& c)
{
  return Max(Max(a,b),c);
}
template <class T> MTL_INLINE static T Max(const T& a, const T& b, const T& c, const T& d)
{
  return Max(Max(a,b),Max(c,d));
}

template <class T> MTL_INLINE static T Min(const T& a, const T& b)
{
  return Conditional(a < b, a, b);
}
template <class T> MTL_INLINE static T Min(const T& a, const T& b, const T& c)
{
  return Min(Min(a,b),c);
}
template <class T> MTL_INLINE static T Min(const T& a, const T& b, const T& c, const T& d)
{
  return Min(Min(a,b),Min(c,d));
}

//
// Recursive helper classes.
//
template <int N, class T> class PowHelper
{
public:  MTL_INLINE static T Compute(const T& x)  { return PowHelper<N-1,T>::Compute(x) * x; }
};
template <class T> class PowHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T& x)  { return x; }
};

// Fixed integer powers.
template <int N, class T> MTL_INLINE static T Pow(const T& x)
{
  return PowHelper<N,T>::Compute(x);
}
template <class T> MTL_INLINE static T Square(const T& x)
{
  return Pow<2>(x);
}
template <class T> MTL_INLINE static T Cube(const T& x)
{
  return Pow<3>(x);
}

template <int N, class T> class Array
{
public:
  MTL_INLINE static T Sum(const T* p)  { return Array<N-1,T>::Sum(p) + p[N-1];          }
  MTL_INLINE static T Min(const T* p)  { return MTL::Min(Array<N-1,T>::Min(p), p[N-1]); }
  MTL_INLINE static T Max(const T* p)  { return MTL::Max(Array<N-1,T>::Max(p), p[N-1]); }

  MTL_INLINE static T MinOfAbsolutes(const T* p)
  {
    return MTL::Min(Array<N-1,T>::MinOfAbsolutes(p), Abs(p[N-1]));
  }
  MTL_INLINE static T MaxOfAbsolutes(const T* p)
  {
    return MTL::Max(Array<N-1,T>::MaxOfAbsolutes(p), Abs(p[N-1]));
  }

  MTL_INLINE static T SumOfSquares(const T* p)
  {
    return Array<N-1,T>::SumOfSquares(p) + Square(p[N-1]);
  }
  MTL_INLINE static T Dot(const T* p1, const T* p2)
  {
    return Array<N-1,T>::Dot(p1, p2) + p1[N-1] * p2[N-1];
  }

  MTL_INLINE static void UnaryMinus(T* a, const T* b)
  {
    Array<N-1,T>::UnaryMinus(a, b);
    a[N-1] = -b[N-1];
  }
};
template <class T> class Array<1,T>
{
public:
  MTL_INLINE static T Sum(const T* p)                   { return p[0];          }
  MTL_INLINE static T Min(const T* p)                   { return p[0];          }
  MTL_INLINE static T Max(const T* p)                   { return p[0];          }
  MTL_INLINE static T MinOfAbsolutes(const T* p)        { return Abs(p[0]);     }
  MTL_INLINE static T MaxOfAbsolutes(const T* p)        { return Abs(p[0]);     }
  MTL_INLINE static T SumOfSquares(const T* p)          { return Square(p[0]);  }
  MTL_INLINE static T Dot(const T* p1, const T* p2)     { return p1[0] * p2[0]; }
  MTL_INLINE static void UnaryMinus(T* a, const T* b)   { a[0] = -b[0];         }
};

// Sum & mean of fixed sized array.
template <int N, class T> MTL_INLINE static T Sum(const T* p)
{
  return Array<N,T>::Sum(p);
}
template <int N, class T> MTL_INLINE static T Mean(const T* p)
{
  return Sum<N>(p) / T(N);
}

// Minimum and maximum.
template <int N, class T> MTL_INLINE static T Minimum(const T* p)
{
  return Array<N,T>::Min(p);
}
template <int N, class T> MTL_INLINE static T Maximum(const T* p)
{
  return Array<N,T>::Max(p);
}
template <int N, class T> MTL_INLINE static T MinimumOfAbsolutes(const T* p)
{
  return Array<N,T>::MinOfAbsolutes(p);
}
template <int N, class T> MTL_INLINE static T MaximumOfAbsolutes(const T* p)
{
  return Array<N,T>::MaxOfAbsolutes(p);
}

template <int N, class T> MTL_INLINE static T SumOfSquares(const T* p)
{
  return Array<N,T>::SumOfSquares(p);
}

// Dot product.
template <int N, class T> MTL_INLINE static T Dot(const T* p1, const T* p2)
{
  return Array<N,T>::Dot(p1, p2);
}

// Simple swap.
template <class T> MTL_INLINE static void Swap(T& a, T& b)
{
  T temp = a;
  a = b;
  b = temp;
}

// A more numerically stable version of sqrt(a*a + b*b).
template <class T>
MTL_INLINE static T Hypotenuse(const T& a, const T& b)
{
  T aa, bb;
  aa = MTL::Abs(a);
  bb = MTL::Abs(b);
  if (aa > bb)
    return aa * Sqrt(T(1.0) + Square(bb/aa));
  else
    return (bb == T(0.0) ? T(0.0) : bb * Sqrt(T(1.0) + Square(aa/bb)));
}

template <class T> MTL_INLINE static T Epsilon();
template <> MTL_INLINE static F32 Epsilon<F32>()  { return FLT_EPSILON; }
template <> MTL_INLINE static F64 Epsilon<F64>()  { return DBL_EPSILON; }

template <class T> MTL_INLINE static T EpsilonSquared()
{
  return Square(Epsilon<T>());
}

// Returns a * b + c
template <class T>
MTL_INLINE static T MultiplyAndAdd(const T& a, const T& b, const T& c)
{
  return a * b + c;
}

}  // namespace MTL

#endif  // MTL_MATH_H
