//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee
//

#ifndef MTL_MATH_H
#define MTL_MATH_H

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
template <int N, class T> class SumHelper
{
public:  MTL_INLINE static T Compute(const T* p)  { return SumHelper<N-1,T>::Compute(p) + p[N-1]; }
};
template <class T> class SumHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T* p)  { return p[0]; }
};
template <int N, class T> class PowHelper
{
public:  MTL_INLINE static T Compute(const T& x)  { return PowHelper<N-1,T>::Compute(x) * x; }
};
template <class T> class PowHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T& x)  { return x; }
};
template <int N, class T> class MinHelper
{
public:
  MTL_INLINE static T Compute(const T* p)  { return Min(MinHelper<N-1,T>::Compute(p), p[N-1]); }
};
template <class T> class MinHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T* p)  { return p[0]; }
};
template <int N, class T> class MaxHelper
{
public:
  MTL_INLINE static T Compute(const T* p)  { return Max(MaxHelper<N-1,T>::Compute(p), p[N-1]); }
};
template <class T> class MaxHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T* p)  { return p[0]; }
};
template <int N, class T> class DotHelper
{
public:  MTL_INLINE static T Compute(const T* p1, const T* p2)
         { return DotHelper<N-1,T>::Compute(p1, p2) + p1[N-1] * p2[N-1]; }
};
template <class T> class DotHelper<1,T>
{
public:  MTL_INLINE static T Compute(const T* p1, const T* p2)  { return p1[0] * p2[0]; }
};

// Sum & mean of fixed sized array.
template <int N, class T> MTL_INLINE static T Sum(const T* p)
{
  return SumHelper<N,T>::Compute(p);
}
template <int N, class T> MTL_INLINE static T Mean(const T* p)
{
  return Sum<N>(p) * (1.0/N);
}

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

// Minimum and maximum.
template <int N, class T> MTL_INLINE static T Minimum(const T* p)
{
  return MinHelper<N,T>::Compute(p);
}
template <int N, class T> MTL_INLINE static T Maximum(const T* p)
{
  return MaxHelper<N,T>::Compute(p);
}

// Dot product.
template <int N, class T> MTL_INLINE static T Dot(const T* p1, const T* p2)
{
  return DotHelper<N,T>::Compute(p1, p2);
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
    return aa * MTL::Sqrt(1.0 + Square(bb/aa));
  else
    return (bb == 0.0 ? 0.0 : bb * Sqrt(1.0 + MTL::Square<2>(aa/bb)));
}

template <class T> MTL_INLINE static T Epsilon();
template <> MTL_INLINE static F32 Epsilon<F32>()  { return FLT_EPSILON; }
template <> MTL_INLINE static F64 Epsilon<F64>()  { return DBL_EPSILON; }

}  // namespace MTL

#endif  // MTL_MATH_H
