//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#include <MTL/Constants.h>
#include <float.h>
#include <cmath>


namespace MTL
{

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

// Returns input value within range.
template <class T> MTL_INLINE static T Limit(const T& input, const T& minValid, const T& maxValid)
{
  return Min(Max(input, minValid), maxValid);
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

// Simple swap.
template <class T> MTL_INLINE static void Swap(T& a, T& b)
{
  T temp = a;
  a = b;
  b = temp;
}

template <class T> MTL_INLINE static T Abs(const T& a)
{
  return std::abs(a);
}

template <class T> MTL_INLINE static T Sqrt(const T& a);
 
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

template <class T> MTL_INLINE T Epsilon();
template <> MTL_INLINE F32 Epsilon<F32>()   { return FLT_EPSILON; }
template <> MTL_INLINE F64 Epsilon<F64>()   { return DBL_EPSILON; }

// A default factor to use with Epsilon() to ensure numerical convergence on some problems.
// User can override this at compile time.
#ifndef MTL_EPSILON_FACTOR
#define MTL_EPSILON_FACTOR 2
#endif
template <class T> MTL_INLINE T EpsilonFactor()
{
  return T(MTL_EPSILON_FACTOR);
}
template <class T> MTL_INLINE T NumericalEpsilon()
{
  return EpsilonFactor<T>() * Epsilon<T>();
}

template <class T> MTL_INLINE static T NumericalEpsilonSquared()
{
  return Square(NumericalEpsilon<T>());
}

template <class T> MTL_INLINE T Infinity();
template <> MTL_INLINE F32 Infinity<F32>()  { return kINF32;      }
template <> MTL_INLINE F64 Infinity<F64>()  { return kINF;        }

template <class T> MTL_INLINE bool IsFinite(const T& a)
{ return a > -Infinity<T>() && a < Infinity<T>(); }

// Returns a * b + c
template <class T>
MTL_INLINE T MultiplyAndAdd(const T& a, const T& b, const T& c)
{
  return std::fma(a, b, c);
}

}  // namespace MTL

#endif  // MTL_MATH_H
