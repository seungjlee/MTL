//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_POLYNOMIAL_H
#define MTL_POLYNOMIAL_H

#include "Point2D.h"
#include "QR.h"

namespace MTL
{

enum PolynomialOrder
{
  // Number of coefficients for each order.
  kLinear    =  2,
  kQuadratic =  3,
  kCubic     =  4,
  kQuartic   =  5,
  kQuintic   =  6,
  kSextic    =  7,
  kSeptic    =  8,
  kOctic     =  9,
  kNonic     = 10,
  kDecic     = 11
};

template<I32 N, class T> class PolynomialHelper;

template<class T>
class Polynomial
{
public:
  // Returns polynomial coefficients.
  // Fit points so that:
  //
  //   y = c[0] * x^(N-1) + c[1] * x^(N-2) + ... + c[N-2] * x + c[N-1]
  //
  template <I32 N>
  MTL_INLINE static ColumnVector<N,T> Fit(const DynamicVector<Point2D<T>>& pts)
  {
    DynamicMatrix<T> At;
    DynamicVector<T> y;

    return Fit<N>(At, y, pts);
  }
  template <I32 N>
  MTL_INLINE static ColumnVector<N,T> Fit(DynamicMatrix<T>& At, DynamicVector<T>& y,
                                          const DynamicVector<Point2D<T>>& pts)
  {
    At.Resize(N, (I32)pts.Size());
    y.Resize(pts.Size());

    FOR_EACH_INDEX(pts)
    {
      T x = pts[ptsIndex].x();

      At[N-2][ptsIndex] = x;
      for (int k = N-3; k >=0; k--)
        At[k][ptsIndex] = At[k+1][ptsIndex] * x;

      y[ptsIndex] = pts[ptsIndex].y();
    }
    OptimizedAssignAll(At[N-1], T(1), pts.Size());

    SolveHouseholderQRTransposed(y, At);
    ColumnVector<N,T> x;
    memcpy(&x[0], y.Begin(), sizeof(x));
    return x;
  }
  template <I32 N>
  MTL_INLINE static ColumnVector<N,T> Fit(const DynamicVector<T>& xs, const DynamicVector<T>& ys)
  {
    DynamicMatrix<T> At;
    DynamicVector<T> y;

    return Fit<N>(At, y, xs, ys);
  }
  template <I32 N>
  MTL_INLINE static ColumnVector<N,T> Fit(DynamicMatrix<T>& At, DynamicVector<T>& y,
                                          const DynamicVector<T>& xs, const DynamicVector<T>& ys)
  {
    assert(xs.Size() == ys.Size());

    At.Resize(N, (I32)xs.Size());
    y = ys;

    FOR_EACH_INDEX(xs)
    {
      T x = xs[xsIndex];

      At[N-2][xsIndex] = x;
      for (int k = N-3; k >=0; k--)
        At[k][xsIndex] = At[k+1][xsIndex] * x;
    }
    OptimizedAssignAll(At[N-1], T(1), xs.Size());

    SolveHouseholderQRTransposed(y, At);
    ColumnVector<N,T> x;
    memcpy(&x[0], y.Begin(), sizeof(x));
    return x;
  }

  template <I32 N> MTL_INLINE static T Evaluate(const T& x, const T* coefficients)
  {
    assert(N >= 3);
    return PolynomialHelper<N-1,T>::Evaluate(x, coefficients);
  }

  template <I32 N>
  MTL_INLINE static T Evaluate(const T& x, const ColumnVector<N,T>& coefficients)
  {
    return Evaluate<N>(x, &coefficients[0]);
  }

  template <I32 N> MTL_INLINE static void ComputeDerivative(T* derivative, const T* coefficients)
  {
    PolynomialHelper<N,T>::ComputeDerivative(derivative, coefficients);
  }

  template <I32 N>
  MTL_INLINE static ColumnVector<N-1,T> ComputeDerivative(const ColumnVector<N,T>& coefficients)
  {
    ColumnVector<N-1,T> derivative;
    ComputeDerivative<N-1>(&derivative[0], &coefficients[0]);
    return derivative;
  }

  template <I32 N>
  MTL_INLINE static void ComputeSecondDerivative(T* derivative, const T* coefficients)
  {
    PolynomialHelper<N,T>::ComputeSecondDerivative(derivative, coefficients);
  }

  template <I32 N> MTL_INLINE static ColumnVector<N-2,T>
  ComputeSecondDerivative(const ColumnVector<N,T>& coefficients)
  {
    ColumnVector<N-2,T> derivative;
    ComputeSecondDerivative<N-2>(&derivative[0], &coefficients[0]);
    return derivative;
  }

  // Solves f(x) = 0 where f(x) is a polynomial.
  template <I32 N, I32 MAX_ITERATIONS>
  MTL_INLINE static T Solve(const ColumnVector<N,T>& coefficients, double initial = 0,
                            const T& tolerance = NumericalEpsilon<T>())
  {
    ColumnVector<N-1> derivative = Polynomial<T>::ComputeDerivative(coefficients);
    T solution = initial;
    T delta;

    for (I32 i = 0; i < MAX_ITERATIONS; i++)
    {
      delta = Polynomial<T>::Evaluate(solution, coefficients) /
              Polynomial<T>::Evaluate(solution, derivative);
      if (Abs(delta) < tolerance)
        break;

      solution -= delta;
    }

    return solution;
  }
};

template<I32 N, class T> class PolynomialHelper
{
public:
  MTL_INLINE static T Evaluate(const T& x, const T* coefficients)
  {
    return PolynomialHelper<N-1,T>::Evaluate(x, coefficients) * x + coefficients[N];
  }
  MTL_INLINE static void ComputeDerivative(T* derivative, const T* coefficients)
  {
    derivative[0] = coefficients[0] * N;
    PolynomialHelper<N-1,T>::ComputeDerivative(derivative + 1, coefficients + 1);
  }

  MTL_INLINE static void ComputeSecondDerivative(T* derivative, const T* coefficients)
  {
    derivative[0] = coefficients[0] * (N * (N + 1));
    PolynomialHelper<N-1,T>::ComputeSecondDerivative(derivative + 1, coefficients + 1);
  }
};

template<class T> class PolynomialHelper<1,T>
{
public:
  MTL_INLINE static T Evaluate(const T& x, const T* coefficients)
  {
    return coefficients[0] * x + coefficients[1];
  }
  MTL_INLINE static void ComputeDerivative(T* derivative, const T* coefficients)
  {
    derivative[0] = coefficients[0];
  }
  MTL_INLINE static void ComputeSecondDerivative(T* derivative, const T* coefficients)
  {
    derivative[0] = coefficients[0] * 2;
  }
};

}  // namespace MTL

#endif // MTL_POLYNOMIAL_H
