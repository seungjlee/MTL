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


#ifndef MTL_POLYNOMIAL_H
#define MTL_POLYNOMIAL_H

#include "Matrix.h"

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

template<class T>
class Polynomial
{
public:
  template <I32 N> MTL_INLINE static T Evaluate(const T& x, const T* coefficients)
  {
    assert(N >= 3);
    return Evaluate<N-1>(x, coefficients) * x + coefficients[N-1];
  }
  template <> MTL_INLINE static T Evaluate<2>(const T& x, const T* coefficients)
  {
    return coefficients[0] * x + coefficients[1];
  }

  template <I32 N>
  MTL_INLINE static T Evaluate(const T& x, const ColumnVector<N,T>& coefficients)
  {
    return Evaluate<N>(x, &coefficients[0]);
  }

  template <I32 N> MTL_INLINE static void ComputeDerivative(T* derivative,
                                                                 const T* coefficients)
  {
    derivative[0] = coefficients[0] * N;
    ComputeDerivative<N-1>(derivative + 1, coefficients + 1);
  }
  template <> MTL_INLINE static void ComputeDerivative<1>(T* derivative,
                                                               const T* coefficients)
  {
    derivative[0] = coefficients[0];
  }

  template <I32 N> MTL_INLINE static ColumnVector<N-1,T>
  ComputeDerivative(const ColumnVector<N,T>& coefficients)
  {
    ColumnVector<N-1,T> derivative;
    ComputeDerivative<N-1>(&derivative[0], &coefficients[0]);
    return derivative;
  }

  template <I32 N> MTL_INLINE static void ComputeSecondDerivative(T* derivative,
                                                                  const T* coefficients)
  {
    derivative[0] = coefficients[0] * (N * (N + 1));
    ComputeSecondDerivative<N-1>(derivative + 1, coefficients + 1);
  }
  template <> MTL_INLINE static void ComputeSecondDerivative<1>(T* derivative,
                                                                const T* coefficients)
  {
    derivative[0] = coefficients[0] * 2;
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
                            const T& tolerance = M * Epsilon<T>())
  {
    ColumnVector<N-1> derivative = FastPolynomial<T>::ComputeDerivative(coefficients);
    T solution = initial;
    T delta;

    for (long i = 0; i < MAX_ITERATIONS; i++)
    {
      delta = FastPolynomial<T>::evaluate(solution, coefficients) /
              FastPolynomial<T>::evaluate(solution, derivative);
      if (Abs(delta) < tolerance)
        break;

      solution -= delta;
    }

    return solution;
  }
};

}  // namespace MTL

#endif // MTL_POLYNOMIAL_H