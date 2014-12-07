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

#include <MTL/Test.h>
#include <MTL/Polynomial.h>

using namespace MTL;

static const double kTol = 1e-12;

TEST(TesPolynomial)
{
  double c[] = { -3, 2, 5, -9 };

  ColumnVector4D cubic(c);
  MTL_EQUAL_FLOAT(Polynomial<F64>::Evaluate(0, cubic), -9, kTol);
  MTL_EQUAL_FLOAT(Polynomial<F64>::Evaluate(1, cubic), -5, kTol);

  ColumnVector3D derivative = Polynomial<F64>::ComputeDerivative(cubic);
  MTL_EQUAL_FLOAT(derivative[0], -9, kTol);
  MTL_EQUAL_FLOAT(derivative[1],  4, kTol);
  MTL_EQUAL_FLOAT(derivative[2],  5, kTol);

  ColumnVector2D secondDerivative = Polynomial<F64>::ComputeSecondDerivative(cubic);
  MTL_EQUAL_FLOAT(secondDerivative[0], -18, kTol);
  MTL_EQUAL_FLOAT(secondDerivative[1],   4, kTol);
}
