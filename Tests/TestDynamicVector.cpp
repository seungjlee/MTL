//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee
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
#include <MTL/DynamicVectorOperators.h>
#include <vector>

using namespace MTL;

static const double kTol = 1e-14;

TEST(TestClassNoDefaultConstructor)
{
  class X
  {
  public:
    X(int i) : test_(i) {}
    int test_;
  };

  DynamicVector<X> v(5, X(7));

  FOR_EACH_INDEX(v)
    MTL_EQUAL(v[vIndex].test_, 7);
}

TEST(TestSetAll)
{
  DynamicVector<double> v(100);
  v.SetAll(42);

  FOR_EACH_INDEX(v)
    MTL_EQUAL_FLOAT(v[vIndex], 42, kTol);
}

TEST(TestReductionsSmallVectors)
{
  DynamicVector<double> a1(1, 3);
  DynamicVector<double> a2(1, -7);

  MTL_EQUAL_FLOAT(Sum(a1), 3, kTol);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 3, kTol);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 9, kTol);
  MTL_EQUAL_FLOAT(Min(a1), 3, kTol);
  MTL_EQUAL_FLOAT(Max(a1), 3, kTol);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 3, kTol);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 3, kTol);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 3, kTol);

  MTL_EQUAL_FLOAT(Sum(a2), -7, kTol);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 7, kTol);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 49, kTol);
  MTL_EQUAL_FLOAT(Min(a2), -7, kTol);
  MTL_EQUAL_FLOAT(Max(a2), -7, kTol);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 7, kTol);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 7, kTol);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 7, kTol);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -21, kTol);

  double aa1[] = { 3.0, 1.2,  4.5};
  double aa2[] = {-2.0, 5.0, -1.0};
  a1.Assign(aa1, 3);
  a2.Assign(aa2, 3);

  MTL_EQUAL_FLOAT(Sum(a1), 8.7, kTol);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a1), 8.7, kTol);
  MTL_EQUAL_FLOAT(SumOfSquares(a1), 30.69, kTol);
  MTL_EQUAL_FLOAT(Min(a1), 1.2, kTol);
  MTL_EQUAL_FLOAT(Max(a1), 4.5, kTol);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a1), 1.2, kTol);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a1), 4.5, kTol);
  MTL_EQUAL_FLOAT(MaxNorm(a1), 4.5, kTol);

  MTL_EQUAL_FLOAT(Sum(a2), 2, kTol);
  MTL_EQUAL_FLOAT(SumOfAbsolutes(a2), 8, kTol);
  MTL_EQUAL_FLOAT(SumOfSquares(a2), 30, kTol);
  MTL_EQUAL_FLOAT(Min(a2), -2, kTol);
  MTL_EQUAL_FLOAT(Max(a2),  5, kTol);
  MTL_EQUAL_FLOAT(MinOfAbsolutes(a2), 1, kTol);
  MTL_EQUAL_FLOAT(MaxOfAbsolutes(a2), 5, kTol);
  MTL_EQUAL_FLOAT(MaxNorm(a2), 5, kTol);

  MTL_EQUAL_FLOAT(DotProduct(a1, a2), -4.5, kTol);
}
