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

#include <MTL/Tools/Test.h>
#include <MTL/Math/OptimizerLevenbergMarquardt.h>
#include <MTL/Math/Polynomial.h>

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

TEST(TestPolynomialFit)
{
  enum
  {
    kSamples = 1024*1024
  };

  Timer t;

  DynamicVector<Point2D<F64>> pts(kSamples);
  DynamicVector<F64> xs(kSamples);
  DynamicVector<F64> ys(kSamples);

  DynamicMatrix<F64> At(4, kSamples);

  for (int i = 0; i < kSamples; i++)
  {
    double x = 2 * i / double(kSamples) - 1;

    At[0][i] = x*x;
    At[1][i] = x;

    xs[i] = x;
    ys[i] = 4.5 * x*x*x - 3.9 * x*x + 2.3 * x - 0.7;

    pts[i] = Point2D<F64>(xs[i], ys[i]);
  }

  class CubicOptimizer : public OptimizerLevenbergMarquardt<4,F64>
  {
  public:
    CubicOptimizer(const DynamicVector<F64>& x, const DynamicVector<F64>& y)
      : OptimizerLevenbergMarquardt(x.Size()), x_(x), y_(y)
    {
      assert(x.Size() == y.Size());
    }

  protected:
    virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& p)
    {
      residuals = x_ * (x_ * (x_ * p[0] + p[1]) + p[2]) + p[3] - y_;
    }

    virtual void ComputeJacobian(DynamicMatrix<F64>& Jt,
                                 const Parameters& currentParameters)
    {
      ParallelComputeJacobianForwardFiniteDifference(Jt, currentParameters);
    }

  private:
    DynamicVector<F64> x_;
    DynamicVector<F64> y_;
  };

  CubicOptimizer optimizer(xs, ys);
  CubicOptimizer::Parameters coeffs;
  coeffs.Zeros();

  t.ResetAndStart();
  optimizer.Optimize(coeffs);
  t.Stop();
  wprintf(L"  Levenberg-Marquardt optimizer finished in %d iterations, %.3f msecs\n",
          (int)optimizer.Iterations(), t.Milliseconds());
  wprintf(L"  Sum of squares of residuals = %e\n", optimizer.SumOfSquaresOfResiduals());

  t.ResetAndStart();
  ColumnVector<4,F64> polyCoeffs1 = Polynomial<F64>::Fit<4>(pts);
  t.Stop();
  wprintf(L"  Polynomial fit time 1: %.3f msecs\n", t.Milliseconds());

  FOR_EACH_INDEX(coeffs)
  {
    MTL_EQUAL_FLOAT(polyCoeffs1[(I32)coeffsIndex], coeffs[(I32)coeffsIndex], kTol);
  }

  t.ResetAndStart();
  ColumnVector<4,F64> polyCoeffs2 = Polynomial<F64>::Fit<4>(xs, ys);
  t.Stop();
  wprintf(L"  Polynomial fit time 2: %.3f msecs\n", t.Milliseconds());

  FOR_EACH_INDEX(coeffs)
  {
    MTL_EQUAL_FLOAT(polyCoeffs2[(I32)coeffsIndex], coeffs[(I32)coeffsIndex], kTol);
  }
}
