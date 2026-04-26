//
// Math Template Library
//
// Copyright (c) 2026: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/Math/Random.h>
#include <MTL/Math/SphereFitLevenbergMarquardt.h>

using namespace MTL;

// SphereFitLevenbergMarquardt is the canonical user of the fixed-size LM
// optimizer; recovering known sphere parameters from sampled points exercises
// CostFunction, the forward-finite-difference Jacobian, the LDLt solve, and
// the LM damping update loop in OptimizerLevenbergMarquardt<4,T>.
TEST(Test_SphereFit_FixedLM)
{
  Point3D<F64> trueCenter(2.0, -1.0, 0.5);
  F64 trueRadius = 3.25;

  Random rng(12345);
  DynamicVector<Point3D<F64>> points;
  points.Reserve(200);
  for (int i = 0; i < 200; i++)
  {
    // Sample roughly uniformly on the sphere using two independent normals.
    F64 u = rng.GetNext(-1.0, 1.0);
    F64 v = rng.GetNext(0.0, 6.283185307179586);
    F64 s = Sqrt(1.0 - u*u);
    Point3D<F64> p(trueCenter[0] + trueRadius * s * std::cos(v),
                   trueCenter[1] + trueRadius * s * std::sin(v),
                   trueCenter[2] + trueRadius * u);
    points.PushBack(p);
  }

  SphereFitLevenbergMarquardt<F64> fitter(points);
  fitter.MaxIterations(50);

  Point3D<F64> center(0.0, 0.0, 0.0);
  F64 radius = 1.0;

  fitter.Optimize(center, radius);

  Out() << L"Iterations: " << fitter.Iterations() << std::endl;
  Out() << L"Final SSR:  " << fitter.SumOfSquaresOfResiduals() << std::endl;
  Out() << L"Center: " << center[0] << L", " << center[1] << L", " << center[2] << std::endl;
  Out() << L"Radius: " << radius << std::endl;

  MTL_EQUAL_FLOAT(center[0], trueCenter[0], 1e-9);
  MTL_EQUAL_FLOAT(center[1], trueCenter[1], 1e-9);
  MTL_EQUAL_FLOAT(center[2], trueCenter[2], 1e-9);
  MTL_EQUAL_FLOAT(radius,    trueRadius,    1e-9);
}

// Fit y = a * exp(b * x) using a dynamic-size LM optimizer. Exercises
// DynamicOptimizerLevenbergMarquardt<T> and DynamicOptimizerNonLinearLeastSquares<T>
// including their forward-finite-difference Jacobian, normal matrix builder,
// and dynamic LDLt solve.
class ExpFit : public DynamicOptimizerLevenbergMarquardt<F64>
{
public:
  ExpFit(const DynamicVector<F64>& xs, const DynamicVector<F64>& ys) : Xs_(xs), Ys_(ys)
  {
    Reset(xs.Size());
  }

protected:
  virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& parameters) override
  {
    F64 a = parameters[0];
    F64 b = parameters[1];
    for (SizeType i = 0; i < Xs_.Size(); i++)
      residuals[i] = a * std::exp(b * Xs_[i]) - Ys_[i];
  }

private:
  DynamicVector<F64> Xs_;
  DynamicVector<F64> Ys_;
};

TEST(Test_DynamicLM_ExpFit)
{
  F64 trueA = 2.5;
  F64 trueB = -0.4;

  DynamicVector<F64> xs(50);
  DynamicVector<F64> ys(50);
  for (SizeType i = 0; i < 50; i++)
  {
    F64 x = -1.0 + 0.05 * F64(i);
    xs[i] = x;
    ys[i] = trueA * std::exp(trueB * x);
  }

  ExpFit fit(xs, ys);
  fit.MaxIterations(100);
  fit.ParametersDeltaTolerance(1e-14);

  DynamicVector<F64> p(2);
  p[0] = 1.0;
  p[1] = 0.0;
  fit.Optimize(p);

  Out() << L"a=" << p[0] << L"  b=" << p[1]
        << L"  iters=" << fit.Iterations()
        << L"  ssr=" << fit.SumOfSquaresOfResiduals() << std::endl;

  MTL_EQUAL_FLOAT(p[0], trueA, 1e-7);
  MTL_EQUAL_FLOAT(p[1], trueB, 1e-7);
}

// Exercise the central-difference and parallel Jacobian helpers to lift their
// coverage. The optimizer below uses central differences instead of forward.
class ParabolaFitCentral : public DynamicOptimizerLevenbergMarquardt<F64>
{
public:
  ParabolaFitCentral(const DynamicVector<F64>& xs, const DynamicVector<F64>& ys)
    : Xs_(xs), Ys_(ys)
  {
    Reset(xs.Size());
  }

protected:
  virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& parameters) override
  {
    F64 a = parameters[0];
    F64 b = parameters[1];
    F64 c = parameters[2];
    for (SizeType i = 0; i < Xs_.Size(); i++)
      residuals[i] = a * Xs_[i] * Xs_[i] + b * Xs_[i] + c - Ys_[i];
  }
  virtual void ComputeJacobian(DynamicMatrix<F64>& Jt, const Parameters& currentParameters) override
  {
    this->ComputeJacobianCentralFiniteDifference(Jt, currentParameters);
  }

private:
  DynamicVector<F64> Xs_;
  DynamicVector<F64> Ys_;
};

TEST(Test_DynamicLM_CentralDifferences)
{
  F64 trueA = 1.5, trueB = -2.0, trueC = 0.75;

  DynamicVector<F64> xs(40), ys(40);
  for (SizeType i = 0; i < 40; i++)
  {
    F64 x = -1.0 + 0.05 * F64(i);
    xs[i] = x;
    ys[i] = trueA * x * x + trueB * x + trueC;
  }

  ParabolaFitCentral fit(xs, ys);
  fit.MaxIterations(100);

  DynamicVector<F64> p(3);
  p[0] = 0.0; p[1] = 0.0; p[2] = 0.0;
  fit.Optimize(p);

  MTL_EQUAL_FLOAT(p[0], trueA, 1e-9);
  MTL_EQUAL_FLOAT(p[1], trueB, 1e-9);
  MTL_EQUAL_FLOAT(p[2], trueC, 1e-9);
}
