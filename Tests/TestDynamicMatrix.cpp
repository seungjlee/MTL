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
#include <MTL/DynamicVectorOperators.h>
#include <MTL/Random.h>
#include <MTL/OptimizerLevenbergMarquardt.h>
#include <MTL/QR.h>

using namespace MTL;

static const double kTol = 1e-12;

TEST(TestMatrixMultiplication)
{
  Random random;

  Matrix<3,5> M1 = random.Matrix<3,5,F64>(-1, 1);
  Matrix<5,7> M2 = random.Matrix<5,7,F64>(-1, 1);

  Matrix<3,7> P0 = M1 * M2;

  DynamicMatrix<F64> A1(M1);
  DynamicMatrix<F64> A2(M2);

  DynamicMatrix<F64> P1 = A1 * A2;

  for (I32 row = 0; row < P1.Rows(); row++)
    for (I32 col = 0; col < P1.Cols(); col++)
      MTL_EQUAL_FLOAT(P1[row][col], P0[row][col], kTol);

  Matrix<3,3> Q0 = M1.MultiplyByTranspose();
  DynamicMatrix<F64> Q1 = A1 .MultiplyByTranspose();
}

TEST(TestHouseholderQR)
{
  enum
  {
    kSamples = 1024*1024+1
  };

  DynamicVector<F64> xs(kSamples);
  DynamicVector<F64> ys(kSamples);

  DynamicMatrix<F64> At(3, kSamples);

  for (int i = 0; i < kSamples; i++)
  {
    double x = 2 * i / double(kSamples) - 1;

    At[0][i] = x*x;
    At[1][i] = x;

    xs[i] = x;
    ys[i] = 2.5 * x*x + 0.3 * x - 1.7;
  }
  OptimizedAssignAll(At[2], 1.0, At.Cols());

  DynamicVector<F64> b = ys;

  I32 rank = SolveHouseholderQRTransposed(b, At);
  MTL_EQUAL(rank, 3);
  MTL_EQUAL_FLOAT(b[0],  2.5, kTol);
  MTL_EQUAL_FLOAT(b[1],  0.3, kTol);
  MTL_EQUAL_FLOAT(b[2], -1.7, kTol);

  class QuadraticOptimizer : public OptimizerLevenbergMarquardt<3,F64>
  {
  public:
    QuadraticOptimizer(const DynamicVector<F64>& x, const DynamicVector<F64>& y)
      : OptimizerLevenbergMarquardt(x.Size()), x_(x), y_(y)
    {
      assert(x.Size() == y.Size());
    }

  protected:
    virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& p)
    {
      residuals = x_ * (x_ * p[0] + p[1]) + p[2] - y_;
    }

  private:
    DynamicVector<F64> x_;
    DynamicVector<F64> y_;
  };

  QuadraticOptimizer optimizer(xs, ys);
  QuadraticOptimizer::Parameters coeffs;
  coeffs.Zeros();

  optimizer.Optimize(coeffs);

  printf("  Levenberg-Marquardt optimizer finished in %d iterations\n", optimizer.Iterations());
  printf("  Sum of squares of residuals = %e\n", optimizer.SumOfSquaresOfResiduals());

  MTL_EQUAL_FLOAT(coeffs[0],  2.5, kTol);
  MTL_EQUAL_FLOAT(coeffs[1],  0.3, kTol);
  MTL_EQUAL_FLOAT(coeffs[2], -1.7, kTol);
}

TEST(TestHouseholderQR_Speed)
{
  enum
  {
    M = 128*1024,
    N = 10,
    kRepeats = 10
  };

  DynamicMatrix<F64> At;
  DynamicVector<F64> b;

  Timer t;
  F64 maxRMS = 0;

  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    At = random.DynamicMatrix<F64>(N,M, -1, 1);
    b = random.DynamicVector(M, -1, 1);

    for (I32 row = 0; row < At.Rows(); row++)
      At[row][row] += 10.0;

    DynamicVector<F64> x = b;

    t.Start();
    I32 rank = SolveHouseholderQRTransposed(x, At);
    t.Stop();
    MTL_EQUAL(rank, N);

    DynamicMatrix<F64> A = At.ComputeTranspose();
    DynamicVector<F64> residuals = A*x - b;
    //printf(" RMS = %e\n", RMS(residuals));

    maxRMS = Max(maxRMS, RMS(residuals));
  }

  printf("  QR solver: %9.3f msecs (%d times), Max RMS = %e\n",
         t.Milliseconds(), kRepeats, maxRMS);
  MTL_LESS_THAN(maxRMS, 1.5);
}
