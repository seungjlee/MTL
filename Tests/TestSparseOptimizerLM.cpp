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
#include <MTL/Math/SparseOptimizerLevenbergMarquardt.h>

using namespace MTL;

// Fit the simple linear model y = a + b_{seg(i)} where parameter 0 (a) is
// shared across all residuals and parameters 1..K (b_k) are per-segment. The
// resulting Jacobian has an "arrow" sparsity pattern: every row touches column
// 0 plus exactly one segment column. J^T*J is then arrow-shaped, which gives
// the COLAMD ordering and SymbolicLDLt analysis non-trivial work.
//
// This concrete subclass exercises:
//   - SparseOptimizerLevenbergMarquardt::Optimize (LM loop on sparse normal eq.)
//   - CompressedSparseMatrix::Create (sparse column construction)
//   - CompressedSparseMatrix::MultiplyTransposeByThisParallel (A = J^T J)
//   - SymbolicLDLt analyze + SolveLDLt with the resulting structure
//   - DynamicOptimizerNonLinearLeastSquares::Reset and CostFunction plumbing
class ArrowFit : public SparseOptimizerLevenbergMarquardt<F64>
{
public:
  ArrowFit(I32 segments, I32 pointsPerSegment, F64 trueA, const DynamicVector<F64>& trueBs)
    : Segments_(segments), PointsPerSegment_(pointsPerSegment)
  {
    SizeType totalPoints = SizeType(segments) * SizeType(pointsPerSegment);
    Ys_.Resize(totalPoints);

    for (I32 s = 0; s < segments; s++)
      for (I32 i = 0; i < pointsPerSegment; i++)
        Ys_[s * pointsPerSegment + i] = trueA + trueBs[s];

    Reset(totalPoints);
  }

protected:
  virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& parameters) override
  {
    F64 a = parameters[0];
    for (I32 s = 0; s < Segments_; s++)
    {
      F64 b = parameters[1 + s];
      for (I32 i = 0; i < PointsPerSegment_; i++)
      {
        I32 r = s * PointsPerSegment_ + i;
        residuals[r] = (a + b) - Ys_[r];
      }
    }
  }

  virtual void ComputeSparsityMatrix(I32 numberOfParams) override
  {
    I32 totalPoints = Segments_ * PointsPerSegment_;
    DynamicVector<DynamicVector<I32>> sparseColumns(numberOfParams);

    // Column 0 (parameter a) touches every residual.
    sparseColumns[0].Reserve(totalPoints);
    for (I32 r = 0; r < totalPoints; r++)
      sparseColumns[0].PushBack(r);

    // Column s+1 (parameter b_s) touches only residuals in segment s.
    for (I32 s = 0; s < Segments_; s++)
    {
      sparseColumns[1 + s].Reserve(PointsPerSegment_);
      for (I32 i = 0; i < PointsPerSegment_; i++)
        sparseColumns[1 + s].PushBack(s * PointsPerSegment_ + i);
    }

    J_.Create(totalPoints, numberOfParams, sparseColumns,
              /*updateOptimizedMultiplyStructure=*/true);
  }

  virtual void ComputeJacobian(CompressedSparseMatrix<F64>& J, const Parameters& /*parameters*/) override
  {
    // Linear model: dr/da = 1, dr/db_seg = 1. The column structure was set up
    // in ComputeSparsityMatrix; we only need to fill values.
    F64* Ax = const_cast<F64*>(J.Ax());
    SizeType nnz = J.NumberOfElements();
    for (SizeType k = 0; k < nnz; k++)
      Ax[k] = 1.0;
  }

private:
  I32 Segments_;
  I32 PointsPerSegment_;
  DynamicVector<F64> Ys_;
};

TEST(Test_SparseLM_ArrowFit)
{
  const I32 K = 8;
  const I32 P = 12;
  const F64 trueA = 0.7;
  DynamicVector<F64> trueBs(K);
  for (I32 s = 0; s < K; s++)
    trueBs[s] = 0.3 * F64(s) - 1.1;

  ArrowFit fitter(K, P, trueA, trueBs);
  fitter.MaxIterations(50);

  // The model y = a + b_s is rank-deficient by one (only a + b_s is observed,
  // not a or b_s separately). LM with LDLt rank check will reject the step,
  // which exercises the rank-deficient branch of the optimizer. We still want
  // the run to complete cleanly.
  DynamicVector<F64> p(K + 1);
  for (SizeType i = 0; i < p.Size(); i++)
    p[i] = 0.0;

  fitter.Optimize(p);

  // Verify the fitted sums match the data even though individual parameters
  // are unidentified.
  F64 a = p[0];
  for (I32 s = 0; s < K; s++)
    MTL_EQUAL_FLOAT(a + p[1 + s], trueA + trueBs[s], 1e-9);

  Out() << L"SparseLM iters=" << fitter.Iterations()
        << L"  ssr=" << fitter.SumOfSquaresOfResiduals() << std::endl;
}

// Fit y = a + b_s + c_s * x with K segments and P points each. Adds a third
// parameter per segment so J^T*J has 1+2K parameters with arrow-plus-block
// structure, which gives COLAMD a more interesting elimination problem.
class ArrowLineFit : public SparseOptimizerLevenbergMarquardt<F64>
{
public:
  ArrowLineFit(I32 segments, I32 pointsPerSegment, F64 trueA,
               const DynamicVector<F64>& trueBs,
               const DynamicVector<F64>& trueCs)
    : Segments_(segments), PointsPerSegment_(pointsPerSegment)
  {
    SizeType totalPoints = SizeType(segments) * SizeType(pointsPerSegment);
    Xs_.Resize(totalPoints);
    Ys_.Resize(totalPoints);

    for (I32 s = 0; s < segments; s++)
      for (I32 i = 0; i < pointsPerSegment; i++)
      {
        I32 r = s * pointsPerSegment + i;
        F64 x = -1.0 + 2.0 * F64(i) / F64(pointsPerSegment - 1);
        Xs_[r] = x;
        Ys_[r] = trueA + trueBs[s] + trueCs[s] * x;
      }

    Reset(totalPoints);
  }

protected:
  virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& parameters) override
  {
    F64 a = parameters[0];
    for (I32 s = 0; s < Segments_; s++)
    {
      F64 b = parameters[1 + 2*s];
      F64 c = parameters[2 + 2*s];
      for (I32 i = 0; i < PointsPerSegment_; i++)
      {
        I32 r = s * PointsPerSegment_ + i;
        residuals[r] = (a + b + c * Xs_[r]) - Ys_[r];
      }
    }
  }

  virtual void ComputeSparsityMatrix(I32 numberOfParams) override
  {
    I32 totalPoints = Segments_ * PointsPerSegment_;
    DynamicVector<DynamicVector<I32>> sparseColumns(numberOfParams);

    sparseColumns[0].Reserve(totalPoints);
    for (I32 r = 0; r < totalPoints; r++)
      sparseColumns[0].PushBack(r);

    for (I32 s = 0; s < Segments_; s++)
    {
      sparseColumns[1 + 2*s].Reserve(PointsPerSegment_);
      sparseColumns[2 + 2*s].Reserve(PointsPerSegment_);
      for (I32 i = 0; i < PointsPerSegment_; i++)
      {
        I32 r = s * PointsPerSegment_ + i;
        sparseColumns[1 + 2*s].PushBack(r);
        sparseColumns[2 + 2*s].PushBack(r);
      }
    }

    J_.Create(totalPoints, numberOfParams, sparseColumns, true);
  }

  virtual void ComputeJacobian(CompressedSparseMatrix<F64>& J, const Parameters& /*parameters*/) override
  {
    // Column 0 (a): all 1.0
    // Column 1+2s (b_s): all 1.0
    // Column 2+2s (c_s): values are Xs_[r] for the segment's residuals.
    F64* Ax = const_cast<F64*>(J.Ax());
    const I32* Ap = J.Ap();
    const I32* Ai = J.Ai();
    I32 N = J.Cols();

    for (I32 col = 0; col < N; col++)
    {
      bool isCColumn = (col >= 1) && (((col - 1) % 2) == 1);
      for (I32 k = Ap[col]; k < Ap[col + 1]; k++)
      {
        I32 row = Ai[k];
        Ax[k] = isCColumn ? Xs_[row] : 1.0;
      }
    }
  }

private:
  I32 Segments_;
  I32 PointsPerSegment_;
  DynamicVector<F64> Xs_;
  DynamicVector<F64> Ys_;
};

TEST(Test_SparseLM_ArrowLineFit)
{
  const I32 K = 6;
  const I32 P = 10;
  const F64 trueA = 0.4;

  DynamicVector<F64> trueBs(K), trueCs(K);
  for (I32 s = 0; s < K; s++)
  {
    trueBs[s] = 0.5 * F64(s) - 1.0;
    trueCs[s] = -0.2 * F64(s) + 0.7;
  }

  ArrowLineFit fitter(K, P, trueA, trueBs, trueCs);
  fitter.MaxIterations(100);

  DynamicVector<F64> p(1 + 2 * K);
  for (SizeType i = 0; i < p.Size(); i++)
    p[i] = 0.0;

  fitter.Optimize(p);

  F64 a = p[0];
  for (I32 s = 0; s < K; s++)
  {
    F64 b = p[1 + 2*s];
    F64 c = p[2 + 2*s];
    // The slope c_s is identifiable. The intercept a + b_s is the only
    // identifiable combination of (a, b_s).
    MTL_EQUAL_FLOAT(c,         trueCs[s],            1e-9);
    MTL_EQUAL_FLOAT(a + b,     trueA + trueBs[s],    1e-9);
  }

  Out() << L"SparseLM iters=" << fitter.Iterations()
        << L"  ssr=" << fitter.SumOfSquaresOfResiduals() << std::endl;
}
