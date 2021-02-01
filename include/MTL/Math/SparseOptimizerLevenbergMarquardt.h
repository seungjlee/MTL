//
// Math Template Library
//
// Copyright (c) 2014-2021: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_SPARSE_OPTIMIZER_LEVENBERG_MARQUARDT_H
#define MTL_SPARSE_OPTIMIZER_LEVENBERG_MARQUARDT_H

#include <MTL/Math/OptimizerLevenbergMarquardt.h>
#include <MTL/Math/SparseMatrix.h>
#include <memory>

namespace MTL
{

template <class T>
class SparseOptimizerLevenbergMarquardt : public DynamicOptimizerNonLinearLeastSquares<T>
{
public:
  SparseOptimizerLevenbergMarquardt()
    : LDLtRankTolerance_(-1.0), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }
  SparseOptimizerLevenbergMarquardt(SizeType residualsSize)
    : DynamicOptimizerNonLinearLeastSquares<T>(residualsSize),
      LDLtRankTolerance_(-1.0), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }

  // Compute A = Jt*J.
  virtual void ComputeNormalMatrix(CompressedSparseMatrix<T>& A, const CompressedSparseMatrix<T>& J)
  {
    J.MultiplyTransposeByThisParallel(A);
  }

  // Solves A*x = b. b is input as x. Returns rank of A.
  virtual I32 Solve(DynamicVector<T>& x, CompressedSparseMatrix<T>& A, T tolerance)
  {
    if (SymbolicLDLt_.get() == NULL)
    {
      SymbolicLDLt_.reset(new SymbolicLDLt<T>(A));
    }

    Temp_ = x;
    return SolveLDLt(x, A, Temp_, *SymbolicLDLt_, tolerance);
  }

  virtual void Optimize(typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters& parameters)
  {
    SymbolicLDLt_.reset();

    if (LDLtRankTolerance_ < 0)
      LDLtRankTolerance_ = NumericalEpsilon<T>();

    I32 N = (I32)parameters.Size();
    ComputeSparsityMatrix(N);

    mu_ = 0;
    T v = T(2);
    Iterations_ = 0;

    this->CostFunction(this->CurrentResiduals_, parameters);
    BestSumOfSquaresOfResiduals_ = SumOfSquares(this->CurrentResiduals_);

    this->ComputeJacobian(J_, parameters);
    ComputeNormalMatrix(A_, J_);

    T maxDiagonal = A_.MaxDiagonal();

    mu_ = maxDiagonal * InitialDampingFactor_;

    typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters G(N);
    typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters delta(N);

    bool done = false;

    do
    {
      Iterations_++;

      T p = 0;
      do
      {
        A_.AddToDiagonalsIfEntryExists(mu_);
        J_.MultiplyTransposed(G, this->CurrentResiduals_);
        OptimizedCopy(delta.Begin(), G.Begin(), delta.Size());

        I32 rank = Solve(delta, A_, LDLtRankTolerance_);

        if (rank == (I32)parameters.Size())
        {
          if (SumOfSquares(delta) > this->SquaredParametersDeltaTolerance_)
          {
            typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters
              newParameters = parameters;
            newParameters -= delta;

            this->CostFunction(this->NewResiduals_, newParameters);

            T newSumOfSquaresOfResiduals = SumOfSquares(this->NewResiduals_);
            if (newSumOfSquaresOfResiduals < BestSumOfSquaresOfResiduals_)
            {
              parameters = newParameters;
              this->CurrentResiduals_ = this->NewResiduals_;

              this->ComputeJacobian(J_, parameters);
              ComputeNormalMatrix(A_, J_);

              p = BestSumOfSquaresOfResiduals_ - newSumOfSquaresOfResiduals;
              p /= DotProduct(delta, delta * mu_ + G);

              mu_ = mu_ * Max(T(kOneThird), T(1) - Cube(T(2)*p - T(1)));

              BestSumOfSquaresOfResiduals_ = newSumOfSquaresOfResiduals;

              if (mu_ < NumericalEpsilon<T>())
                done = true;

              v = T(2);
            }
            else
            {
              mu_ *= v; 
              v *= T(2);
            }
          }
          else
          {
            done = true;
          }
        }
        else
        {
          mu_ *= v; 
          v *= T(2);
        }

        if (mu_ > kMaxDamping)
          done = true;
      }
      while(!done && p < 0);
    }
    while(!done && Iterations_ < this->MaxIterations_);
  }

  U32 Iterations() const             { return Iterations_;                  }
  T SumOfSquaresOfResiduals() const  { return BestSumOfSquaresOfResiduals_; }

protected:
  CompressedSparseMatrix<T> J_;  // Jacobian matrix.
  CompressedSparseMatrix<T> A_;
  std::shared_ptr<SymbolicLDLt<T>> SymbolicLDLt_;
  DynamicVector<T> Temp_;
  T LDLtRankTolerance_;
  T BestSumOfSquaresOfResiduals_;
  T InitialDampingFactor_;
  T mu_;
  U32 Iterations_;

  virtual void ComputeSparsityMatrix(I32 numberOfParams) = 0;
  virtual void ComputeJacobian(CompressedSparseMatrix<T>& J,
                               const typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters&
                               currentParameters) = 0;
};

}  // namespace MTL

#endif  // MTL_SPARSE_OPTIMIZER_LEVENBERG_MARQUARDT_H
