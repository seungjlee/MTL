//
// Math Template Library
//
// Copyright (c) 2014-2015: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H
#define MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H

#include "DynamicVectorOperators.h"
#include "OptimizerNonLinearLeastSquares.h"
#include "SparseMatrix.h"
#include "LDLt.h"

namespace MTL
{

static const double kDefaultInitialDampingFactor = 1e-4;
static const double kMaxDamping = 1e35;

template <I32 N, class T>
class OptimizerLevenbergMarquardt : public OptimizerNonLinearLeastSquares<N,T>
{
public:
  OptimizerLevenbergMarquardt()
    : LDLtRankTolerance_(N * Epsilon<T>()), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }
  OptimizerLevenbergMarquardt(SizeType residualsSize)
    : OptimizerNonLinearLeastSquares<N,T>(residualsSize),
      LDLtRankTolerance_(N * Epsilon<T>()), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }

  // Compute A = Jt*J.
  virtual void ComputeNormalMatrix(SquareMatrix<N,T>& A, const DynamicMatrix<T>& Jt)
  {
    MultiplyByTranspose(A[0], Jt[0], N, Jt.Cols(), N, Jt.RowSize());
  }

  // Solves A*x = b. b is input as x. Returns rank of A.
  virtual I32 Solve(ColumnVector<N,T>& x, SquareMatrix<N,T>& A, T tolerance)
  {
    return SolveLDLt<N,T>(x, A, tolerance);
  }

  virtual void Optimize(typename OptimizerNonLinearLeastSquares<N,T>::Parameters& parameters)
  {
    T v = T(2);
    Iterations_ = 0;

    this->CostFunction(this->CurrentResiduals_, parameters);
    BestSumOfSquaresOfResiduals_ = SumOfSquares(this->CurrentResiduals_);

    Jt_.Resize(N, (I32)this->CurrentResiduals_.Size());
    Jt_.Zeros();
    A_.Zeros();

    this->ComputeJacobian(Jt_, parameters);
    ComputeNormalMatrix(A_, Jt_);

    T maxDiagonal = A_[0][0];
    for (I32 i = 1; i < N; i++)
      maxDiagonal = Max(maxDiagonal, A_[i][i]);

    mu_ = maxDiagonal * InitialDampingFactor_;

    DynamicVector<T> G;
    typename OptimizerNonLinearLeastSquares<N,T>::Parameters delta;

    bool done = false;

    do
    {
      Iterations_++;

      T p = 0;
      do
      {
        A_.AddToDiagonals(mu_);
        G = Jt_ * this->CurrentResiduals_;
        memcpy(&delta[0], G.Begin(), N*sizeof(T));

        I32 rank = Solve(delta, A_, LDLtRankTolerance_);

        if (rank == N)
        {
          if (delta.SumOfSquares() > this->SquaredParametersDeltaTolerance_)
          {
            typename OptimizerNonLinearLeastSquares<N,T>::Parameters newParameters = parameters;
            newParameters -= delta;

            this->CostFunction(this->NewResiduals_, newParameters);

            T newSumOfSquaresOfResiduals = SumOfSquares(this->NewResiduals_);
            if (newSumOfSquaresOfResiduals < BestSumOfSquaresOfResiduals_)
            {
              parameters = newParameters;
              this->CurrentResiduals_ = this->NewResiduals_;

              this->ComputeJacobian(Jt_, parameters);
              ComputeNormalMatrix(A_, Jt_);

              p = BestSumOfSquaresOfResiduals_ - newSumOfSquaresOfResiduals;
              p /= delta.Dot(delta * mu_ + typename OptimizerNonLinearLeastSquares<N,T>::Parameters(G.Begin()));

              mu_ = mu_ * Max(T(kOneThird), T(1) - Cube(T(2)*p - T(1)));

              BestSumOfSquaresOfResiduals_ = newSumOfSquaresOfResiduals;

              if (mu_ < Epsilon<T>())
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
  DynamicMatrix<T> Jt_;  // Jacobian matrix transposed.
  SquareMatrix<N,T> A_;
  T LDLtRankTolerance_;
  T BestSumOfSquaresOfResiduals_;
  T InitialDampingFactor_;
  T mu_;
  U32 Iterations_;
};

template <class T>
class DynamicOptimizerLevenbergMarquardt : public DynamicOptimizerNonLinearLeastSquares<T>
{
public:
  DynamicOptimizerLevenbergMarquardt()
    : LDLtRankTolerance_(-1.0), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }
  DynamicOptimizerLevenbergMarquardt(SizeType residualsSize)
    : DynamicOptimizerNonLinearLeastSquares<T>(residualsSize),
      LDLtRankTolerance_(-1.0), InitialDampingFactor_(kDefaultInitialDampingFactor)
  {
  }

  // Compute A = Jt*J.
  virtual void ComputeNormalMatrix(DynamicMatrix<T>& A, const DynamicMatrix<T>& Jt)
  {
    MultiplyByTranspose(A[0], Jt[0], Jt.Rows(), Jt.Cols(), A.RowSize(), Jt.RowSize());
  }

  // Solves A*x = b. b is input as x. Returns rank of A.
  virtual I32 Solve(DynamicVector<T>& x, DynamicMatrix<T>& A, T tolerance)
  {
    return SolveLDLt(x, A, tolerance);
  }

  virtual void Optimize(typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters& parameters)
  {
    if (LDLtRankTolerance_ < 0)
      LDLtRankTolerance_ = Epsilon<T>() * parameters.Size();

    I32 N = (I32)parameters.Size();

    T v = T(2);
    Iterations_ = 0;

    this->CostFunction(this->CurrentResiduals_, parameters);
    BestSumOfSquaresOfResiduals_ = SumOfSquares(this->CurrentResiduals_);

    Jt_.Resize(N, (I32)this->CurrentResiduals_.Size());
    A_.Resize(N, N);

    Jt_.Zeros();
    A_.Zeros();

    this->ComputeJacobian(Jt_, parameters);
    ComputeNormalMatrix(A_, Jt_);

    T maxDiagonal = A_[0][0];
    for (I32 i = 1; i < N; i++)
      maxDiagonal = Max(maxDiagonal, A_[i][i]);

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
        A_.AddToDiagonals(mu_);
        G = Jt_ * this->CurrentResiduals_;
        OptimizedCopy(delta.Begin(), G.Begin(), delta.Size());

        I32 rank = Solve(delta, A_, LDLtRankTolerance_);

        if (rank == parameters.Size())
        {
          if (SumOfSquares(delta) > this->SquaredParametersDeltaTolerance_)
          {
            typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters newParameters = parameters;
            newParameters -= delta;

            this->CostFunction(this->NewResiduals_, newParameters);

            T newSumOfSquaresOfResiduals = SumOfSquares(this->NewResiduals_);
            if (newSumOfSquaresOfResiduals < BestSumOfSquaresOfResiduals_)
            {
              parameters = newParameters;
              this->CurrentResiduals_ = this->NewResiduals_;

              this->ComputeJacobian(Jt_, parameters);
              ComputeNormalMatrix(A_, Jt_);

              p = BestSumOfSquaresOfResiduals_ - newSumOfSquaresOfResiduals;
              p /= DotProduct(delta, delta * mu_ + G);

              mu_ = mu_ * Max(T(kOneThird), T(1) - Cube(T(2)*p - T(1)));

              BestSumOfSquaresOfResiduals_ = newSumOfSquaresOfResiduals;

              if (mu_ < Epsilon<T>())
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
  DynamicMatrix<T> Jt_;  // Jacobian matrix transposed.
  DynamicMatrix<T> A_;
  T LDLtRankTolerance_;
  T BestSumOfSquaresOfResiduals_;
  T InitialDampingFactor_;
  T mu_;
  U32 Iterations_;
};

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
  virtual void ComputeNormalMatrix(DynamicMatrix<T>& A, const CompressedSparseMatrix<T>& J)
  {
    J.MultiplyTransposeByThisParallel(A);
  }

  // Solves A*x = b. b is input as x. Returns rank of A.
  virtual I32 Solve(DynamicVector<T>& x, DynamicMatrix<T>& A, T tolerance)
  {
    return SolveLDLt(x, A, tolerance);
  }

  virtual void Optimize(typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters& parameters)
  {
    if (LDLtRankTolerance_ < 0)
      LDLtRankTolerance_ = Epsilon<T>() * parameters.Size();

    I32 N = (I32)parameters.Size();
    ComputeSparsityMatrix(N);

    T v = T(2);
    Iterations_ = 0;

    this->CostFunction(this->CurrentResiduals_, parameters);
    BestSumOfSquaresOfResiduals_ = SumOfSquares(this->CurrentResiduals_);

    A_.Resize(N, N);

    A_.Zeros();

    this->ComputeJacobian(J_, parameters);
    ComputeNormalMatrix(A_, J_);

    T maxDiagonal = A_[0][0];
    for (I32 i = 1; i < N; i++)
      maxDiagonal = Max(maxDiagonal, A_[i][i]);

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
        A_.AddToDiagonals(mu_);
        J_.MultiplyTransposed(G, this->CurrentResiduals_);
        OptimizedCopy(delta.Begin(), G.Begin(), delta.Size());

        I32 rank = Solve(delta, A_, LDLtRankTolerance_);

        if (rank == parameters.Size())
        {
          if (SumOfSquares(delta) > this->SquaredParametersDeltaTolerance_)
          {
            typename DynamicOptimizerNonLinearLeastSquares<T>::Parameters newParameters = parameters;
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

              if (mu_ < Epsilon<T>())
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
  DynamicMatrix<T> A_;  // This should be sparse as well but for now...
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

#endif  // MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H
