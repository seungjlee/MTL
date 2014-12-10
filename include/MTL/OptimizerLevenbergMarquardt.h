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


#ifndef MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H
#define MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H

#include "OptimizerNonLinearLeastSquares.h"
#include "LDLt.h"

namespace MTL
{

template <I32 N, class T>
class OptimizerLevenbergMarquardt : public OptimizerNonLinearLeastSquares<N,T>
{
public:
  OptimizerLevenbergMarquardt(SizeType inputDataSize)
    : OptimizerNonLinearLeastSquares<N,T>(inputDataSize),
      LDLtRankTolerance_(N * Epsilon<T>())
  {
  }

  virtual void Optimize(Parameters& parameters)
  {
    T v = T(2);
    SquareMatrix<N,T> A;

    CostFunction(currentResiduals_, parameters);
    BestSumOfSquaresOfResiduals_ = SumOfSquares(currentResiduals_);

    DynamicMatrix<T> Jt(N, (I32)currentResiduals_.Size());

    ComputeJacobian(Jt, parameters);
    MultiplyByTranspose(A[0], Jt[0], N, Jt.Cols(), N, Jt.RowSize());

    T maxDiagonal = A[0][0];
    for (I32 i = 1; i < N; i++)
      maxDiagonal = Max(maxDiagonal, A[i][i]);

    T mu = T(1e-3);
    T damping = maxDiagonal * mu;

    DynamicVector<T> G;

    Iterations_ = 0;
    bool done = false;

    Parameters delta;

    do
    {
      Iterations_++;

      T p = 0;
      do
      {
        A.AddToDiagonals(mu);
        G = Jt * currentResiduals_;
        memcpy(&delta[0], G.Begin(), N*sizeof(T));

        I32 rank = SolveLDLt(delta, A, LDLtRankTolerance_);

        if (rank == N)
        {
          if (delta.SumOfSquares() > ParametersDeltaTolerance_)
          {
            Parameters newParameters = parameters;
            newParameters -= delta;

            CostFunction(newResiduals_, newParameters);

            double newSumOfSquaresOfResiduals = SumOfSquares(newResiduals_);
            if (newSumOfSquaresOfResiduals <= BestSumOfSquaresOfResiduals_)
            {
              BestSumOfSquaresOfResiduals_ = newSumOfSquaresOfResiduals;
              parameters = newParameters;
              currentResiduals_ = newResiduals_;

              ComputeJacobian(Jt, parameters);
              MultiplyByTranspose(A[0], Jt[0], N, Jt.Cols(), N, Jt.RowSize());

              mu = mu * Max(T(kOneThird), T(1) - Cube(T(2)*p - T(1)));

              if (Abs(mu) < Epsilon<T>())
                done = true;

              v = T(2);
            }
            else
            {
              mu *= v; 
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
          mu *= v; 
          v *= T(2);
        }
      }
      while(!done && p < 0);
    }
    while(!done && Iterations_ < MaxIterations_);
  }

  U32 Iterations() const             { return Iterations_;                  }
  T SumOfSquaresOfResiduals() const  { return BestSumOfSquaresOfResiduals_; }

private:
  T LDLtRankTolerance_;
  T BestSumOfSquaresOfResiduals_;
  U32 Iterations_;
};

}  // namespace MTL

#endif  // MTL_OPTIMIZER_LEVENBERG_MARQUARDT_H
