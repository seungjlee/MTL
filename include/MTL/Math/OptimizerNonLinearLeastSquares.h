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



#ifndef MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H
#define MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H

#include "DynamicMatrix.h"

namespace MTL
{

template <I32 N, class T>
class OptimizerNonLinearLeastSquares
{
public:
  typedef ColumnVector<N,T> Parameters;

  OptimizerNonLinearLeastSquares()
    : MaxIterations_(100), SquaredParametersDeltaTolerance_(Square(NumericalEpsilon<T>())),
      FiniteDifferenceDelta_(Sqrt(NumericalEpsilon<T>())),
      FiniteDifferenceTwoDelta_(FiniteDifferenceDelta_ * T(2))
  {
  }
  OptimizerNonLinearLeastSquares(SizeType residualsSize)
    : MaxIterations_(100), SquaredParametersDeltaTolerance_(Square(NumericalEpsilon<T>())),
      FiniteDifferenceDelta_(Sqrt(NumericalEpsilon<T>())),
      FiniteDifferenceTwoDelta_(FiniteDifferenceDelta_ * T(2))
  {
    Reset(residualsSize);
  }

  void MaxIterations(U32 i)             { MaxIterations_ = i;                             }
  void ParametersDeltaTolerance(T tol)  { SquaredParametersDeltaTolerance_ = Square(tol); }

  virtual void Optimize(Parameters& parameters) = 0;
  virtual void Reset(SizeType residualsSize)
  {
    NewResiduals_.Resize(residualsSize);
    CurrentResiduals_.Resize(residualsSize);
  }

protected:
  virtual void CostFunction(DynamicVector<T>& residuals, const Parameters& parameters) = 0;

  // User might need to override this if CostFunction is not thread safe.
  virtual void ParallelCostFunction(DynamicVector<T>& residuals, const Parameters& parameters)
  {
    CostFunction(residuals, parameters);
  }

  // User can override to use a more accurate Jacobian with partial derivatives or a different
  // finite difference method such as ComputeJacobianCentralFiniteDifference.
  virtual void ComputeJacobian(DynamicMatrix<T>& Jt,
                               const Parameters& currentParameters)
  {
    ComputeJacobianForwardFiniteDifference(Jt, currentParameters);
  }

  void ComputeJacobianForwardFiniteDifference(DynamicMatrix<T>& Jt,
                                              const Parameters& parameters)
  {
    for (I32 i = 0; i < N; i++)
    {
      Parameters forwardDifferenceParameters = parameters;
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      forwardDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta, forwardDifferenceParameters);
      forwardDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= CurrentResiduals_;
      residualsPlusDelta /= FiniteDifferenceDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void ComputeJacobianCentralFiniteDifference(DynamicMatrix<T>& Jt,
                                              const Parameters& parameters)
  {
    for (I32 i = 0; i < N; i++)
    {
      Parameters centralDifferenceParameters = parameters;
      DynamicVector<T> residualsMinusDelta(CurrentResiduals_.Size());
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      centralDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i] - FiniteDifferenceDelta_;
      CostFunction(residualsMinusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= residualsMinusDelta;
      residualsPlusDelta /= FiniteDifferenceTwoDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void ParallelComputeJacobianForwardFiniteDifference(DynamicMatrix<T>& Jt,
                                                      const Parameters& parameters)
  {
    MTL_PARALLEL_FOR_BLOCKS(N)
    for (I32 i = 0; i < N; i++)
    {
      Parameters forwardDifferenceParameters = parameters;
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      forwardDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      ParallelCostFunction(residualsPlusDelta, forwardDifferenceParameters);
      forwardDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= CurrentResiduals_;
      residualsPlusDelta /= FiniteDifferenceDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void ParallelComputeJacobianCentralFiniteDifference(DynamicMatrix<T>& Jt,
                                                      const Parameters& parameters)
  {
    MTL_PARALLEL_FOR_BLOCKS(N)
    for (I32 i = 0; i < N; i++)
    {
      Parameters centralDifferenceParameters = parameters;
      DynamicVector<T> residualsMinusDelta(CurrentResiduals_.Size());
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      centralDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      ParallelCostFunction(residualsPlusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i] - FiniteDifferenceDelta_;
      ParallelCostFunction(residualsMinusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= residualsMinusDelta;
      residualsPlusDelta /= FiniteDifferenceTwoDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  U32 MaxIterations_;
  T SquaredParametersDeltaTolerance_;
  T FiniteDifferenceDelta_;
  T FiniteDifferenceTwoDelta_;

  DynamicVector<T> CurrentResiduals_;
  DynamicVector<T> NewResiduals_;
};

// Dynamic version.
template <class T>
class DynamicOptimizerNonLinearLeastSquares
{
public:
  typedef DynamicVector<T> Parameters;

  DynamicOptimizerNonLinearLeastSquares()
    : MaxIterations_(100), SquaredParametersDeltaTolerance_(Square(NumericalEpsilon<T>())),
      FiniteDifferenceDelta_(Sqrt(NumericalEpsilon<T>())),
      FiniteDifferenceTwoDelta_(FiniteDifferenceDelta_ * T(2))
  {
  }
  DynamicOptimizerNonLinearLeastSquares(SizeType residualsSize)
    : MaxIterations_(100), SquaredParametersDeltaTolerance_(Square(NumericalEpsilon<T>())),
      FiniteDifferenceDelta_(Sqrt(NumericalEpsilon<T>())),
      FiniteDifferenceTwoDelta_(FiniteDifferenceDelta_ * T(2))
  {
    Reset(residualsSize);
  }

  void MaxIterations(U32 i)             { MaxIterations_ = i;                             }
  void ParametersDeltaTolerance(T tol)  { SquaredParametersDeltaTolerance_ = Square(tol); }

  virtual void Optimize(Parameters& parameters) = 0;

  virtual void Reset(SizeType residualsSize)
  {
    NewResiduals_.Resize(residualsSize);
    CurrentResiduals_.Resize(residualsSize);
  }

protected:
  virtual void CostFunction(DynamicVector<T>& residuals, const Parameters& parameters) = 0;

  // User might need to override this if CostFunction is not thread safe.
  virtual void ParallelCostFunction(DynamicVector<T>& residuals, const Parameters& parameters)
  {
    CostFunction(residuals, parameters);
  }

  // User can override to use a more accurate Jacobian with partial derivatives or a different
  // finite difference method such as ComputeJacobianCentralFiniteDifference.
  virtual void ComputeJacobian(DynamicMatrix<T>& Jt,
                               const Parameters& currentParameters)
  {
    ComputeJacobianForwardFiniteDifference(Jt, currentParameters);
  }

  void ComputeJacobianForwardFiniteDifference(DynamicMatrix<T>& Jt,
                                              const Parameters& parameters)
  {
    for (I32 i = 0; i < (I32)parameters.Size(); i++)
    {
      Parameters forwardDifferenceParameters = parameters;
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      forwardDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta, forwardDifferenceParameters);
      forwardDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= CurrentResiduals_;
      residualsPlusDelta /= FiniteDifferenceDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void ComputeJacobianCentralFiniteDifference(DynamicMatrix<T>& Jt,
                                              const Parameters& parameters)
  {
    for (I32 i = 0; i < (I32)parameters.Size(); i++)
    {
      Parameters centralDifferenceParameters = parameters;
      DynamicVector<T> residualsMinusDelta(CurrentResiduals_.Size());
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      centralDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i] - FiniteDifferenceDelta_;
      CostFunction(residualsMinusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= residualsMinusDelta;
      residualsPlusDelta /= FiniteDifferenceTwoDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void ParallelComputeJacobianForwardFiniteDifference(DynamicMatrix<T>& Jt,
                                                      const Parameters& parameters)
  {
    MTL_PARALLEL_FOR_BLOCKS(parameters.Size())
    for (I32 i = 0; i < (I32)parameters.Size(); i++)
    {
      Parameters forwardDifferenceParameters = parameters;
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      forwardDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      ParallelCostFunction(residualsPlusDelta, forwardDifferenceParameters);
      forwardDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= CurrentResiduals_;
      residualsPlusDelta /= FiniteDifferenceDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  void  ParallelComputeJacobianCentralFiniteDifference(DynamicMatrix<T>& Jt,
                                                       const Parameters& parameters)
  {
    MTL_PARALLEL_FOR_BLOCKS(parameters.Size())
    for (I32 i = 0; i < (I32)parameters.Size(); i++)
    {
      Parameters centralDifferenceParameters = parameters;
      DynamicVector<T> residualsMinusDelta(CurrentResiduals_.Size());
      DynamicVector<T> residualsPlusDelta(CurrentResiduals_.Size());

      centralDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      ParallelCostFunction(residualsPlusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i] - FiniteDifferenceDelta_;
      ParallelCostFunction(residualsMinusDelta, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i];

      residualsPlusDelta -= residualsMinusDelta;
      residualsPlusDelta /= FiniteDifferenceTwoDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta.Begin(), residualsPlusDelta.Size());
    }
  }

  U32 MaxIterations_;
  T SquaredParametersDeltaTolerance_;
  T FiniteDifferenceDelta_;
  T FiniteDifferenceTwoDelta_;

  DynamicVector<T> CurrentResiduals_;
  DynamicVector<T> NewResiduals_;
};

}  // namespace MTL

#endif  // MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H
