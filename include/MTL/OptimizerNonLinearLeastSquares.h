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



#ifndef MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H
#define MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H

#include "DynamicMatrix.h"

namespace MTL
{

template <I32 N, class T>
class OptimizerNonLinearLeastSquares
{
public:
  typedef ColumnVector<N> Parameters;

  OptimizerNonLinearLeastSquares(SizeType inputDataSize)
    : MaxIterations_(100), SquaredParametersDeltaTolerance_(N * Square(Epsilon<T>())),
      FiniteDifferenceDelta_(Sqrt(Epsilon<T>())),
      currentResiduals_(inputDataSize), residualsPlusDelta_(inputDataSize),
      newResiduals_(inputDataSize)
  {
    FiniteDifferenceTwoDelta_ = FiniteDifferenceDelta_ * T(2);
  }

  void MaxIterations(U32 i)             { MaxIterations_ = i;                             }
  void ParametersDeltaTolerance(T tol)  { SquaredParametersDeltaTolerance_ = Square(tol); }

  virtual void Optimize(Parameters& parameters) = 0;

protected:
  virtual void CostFunction(DynamicVector<T>& residuals, const Parameters& parameters) = 0;

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
    Parameters forwardDifferenceParameters = parameters;

    for (I32 i = 0; i < N; i++)
    {
      forwardDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta_, forwardDifferenceParameters);
      forwardDifferenceParameters[i] = parameters[i];

      residualsPlusDelta_ -= currentResiduals_;
      residualsPlusDelta_ /= FiniteDifferenceDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta_.Begin(), residualsPlusDelta_.Size());
    }
  }

  void ComputeJacobianCentralFiniteDifference(DynamicMatrix<T>& Jt,
                                              const Parameters& currentParameters)
  {
    Params centralDifferenceParameters = parameters;

    for (I32 i = 0; i < N; i++)
    {
      centralDifferenceParameters[i] = parameters[i] + FiniteDifferenceDelta_;
      CostFunction(residualsPlusDelta_, centralDifferenceParameters);
      centralDifferenceParameters[i] = parameters[i] - FiniteDifferenceDelta_;
      CostFunction(residualsMinusDelta_, centralDifferenceParameters, );
      centralDifferenceParameters[i] = parameters[i];

      residualsPlusDelta_ -= residualsMinusDelta_;
      residualsPlusDelta_ /= FiniteDifferenceTwoDelta_;

      OptimizedCopy(Jt[i], residualsPlusDelta_.Begin(), residualsPlusDelta_.Size());
    }
  }

  U32 MaxIterations_;
  T SquaredParametersDeltaTolerance_;
  T FiniteDifferenceDelta_;
  T FiniteDifferenceTwoDelta_;

  DynamicVector<T> currentResiduals_;
  DynamicVector<T> newResiduals_;
  DynamicVector<T> residualsPlusDelta_;
  DynamicVector<T> residualsMinusDelta_;
};

}  // namespace MTL

#endif  // MTL_OPTIMIZER_NON_LINEAR_LEAST_SQUARES_H
