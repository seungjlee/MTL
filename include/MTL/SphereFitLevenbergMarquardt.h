//
// Math Template Library
//
// Copyright (c) 2016: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_SPHERE_FIT_LEVENBERG_MARQUARDT_H
#define MTL_SPHERE_FIT_LEVENBERG_MARQUARDT_H

#include "OptimizerLevenbergMarquardt.h"
#include "Point3D.h"

namespace MTL
{

template<class T>
class SphereFitLevenbergMarquardt : public OptimizerLevenbergMarquardt<4,T>
{
public:
  SphereFitLevenbergMarquardt(const DynamicVector<Point3D<T>>& points)
    : Points_(points)
  {
    Reset(Points_.Size());
  }

  void Optimize(Point3D<T>& center, double& radius)
  {
    Parameters variables;
    memcpy(&variables[0], &center[0], sizeof(center));
    variables[3] = radius;

    OptimizerLevenbergMarquardt<4,T>::Optimize(variables);

    memcpy(&center[0], &variables[0], sizeof(center));
    radius = variables[3];
  }

protected:
  virtual void CostFunction(MTL::DynamicVector<T>& residuals, const Parameters& parameters)
  {
    for (U32 i = 0; i < Points_.Size(); i++)
    {
      Point3D<T> center(&parameters[0]);
      T radius = parameters[3];

      residuals[i] = radius - center.Distance(Points_[i]);
    }
  }

  virtual void ComputeJacobian(DynamicMatrix<T>& Jt, const Parameters& currentParameters)
  {
    ComputeJacobianForwardFiniteDifference(Jt, currentParameters);
  }

private:
  DynamicVector<Point3D<T>> Points_;
};

}  // namespace MTL

#endif // MTL_SPHERE_FIT_LEVENBERG_MARQUARDT_H
