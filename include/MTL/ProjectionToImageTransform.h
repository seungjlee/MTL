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


#ifndef MTL_PROJECTION_TO_IMAGE_TRANSFORM_3D_H
#define MTL_PROJECTION_TO_IMAGE_TRANSFORM_3D_H

#include "AffineTransform3D.h"
#include "Point2D.h"

namespace MTL
{

template <class T>
class ProjectionToImageTransform : public AffineTransform3D<T>
{
public:
  ProjectionToImageTransform(const typename AffineTransform3D<T>::MatrixType& m =
                             AffineTransform3D<T>::MatrixType::eIdentity,
                             const Vector3D<T>& v = Vector3D<T>(0,0,0))
    : AffineTransform3D<T>(m, v)
  {
  }

  MTL_INLINE Point2D<T> operator*(const Point3D<T>& point) const
  {
    Point3D<T> pt3D = AffineTransform3D<T>::operator*(point);

    return Point2D<T>(pt3D.x() / pt3D.z(), pt3D.y() / pt3D.z());
  }
};

/*

From http://www.robots.ox.ac.uk/~vgg/hzbook/code/

%F = vgg_F_from_P(P)  Compute fundamental matrix from two camera matrices.
%   P is cell (2), P = {P1 P2}. F has size (3,3). It is x2'*F*x1 = 0
%
%   Overall scale of F is unique and such that, for any X, P1, P2, it is
%   F*x1 = vgg_contreps(e2)*x2, where
%   x1 = P1*X, x2 = P2*X, e2 = P2*C1, C1 = vgg_wedge(P1).

function F = vgg_F_from_P(P, P2)

if nargin == 1
  P1 = P{1};
  P2 = P{2};
else
  P1 = P;
end

X1 = P1([2 3],:);
X2 = P1([3 1],:);
X3 = P1([1 2],:);
Y1 = P2([2 3],:);
Y2 = P2([3 1],:);
Y3 = P2([1 2],:);

F = [det([X1; Y1]) det([X2; Y1]) det([X3; Y1])
     det([X1; Y2]) det([X2; Y2]) det([X3; Y2])
     det([X1; Y3]) det([X2; Y3]) det([X3; Y3])];

return
*/
#ifdef WIN32  // Disabled for g++ which has trouble with templates.
template <class T>
static SquareMatrix<3,T> ComputeFundamentalMatrix(const ProjectionToImageTransform<T>& P1,
                                                  const ProjectionToImageTransform<T>& P2)
{
  Matrix<2,4,T> X1 =  P1.Matrix().SubMatrix<1,0,2,3>() || P1.Vector().SubColumn<1,2>();
  Matrix<2,4,T> X2 = (P1.Matrix().SubMatrix<2,0,1,3>() || P1.Vector().SubColumn<2,1>()) &&
                     (P1.Matrix().SubMatrix<0,0,1,3>() || P1.Vector().SubColumn<0,1>());
  Matrix<2,4,T> X3 =  P1.Matrix().SubMatrix<0,0,2,3>() || P1.Vector().SubColumn<0,2>();
  Matrix<2,4,T> Y1 =  P2.Matrix().SubMatrix<1,0,2,3>() || P2.Vector().SubColumn<1,2>();
  Matrix<2,4,T> Y2 = (P2.Matrix().SubMatrix<2,0,1,3>() || P2.Vector().SubColumn<2,1>()) &&
                     (P2.Matrix().SubMatrix<0,0,1,3>() || P2.Vector().SubColumn<0,1>());
  Matrix<2,4,T> Y3 =  P2.Matrix().SubMatrix<0,0,2,3>() || P2.Vector().SubColumn<0,2>();

  SquareMatrix<3,T> F;

  F[0][0] = StableDeterminant(X1 && Y1);
  F[0][1] = StableDeterminant(X2 && Y1);
  F[0][2] = StableDeterminant(X3 && Y1);
  F[1][0] = StableDeterminant(X1 && Y2);
  F[1][1] = StableDeterminant(X2 && Y2);
  F[1][2] = StableDeterminant(X3 && Y2);
  F[2][0] = StableDeterminant(X1 && Y3);
  F[2][1] = StableDeterminant(X2 && Y3);
  F[2][2] = StableDeterminant(X3 && Y3);

  return F;
}
#endif

}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ProjectionToImageTransform<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ProjectionToImageTransform<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ProjectionToImageTransform<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ProjectionToImageTransform<F64>,F64);

#endif  // MTL_PROJECTION_TO_IMAGE_TRANSFORM_3D_H
