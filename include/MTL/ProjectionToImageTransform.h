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
  ProjectionToImageTransform(const MatrixType& m = MatrixType::eIdentity,
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

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ProjectionToImageTransform<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ProjectionToImageTransform<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ProjectionToImageTransform<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ProjectionToImageTransform<F64>,F64);

}  // namespace MTL

#endif  // MTL_PROJECTION_TO_IMAGE_TRANSFORM_3D_H
