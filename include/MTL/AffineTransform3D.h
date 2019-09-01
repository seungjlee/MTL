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


#ifndef MTL_AFFINE_TRANSFORM_3D_H
#define MTL_AFFINE_TRANSFORM_3D_H

#include "Point3D.h"
#include "Vector3D.h"

namespace MTL
{

template <class T>
class AffineTransform3D
{
public:
  typedef SquareMatrix<3,T> MatrixType;
  typedef Vector3D<T> VectorType;

  AffineTransform3D(const MatrixType& m, const Vector3D<T>& v)
    : Matrix_(m), Vector_(v)
  {
  }

  MTL_INLINE AffineTransform3D operator*(const AffineTransform3D& transform) const
  {
    return AffineTransform3D(Matrix_ * transform.Matrix_, Matrix_ * transform.Vector_ + Vector_);
  }
  MTL_INLINE AffineTransform3D& operator*=(const AffineTransform3D& transform)
  {
    *this = *this * transform;
    return *this;
  }

  AffineTransform3D Inverse() const
  {
    MatrixType inv = Matrix_.Inverse();
    return AffineTransform3D(inv, inv * -Vector_);
  }

  MTL_INLINE Point3D<T> operator*(const Point3D<T>& point) const
  {
    return Matrix_ * point + Vector_;
  }
  MTL_INLINE Vector3D<T> operator*(const Vector3D<T>& vector) const
  {
    return Matrix_ * vector;
  }

  MTL_INLINE void ScaleVector(T factor)  { Vector_ *= factor;    }

  MTL_INLINE const MatrixType& Matrix() const  { return Matrix_; }
  MTL_INLINE void Matrix(const MatrixType& m)  { Matrix_ = m;    }
  MTL_INLINE const VectorType& Vector() const  { return Vector_; }
  MTL_INLINE void Vector(const VectorType& v)  { Vector_ = v;    }

protected:
  MatrixType Matrix_;
  VectorType Vector_;  // Translation vector.
};

}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(AffineTransform3D<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(AffineTransform3D<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(AffineTransform3D<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(AffineTransform3D<F64>,F64);

#endif  // MTL_AFFINE_TRANSFORM_3D_H
