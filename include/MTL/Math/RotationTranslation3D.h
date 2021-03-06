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


#ifndef MTL_ROTATION_TRANSLATION_3D_H
#define MTL_ROTATION_TRANSLATION_3D_H

#include "AffineTransform3D.h"
#include "Rotation3D.h"

namespace MTL
{

template <class T>
class RotationTranslation3D : public AffineTransform3D<T>
{
public:
  RotationTranslation3D(const Rotation3D<T>& r = Rotation3D<T>(Rotation3D<T>::eIdentity),
                        const Vector3D<T>& v = Vector3D<T>(0,0,0))
    : AffineTransform3D<T>(r, v)
  {
  }

  MTL_INLINE RotationTranslation3D operator*(const RotationTranslation3D& transform) const
  {
    return RotationTranslation3D(this->Matrix_ * transform.Matrix_,
                                 this->Matrix_ * transform.Vector_ + this->Vector_);
  }
  MTL_INLINE RotationTranslation3D& operator*=(const RotationTranslation3D& transform)
  {
    *this = *this * transform;
    return *this;
  }

  MTL_INLINE Point3D<T> operator*(const Point3D<T>& point) const
  {
    return this->Matrix_ * point + this->Vector_;
  }
  MTL_INLINE Vector3D<T> operator*(const Vector3D<T>& vector) const
  {
    return this->Matrix_ * vector;
  }

  RotationTranslation3D<T> Inverse() const
  {
    typename AffineTransform3D<T>::MatrixType
      inv = AffineTransform3D<T>::Matrix_.ComputeTranspose();

    return RotationTranslation3D<T>(inv, inv * -AffineTransform3D<T>::Vector_);
  }
};

}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(RotationTranslation3D<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(RotationTranslation3D<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(RotationTranslation3D<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(RotationTranslation3D<F64>,F64);

#endif // MTL_ROTATION_TRANSLATION_3D_H
