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


#ifndef MTL_ROTATION_3D_H
#define MTL_ROTATION_3D_H

#include "Point3D.h"
#include "Vector3D.h"
#include "Trigonometry.h"

namespace MTL
{

template <class T>
class Rotation3D : public SquareMatrix<3,T>
{
  typedef SquareMatrix<3,T> Base;

public:
  MTL_SQUARE_MATRIX_COMMON_DEFINITIONS(Rotation3D, SquareMatrix, 3, T);

  Rotation3D() : SquareMatrix<3,T>(eIdentity) {}

  MTL_INLINE Rotation3D(T roll, T pitch, T yaw) 
  {
    *this = RotationX(roll) * RotationY(pitch) * RotationZ(yaw);
  }

  MTL_INLINE static Rotation3D RotationX(T angle)
  {
    T c, s;
    ComputeCosineSine(c, s, angle);
    return RotationX(c, s);
  }
  MTL_INLINE static Rotation3D RotationY(T angle)
  {
    T c, s;
    ComputeCosineSine(c, s, angle);
    return RotationY(c, s);
  }
  MTL_INLINE static Rotation3D RotationZ(T angle)
  {
    T c, s;
    ComputeCosineSine(c, s, angle);
    return RotationZ(c, s);
  }

  MTL_INLINE static Rotation3D RotationX(T c, T s)
  {
    T rotX[3][3] = { { 1, 0,  0 },
                     { 0, c, -s },
                     { 0, s,  c } };

    return Rotation3D(rotX);
  }
  MTL_INLINE static Rotation3D RotationY(T c, T s)
  {
    T rotY[3][3] = { {  c, 0, s },
                     {  0, 1, 0 },
                     { -s, 0, c } };

    return Rotation3D(rotY);
  }
  MTL_INLINE static Rotation3D RotationZ(T c, T s)
  {
    T rotZ[3][3] = { { c, -s, 0 },
                     { s,  c, 0 },
                     { 0,  0, 1 } };

    return Rotation3D(rotZ);
  }

  MTL_INLINE Rotation3D Inverse() const  { return ComputeTranspose(); }
  MTL_INLINE void Invert()               { Transpose();               }

  MTL_INLINE Point3D<T> operator*(const Point3D<T>& pt) const
  { return Base::operator*(pt); }

  MTL_INLINE Vector3D<T> operator*(const Vector3D<T>& pt) const
  { return Base::operator*(pt); }
};

}  // namespace MTL

#endif // MTL_ROTATION_3D_H
