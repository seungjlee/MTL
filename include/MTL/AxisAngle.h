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


#ifndef MTL_AXIS_ANGLE_H
#define MTL_AXIS_ANGLE_H

#include "Rotation3D.h"
#include "SVD.h"
#include "Trigonometry.h"

namespace MTL
{

// Compact 3D axis angle rotation representation class.
template<class T>
class AxisAngle
{
public:
  AxisAngle()
    : RotationVector_(0,0,0), DirtyCachedValues_(true) {}
  AxisAngle(const Vector3D<T>& rotationVector)
    : RotationVector_(rotationVector), DirtyCachedValues_(true) {}

  AxisAngle(const Rotation3D<T>& R)
    : DirtyCachedValues_(true)
  {
    // From Multiple View Geometry (2nd Edition), R. Hartley and A. Zisserman (A4.10).
    Rotation3D<T> A = R - Rotation3D<T>();

    long rank;
    T conditionNumber;
    SolveJacobiSVDHomogeneous(A, RotationVector_, rank, conditionNumber);

    Vector3D<T> v(R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]);

    T angle = atan2(RotationVector_.Dot(v), R[0][0] + R[1][1] + R[2][2] - 1);

    RotationVector_ *= angle;
  }

  MTL_INLINE AxisAngle Inverse() const  { return AxisAngle(-RotationVector_); }

  MTL_INLINE Vector3D<T> operator*(const Vector3D<T>& pt) const
  {
    const_cast<AxisAngle*>(this)->ComputeCachedValues();

    // Use Rodrigues' rotation formula.
    return pt * CosAngle_ + Cross(UnitRotationVector_, pt) * SinAngle_ + 
           UnitRotationVector_ * (UnitRotationVector_.Dot(pt) * (1 - CosAngle_));
  }

  MTL_INLINE AxisAngle operator*(const AxisAngle& other) const
  {
    return AxisAngle(Multiply(other));
  }

  MTL_INLINE AxisAngle& operator*=(const AxisAngle& other)
  {
    *this = *this * other;
    return *this;
  }

  MTL_INLINE void GetUnitQuaternions(T& qs, Vector3D<T>& qv) const
  {
    const_cast<AxisAngle*>(this)->computeCachedValues();
    qs = CosHalfAngle_;
    qv = UnitRotationVector_ * SinHalfAngle_;
  }

  MTL_INLINE bool IsIdentity() const  { return RotationVector_ == vector3D<T>(0,0,0); }

  MTL_INLINE Rotation3D<T> GetRotationMatrix() const
  {
    const_cast<AxisAngle*>(this)->ComputeCachedValues();

    if (SinAngle_ == 0 && CosAngle_ > 0)
      return Rotation3D<T>(Rotation3D<T>::eIdentity);

    const Vector3D<T>& u = UnitRotationVector_;
    T m1[3][3] = {{          cosAngle_, -u.z() * sinAngle_,  u.y() * sinAngle_ },
                  {  u.z() * sinAngle_,          cosAngle_, -u.x() * sinAngle_ },
                  { -u.y() * sinAngle_,  u.x() * sinAngle_,          cosAngle_ }};

    T xx = u.x() * u.x();
    T xy = u.x() * u.y();
    T xz = u.x() * u.z();
    T yy = u.y() * u.y();
    T yz = u.y() * u.z();
    T zz = u.z() * u.z();
    T m2[3][3] = {{ xx, xy, xz },
                  { xy, yy, yz },
                  { xz, yz, zz }};

    return Rotation3D<T>(m1) + Rotation3D<T>(m2) * (1 - cosAngle_);
  }

private:
  Vector3D<T> RotationVector_;

  // Cached values.
  T SinHalfAngle_;
  T CosHalfAngle_;
  T SinAngle_;
  T CosAngle_;
  T Angle_;
  Vector3D<T> UnitRotationVector_;
  bool DirtyCachedValues_;

  MTL_INLINE Vector3D<T> Multiply(const AxisAngle& other) const
  {
    // Compute product of unit quaternions.
    T    qs1, qs2;
    Vector3DD qv1, qv2;
    GetUnitQuaternions(qs1, qv1);
    other.GetUnitQuaternions(qs2, qv2);

    T qs = qs1*qs2 - qv1.dot(qv2);

    if (Abs(qs) >= 1)
    {
      return Vector3DD(0, 0, 0);
    }
    else
    {
      T angle = 2*acos(qs);
      T sinHalfAngle = sin(0.5*angle);
      Vector3D<T> qv = qv1*qs2 + qv2*qs1 + qv1.Cross(qv2);
      Vector3D<T> rotationVector = qv * angle / sinHalfAngle;

      return rotationVector;
    }
  }

  MTL_INLINE void ComputeCachedValues()
  {
    if (DirtyCachedValues_)
    {
      Angle_ = RotationVector_.FrobeniusNorm();
      if (Angle_ == T(0.0))
      {
        UnitRotationVector_ = Vector3D<T>(0,0,0);
        SinAngle_ = T(0.0);
        CosAngle_ = T(1.0);
        SinHalfAngle_ = T(0.0);
        CosHalfAngle_ = T(1.0);
      }
      else
      {
        UnitRotationVector_ = RotationVector_ / Angle_;

        T halfAngle = T(0.5) * Angle_;
        ComputeCosineSine(CosHalfAngle_, SinHalfAngle_, halfAngle);

        SinAngle_ = T(2.0) * CosHalfAngle_ * SinHalfAngle_;
        CosAngle_ = Square(CosHalfAngle_) - Square(SinHalfAngle_);
      }

      DirtyCachedValues_ = false;
    }
  }
};

}  // namespace MTL

#endif // MTL_AXIS_ANGLE_H
