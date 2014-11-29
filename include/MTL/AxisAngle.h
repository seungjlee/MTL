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

#include "Vector3D.h"

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

  MTL_INLINE AxisAngle Inverse() const  { return AxisAngle(-RotationVector_); }

private:
  Vector3D<T> RotationVector_;

  // Cached values.
  T Angle_;
  T SinHalfAngle_;
  T CosHalfAngle_;
  T SinAngle_;
  T CosAngle_;
  Vector3D<T> UnitRotationVector_;
  bool DirtyCachedValues_;
};

}  // namespace MTL

#endif // MTL_AXIS_ANGLE_H
