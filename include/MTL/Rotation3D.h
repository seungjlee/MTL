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
  MTL_SQUARE_MATRIX_COMMON_DEFINITIONS(Rotation3D, Base, 3, T);

  Rotation3D() {}

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

  MTL_INLINE Rotation3D Inverse() const  { return Base::ComputeTranspose(); }
  MTL_INLINE void Invert()               { Base::Transpose();               }

  MTL_INLINE Point3D<T> operator*(const Point3D<T>& pt) const
  { return Base::operator*(pt); }

  MTL_INLINE Vector3D<T> operator*(const Vector3D<T>& pt) const
  { return Base::operator*(pt); }

  //
  // Computing rotations to/from a different axis using Given rotations.
  //
  MTL_INLINE static Rotation3D RotationToXAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToXAxis(c1, c2, s1, s2, v);

    return Rotation3D( c1*c2, -s2, -c2*s1,
                       c1*s2,  c2, -s1*s2,
                          s1,   0,     c1 );
  }

  MTL_INLINE static Rotation3D RotationFromXAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToXAxis(c1, c2, s1, s2, v);

    return Rotation3D(  c1*c2,  c1*s2, s1,
                          -s2,     c2,  0,
                       -c2*s1, -s1*s2, c1 );
  }

  MTL_INLINE static Rotation3D RotationToYAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToYAxis(c1, c2, s1, s2, v);

    return Rotation3D(-s2, -c1*c2, s1*c2,
                       c2, -c1*s2, s1*s2,
                        0,     s1,    c1 );
  }

  MTL_INLINE static Rotation3D RotationFromYAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToYAxis(c1, c2, s1, s2, v);

    return Rotation3D(    -s2,     c2, 0,
                       -c1*c2, -c1*s2, s1,
                        s1*c2,  s1*s2, c1 );
  }

  MTL_INLINE static Rotation3D RotationToZAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToZAxis(c1, c2, s1, s2, v);

    return Rotation3D(    -s1,   0,   -c1,
                       -c1*c2, -s2, c2*s1,
                       -c1*s2,  c2, s1*s2 );
  }

  MTL_INLINE static Rotation3D RotationFromZAxis(const Vector3D<T>& v)
  {
    if (v.SumOfSquares() < EpsilonSquared<T>())
      return Rotation3D();

    T c1, c2, s1, s2;
    RotationToZAxis(c1, c2, s1, s2, v);

    return Rotation3D( -s1, -c1*c2, -c1*s2,
                         0,    -s2,     c2,
                       -c1,  c2*s1,  s1*s2 );
  }

private:
  MTL_INLINE static void RotationToXAxis(T& c1, T& c2, T& s1, T& s2,
                                         const Vector3D<T>& v)
  {
    if (v[0] == 0 && v[2] == 0)
    {
      assert(v[1] != 0);

      s1 =  0;
      s2 = -1;
      c1 =  1;
      c2 =  0;
    }
    else
    {
      T r = Hypotenuse(v[0], v[2]);
      c1 =  v[0]/r;
      s1 = -v[2]/r;

      T w0 = c1 * v[0] - s1 * v[2];

      r = Hypotenuse(w0, v[1]);
      c2 = w0/r;
      s2 =  -v[1]/r;
    }
  }

  MTL_INLINE static void RotationToYAxis(T& c1, T& c2, T& s1, T& s2,
                                         const Vector3D<T>& v)
  {
    if (v[1] == 0 && v[2] == 0)
    {
      assert(v[0] != 0);

      s1 = 0;
      s2 = 0;
      c1 = 1;
      c2 = 1;
    }
    else
    {
      T r = Hypotenuse(v[1], v[2]);
      c1 =  v[1]/r;
      s1 = -v[2]/r;

      T w1 = c1 * v[1] - s1 * v[2];

      r = Hypotenuse(w1, v[0]);
      c2 = v[0]/r;
      s2 =  -w1/r;
    }
  }

  MTL_INLINE static void RotationToZAxis(T& c1, T& c2, T& s1, T& s2,
                                         const Vector3D<T>& v)
  {
    if (v[0] == 0 && v[2] == 0)
    {
      assert(v[1] != 0);

      s1 = -1;
      s2 =  0;
      c1 =  0;
      c2 =  1;
    }
    else
    {
      T r = Hypotenuse(v[0], v[2]);
      c1 =  v[0]/r;
      s1 = -v[2]/r;

      T w2 = c1 * v[0] - s1 * v[2];

      r = Hypotenuse(w2, v[1]);
      c2 = v[1]/r;
      s2 =  -w2/r;
    }
  }

  MTL_INLINE Rotation3D(double e00, double e01, double e02,
                        double e10, double e11, double e12,
                        double e20, double e21, double e22)
  {
    Base::Data_[0][0] = e00;  Base::Data_[0][1] = e01;  Base::Data_[0][2] = e02;
    Base::Data_[1][0] = e10;  Base::Data_[1][1] = e11;  Base::Data_[1][2] = e12;
    Base::Data_[2][0] = e20;  Base::Data_[2][1] = e21;  Base::Data_[2][2] = e22;
  }
};

}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(MTL::Rotation3D<MTL::F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(MTL::Rotation3D<MTL::F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(MTL::Rotation3D<MTL::F32>,MTL::F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(MTL::Rotation3D<MTL::F64>,MTL::F64);

#endif // MTL_ROTATION_3D_H