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


#ifndef MTL_POINT_3D_H
#define MTL_POINT_3D_H

#include "Matrix.h"
#include "DynamicVectorOperators.h"

namespace MTL
{

template<class T>
class Point3D : public ColumnVector<3,T>
{
  typedef ColumnVector<3,T> Base;

public:
  MTL_COLUMN_VECTOR_COMMON_DEFINITIONS(Point3D, Base, 3, T);

  MTL_INLINE Point3D() : ColumnVector<3,T>() {}
  MTL_INLINE Point3D(T xx, T yy, T zz)
  {
    x(xx);
    y(yy);
    z(zz);
  }

  T DistanceSquared(const Point3D& pt) const
  {
    Point3D delta = pt - *this;
    return delta.SumOfSquares();
  }

  T Distance(const Point3D& pt) const
  {
    return Sqrt(DistanceSquared(pt));
  }

  T Length() const  { return Base::FrobeniusNorm(); }

  const T& x() const   { return (*this)[0]; }
  const T& y() const   { return (*this)[1]; }
  const T& z() const   { return (*this)[2]; }

  void x(const T& xx)  { (*this)[0] = xx;   }
  void y(const T& yy)  { (*this)[1] = yy;   }
  void z(const T& yy)  { (*this)[2] = yy;   }
};

template<class T>
MTL_INLINE static Point3D<T> Mean(const DynamicVector<Point3D<T>>& v)
{
  return Sum(v) / T(v.Size());
}


}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point3D<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point3D<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Point3D<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Point3D<F64>,F64);

#endif // MTL_POINT_3D_H