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


#ifndef MTL_POINT_2D_H
#define MTL_POINT_2D_H

#include "Matrix.h"
#include "DynamicVectorOperators.h"

namespace MTL
{

template<class T>
class Point2D : public ColumnVector<2,T>
{
  typedef ColumnVector<2,T> Base;

public:
  MTL_COLUMN_VECTOR_COMMON_DEFINITIONS(Point2D, Base, 2, T);

  MTL_INLINE Point2D() : ColumnVector<2,T>() {}
  MTL_INLINE Point2D(T xx, T yy)
  {
    x(xx);
    y(yy);
  }

  T DistanceSquared(const Point2D& pt) const
  {
    Point2D delta = pt - *this;
    return delta.SumOfSquares();
  }

  T Distance(const Point2D& pt) const
  {
    return Sqrt(DistanceSquared(pt));
  }

  T Length() const  { return Base::FrobeniusNorm(); }

  const T& x() const   { return (*this)[0]; }
  const T& y() const   { return (*this)[1]; }

  void x(const T& xx)  { (*this)[0] = xx;   }
  void y(const T& yy)  { (*this)[1] = yy;   }
};

template<class T>
MTL_INLINE static Point2D<T> Mean(const DynamicVector<Point2D<T>>& v)
{
  return Sum(v) / T(v.Size());
}


}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Point2D<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Point2D<F64>,F64);

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<I32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<U32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<I16>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Point2D<U16>);

#endif // MTL_POINT_2D_H
