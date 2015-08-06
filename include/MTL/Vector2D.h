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


#ifndef MTL_VECTOR_2D_H
#define MTL_VECTOR_2D_H

#include "Matrix.h"
#include "DynamicVectorOperators.h"

namespace MTL
{

template<class T>
class Vector2D : public ColumnVector<2,T>
{
public:
  MTL_COLUMN_VECTOR_COMMON_DEFINITIONS(Vector2D, ColumnVector, 2, T);

  MTL_INLINE Vector2D() : ColumnVector<2,T>() {}
  MTL_INLINE Vector2D(T xx, T yy)
  {
    x(xx);
    y(yy);
  }

  T Length() const  { return FrobeniusNorm(); }

  // Returns unit vector.
  Vector2D unit() const  { return *this / Length(); }

  const T& x() const   { return (*this)[0]; }
  const T& y() const   { return (*this)[1]; }

  void x(const T& xx)  { (*this)[0] = xx;   }
  void y(const T& yy)  { (*this)[1] = yy;   }
};

}  // namespace MTL

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Vector2D<F32>);
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(Vector2D<F64>);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Vector2D<F32>,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(Vector2D<F64>,F64);

#endif // MTL_VECTOR_2D_H
