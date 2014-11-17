//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee
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


#ifndef MTL_DYNAMIC_VECTOR_OPERATORS_H
#define MTL_DYNAMIC_VECTOR_OPERATORS_H

#include "DynamicVector.h"

namespace MTL
{
// Helper class to override with optimized code.
template <class T>
class DynamicVectorParallelOperations
{
public:
  MTL_INLINE static DynamicVector<T> UnaryMinus(const DynamicVector<T>& v)
  {
    DynamicVector<T> negated(v.Size());
    FOR_EACH_INDEX(v)
      negated[vIndex] = -v[vIndex];

    return negated;
  }
  MTL_INLINE static void Addition(DynamicVector<T>& v1, const DynamicVector<T>& v2)
  {
    assert(v1.Size() == v2.Size());
    FOR_EACH_INDEX(v1)
      v1[v1Index] += v2[v1Index];
  }
  MTL_INLINE static void Subtraction(DynamicVector<T>& v1, const DynamicVector<T>& v2)
  {
    assert(v1.Size() == v2.Size());
    FOR_EACH_INDEX(v1)
      v1[v1Index] -= v2[v1Index];
  }
  MTL_INLINE static void Multiplication(DynamicVector<T>& v1, const DynamicVector<T>& v2)
  {
    assert(v1.Size() == v2.Size());
    FOR_EACH_INDEX(v1)
      v1[v1Index] *= v2[v1Index];
  }
  MTL_INLINE static void Division(DynamicVector<T>& v1, const DynamicVector<T>& v2)
  {
    assert(v1.Size() == v2.Size());
    FOR_EACH_INDEX(v1)
      v1[v1Index] /= v2[v1Index];
  }
  MTL_INLINE static void ScalarAddition(DynamicVector<T>& v, double s)
  {
    FOR_EACH_INDEX(v)
      v[vIndex] = T(v[vIndex] + s);
  }
  MTL_INLINE static void ScalarSubtraction(DynamicVector<T>& v, double s)
  {
    FOR_EACH_INDEX(v)
      v[vIndex] = T(v[vIndex] - s);
  }
  MTL_INLINE static void ScalarMultiplication(DynamicVector<T>& v, double s)
  {
    FOR_EACH_INDEX(v)
      v[vIndex] = T(v[vIndex] * s);
  }
  MTL_INLINE static void ScalarDivision(DynamicVector<T>& v, double s)
  {
    FOR_EACH_INDEX(v)
      v[vIndex] = T(v[vIndex] / s);
  }
};

template <class T> MTL_INLINE static DynamicVector<T> operator-(const DynamicVector<T>& v)
{
  return DynamicVectorParallelOperations<T>::UnaryMinus(v);
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator+=(DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVectorParallelOperations<T>::Addition(v1, v2);
  return v1;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator+(const DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVector<T> sum(v1);
  sum += v2;
  return sum;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator-=(DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVectorParallelOperations<T>::Subtraction(v1, v2);
  return v1;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator-(const DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVector<T> sum(v1);
  sum -= v2;
  return sum;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator*=(DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVectorParallelOperations<T>::Multiplication(v1, v2);
  return v1;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator*(const DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVector<T> result(v1);
  result *= v2;
  return result;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator/=(DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVectorParallelOperations<T>::Division(v1, v2);
  return v1;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator/(const DynamicVector<T>& v1, const DynamicVector<T>& v2)
{
  DynamicVector<T> result(v1);
  result /= v2;
  return result;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator+=(DynamicVector<T>& v, double s)
{
  DynamicVectorParallelOperations<T>::ScalarAddition(v, s);
  return v;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator+(const DynamicVector<T>& v, double s)
{
  DynamicVector<T> u(v);
  u += s;
  return u;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator+(double& s, const DynamicVector<T>& v)
{
  return v + s;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator-=(DynamicVector<T>& v, double s)
{
  DynamicVectorParallelOperations<T>::ScalarSubtraction(v, s);
  return v;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator-(const DynamicVector<T>& v, double s)
{
  DynamicVector<T> u(v);
  u -= s;
  return u;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator*=(DynamicVector<T>& v, double s)
{
  DynamicVectorParallelOperations<T>::ScalarMultiplication(v, s);
  return v;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator*(const DynamicVector<T>& v, double s)
{
  DynamicVector<T> u(v);
  u *= s;
  return u;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator*(double s, const DynamicVector<T>& v)
{
  return v * s;
}

template <class T> MTL_INLINE static
DynamicVector<T>& operator/=(DynamicVector<T>& v, double s)
{
  DynamicVectorParallelOperations<T>::ScalarDivision(v, s);
  return v;
}
template <class T> MTL_INLINE static
DynamicVector<T> operator/(const DynamicVector<T>& v, double s)
{
  DynamicVector<T> u(v);
  u /= s;
  return u;
}

// Define stream optimizations.
#define MTL_DYNAMIC_VECTOR_SSE_PARALLEL_OPERATIONS(T, CastT)                                     \
template<> class MTL::DynamicVectorParallelOperations<T>                                         \
{                                                                                                \
public:                                                                                          \
  MTL_INLINE static MTL::DynamicVector<T> UnaryMinus(const MTL::DynamicVector<T>& v)             \
  {                                                                                              \
    MTL::DynamicVector<T> negated = v;                                                           \
    MTL::UnaryMinus_StreamAligned_Parallel((CastT*)negated.Begin(),                              \
                                           negated.Size()*sizeof(T)/sizeof(CastT));              \
    return negated;                                                                              \
  }                                                                                              \
  MTL_INLINE static void Addition(MTL::DynamicVector<T>& v1,                                     \
                                  const MTL::DynamicVector<T>& v2)                               \
  {                                                                                              \
    assert(v1.Size() == v2.Size());                                                              \
    MTL::Addition_StreamAligned_Parallel((CastT*)v1.Begin(), (CastT*)v2.Begin(),                 \
                                         v2.Size()*sizeof(T)/sizeof(CastT));                     \
  }                                                                                              \
  MTL_INLINE static void Subtraction(MTL::DynamicVector<T>& v1,                                  \
                                     const MTL::DynamicVector<T>& v2)                            \
  {                                                                                              \
    assert(v1.Size() == v2.Size());                                                              \
    MTL::Subtraction_StreamAligned_Parallel((CastT*)v1.Begin(), (CastT*)v2.Begin(),              \
                                            v2.Size()*sizeof(T)/sizeof(CastT));                  \
  }                                                                                              \
  MTL_INLINE static void Multiplication(MTL::DynamicVector<T>& v1,                               \
                                        const MTL::DynamicVector<T>& v2)                         \
  {                                                                                              \
    assert(v1.Size() == v2.Size());                                                              \
    MTL::Multiplication_StreamAligned_Parallel((CastT*)v1.Begin(), (CastT*)v2.Begin(),           \
                                                v2.Size()*sizeof(T)/sizeof(CastT));              \
  }                                                                                              \
  MTL_INLINE static void Division(MTL::DynamicVector<T>& v1,                                     \
                                  const MTL::DynamicVector<T>& v2)                               \
  {                                                                                              \
    assert(v1.Size() == v2.Size());                                                              \
    MTL::Division_StreamAligned_Parallel((CastT*)v1.Begin(), (CastT*)v2.Begin(),                 \
                                         v2.Size()*sizeof(T)/sizeof(CastT));                     \
  }                                                                                              \
  MTL_INLINE static void ScalarAddition(MTL::DynamicVector<T>& v, double s)                      \
  {                                                                                              \
    MTL::ScalarAddition_StreamAligned_Parallel((CastT*)v.Begin(), (CastT)s,                      \
                                               v.Size()*sizeof(T)/sizeof(CastT));                \
  }                                                                                              \
  MTL_INLINE static void ScalarSubtraction(MTL::DynamicVector<T>& v, double s)                   \
  {                                                                                              \
    MTL::ScalarSubtraction_StreamAligned_Parallel((CastT*)v.Begin(), (CastT)s,                   \
                                                  v.Size()*sizeof(T)/sizeof(CastT));             \
  }                                                                                              \
  MTL_INLINE static void ScalarMultiplication(MTL::DynamicVector<T>& v, double s)                \
  {                                                                                              \
    MTL::ScalarMultiplication_StreamAligned_Parallel((CastT*)v.Begin(), (CastT)s,                \
                                                     v.Size()*sizeof(T)/sizeof(CastT));          \
  }                                                                                              \
  MTL_INLINE static void ScalarDivision(MTL::DynamicVector<T>& v, double s)                      \
  {                                                                                              \
    MTL::ScalarDivision_StreamAligned_Parallel((CastT*)v.Begin(), (CastT)s,                      \
                                                v.Size()*sizeof(T)/sizeof(CastT));               \
  }                                                                                              \
};

MTL_DYNAMIC_VECTOR_SSE_PARALLEL_OPERATIONS(F32,F32);
MTL_DYNAMIC_VECTOR_SSE_PARALLEL_OPERATIONS(F64,F64);

}  // namespace MTL


#endif // MTL_DYNAMIC_VECTOR_OPERATORS_H
