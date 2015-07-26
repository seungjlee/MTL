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

//
// Need to optimize these.
//
template <class T> MTL_INLINE static
void SquareAll(DynamicVector<T>& v)
{
  FOR_EACH_INDEX(v)
    v[vIndex] = Square(v[vIndex]);
}
template <class T> MTL_INLINE static
void SquareRootAll(DynamicVector<T>& v)
{
  FOR_EACH_INDEX(v)
    v[vIndex] = Sqrt(v[vIndex]);
}

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
// Define stream optimizations.
#define MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(T, CastT)                                  \
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
#else
#define MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(T, CastT)
#endif

//
// Reduction operations.
//
template <class T> MTL_INLINE static T Sum(const DynamicVector<T>& v)
{
  return Sum_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T SumOfAbsolutes(const DynamicVector<T>& v)
{
  return SumOfAbsolutes_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T SumOfSquares(const DynamicVector<T>& v)
{
  return SumOfSquares_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T Min(const DynamicVector<T>& v)
{
  return Min_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T Max(const DynamicVector<T>& v)
{
  return Max_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T MinOfAbsolutes(const DynamicVector<T>& v)
{
  return MinOfAbsolutes_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T MaxOfAbsolutes(const DynamicVector<T>& v)
{
  return MaxOfAbsolutes_Sequential(v.Begin(), v.End());
}
template <class T> MTL_INLINE static T MaxNorm(const DynamicVector<T>& v)
{
  return MaxOfAbsolutes(v);
}
template <class T> MTL_INLINE static T DotProduct(const DynamicVector<T>& v1,
                                                  const DynamicVector<T>& v2)
{
  return DotProduct_Sequential(v1.Begin(), v2.Begin(), v1.End());
}

//
// Some extra operations.
//

template <class T> MTL_INLINE static T Mean(const DynamicVector<T>& v)
{
  return Sum(v) / v.Size();
}

template <class T> MTL_INLINE static T Variance(const DynamicVector<T>& v, const T& mean)
{
  if (v.Size() < 2)
    return T();

  DynamicVector<T> deltas = v;
  deltas -= mean;
  return SumOfSquares(deltas) / (deltas.Size()-1);
}

template <class T> MTL_INLINE static T RMS(const DynamicVector<T>& v)
{
  if (v.Size() > 0)
    return Sqrt(SumOfSquares(v) / v.Size());
  else
    return T(kINF);
}

template <class T> MTL_INLINE static T FrobeniusNorm(const DynamicVector<T>& v)
{
  return Sqrt(SumOfSquares(v));
}


#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
#define MTL_DYNAMIC_VECTOR_STREAM_SUM(T)                                                          \
MTL_INLINE static T Sum(const MTL::DynamicVector<T>& v)                                           \
{                                                                                                 \
  return MTL::Sum_StreamAligned_Parallel(v.Begin(), v.Size());                                    \
}

#define MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_ABSOLUTES(T)                                             \
MTL_INLINE static T SumOfAbsolutes(const MTL::DynamicVector<T>& v)                                \
{                                                                                                 \
  return MTL::SumOfAbsolutes_StreamAligned_Parallel(v.Begin(), v.Size());                         \
}

#define MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_SQUARES(T)                                               \
MTL_INLINE static T SumOfSquares(const MTL::DynamicVector<T>& v)                                  \
{                                                                                                 \
  return MTL::SumOfSquares_StreamAligned_Parallel(v.Begin(), v.Size());                           \
}

#define MTL_DYNAMIC_VECTOR_STREAM_MIN(T)                                                          \
MTL_INLINE static T Min(const MTL::DynamicVector<T>& v)                                           \
{                                                                                                 \
  return MTL::Min_StreamAligned_Parallel(v.Begin(), v.Size());                                    \
}

#define MTL_DYNAMIC_VECTOR_STREAM_MAX(T)                                                          \
MTL_INLINE static T Max(const MTL::DynamicVector<T>& v)                                           \
{                                                                                                 \
  return MTL::Max_StreamAligned_Parallel(v.Begin(), v.Size());                                    \
}

#define MTL_DYNAMIC_VECTOR_STREAM_MIN_OF_ABSOLUTES(T)                                             \
MTL_INLINE static T MinOfAbsolutes(const MTL::DynamicVector<T>& v)                                \
{                                                                                                 \
  return MTL::MinOfAbsolutes_StreamAligned_Parallel(v.Begin(), v.Size());                         \
}

#define MTL_DYNAMIC_VECTOR_STREAM_MAX_OF_ABSOLUTES(T)                                             \
MTL_INLINE static T MaxOfAbsolutes(const MTL::DynamicVector<T>& v)                                \
{                                                                                                 \
  return MTL::MaxOfAbsolutes_StreamAligned_Parallel(v.Begin(), v.Size());                         \
}

#define MTL_DYNAMIC_VECTOR_STREAM_MAX_NORM(T)                                                     \
MTL_INLINE static T MaxNorm(const MTL::DynamicVector<T>& v)                                       \
{                                                                                                 \
  return MTL::MaxOfAbsolutes(v);                                                                  \
}

#define MTL_DYNAMIC_VECTOR_STREAM_DOT_PRODUCT(T)                                                  \
MTL_INLINE static T DotProduct(const MTL::DynamicVector<T>& v1, const MTL::DynamicVector<T>& v2)  \
{                                                                                                 \
  assert(v1.Size() == v2.Size());                                                                 \
  return MTL::DotProduct_StreamAligned_Parallel(v1.Begin(), v2.Begin(), v2.Size());               \
}
#else // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX

#define MTL_DYNAMIC_VECTOR_STREAM_SUM(T)
#define MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_ABSOLUTES(T)
#define MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_SQUARES(T)
#define MTL_DYNAMIC_VECTOR_STREAM_MIN(T)
#define MTL_DYNAMIC_VECTOR_STREAM_MAX(T)
#define MTL_DYNAMIC_VECTOR_STREAM_MIN_OF_ABSOLUTES(T)
#define MTL_DYNAMIC_VECTOR_STREAM_MAX_OF_ABSOLUTES(T)
#define MTL_DYNAMIC_VECTOR_STREAM_MAX_NORM(T)
#define MTL_DYNAMIC_VECTOR_STREAM_DOT_PRODUCT(T)

#endif // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX


#define MTL_DYNAMIC_VECTOR_STREAM_REDUCTION_OPERATIONS(T)  \
MTL_DYNAMIC_VECTOR_STREAM_SUM(T)                           \
MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_SQUARES(T)                \
MTL_DYNAMIC_VECTOR_STREAM_SUM_OF_ABSOLUTES(T)              \
MTL_DYNAMIC_VECTOR_STREAM_MIN(T)                           \
MTL_DYNAMIC_VECTOR_STREAM_MAX(T)                           \
MTL_DYNAMIC_VECTOR_STREAM_MIN_OF_ABSOLUTES(T)              \
MTL_DYNAMIC_VECTOR_STREAM_MAX_OF_ABSOLUTES(T)              \
MTL_DYNAMIC_VECTOR_STREAM_MAX_NORM(T)                      \
MTL_DYNAMIC_VECTOR_STREAM_DOT_PRODUCT(T)

// Define optimizations for valid types.
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(F32,F32);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(F64,F64);
MTL_DYNAMIC_VECTOR_STREAM_REDUCTION_OPERATIONS(F32);
MTL_DYNAMIC_VECTOR_STREAM_REDUCTION_OPERATIONS(F64);

MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector1D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector2D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector3D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector4D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector5D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(ColumnVector6D,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix1x1,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix2x2,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix3x3,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix4x4,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix5x5,F64);
MTL_DYNAMIC_VECTOR_STREAM_PARALLEL_OPERATIONS(SquareMatrix6x6,F64);

}  // namespace MTL


#endif // MTL_DYNAMIC_VECTOR_OPERATORS_H