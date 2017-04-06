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

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX

#ifndef MTL_STREAM_ARRAY_H
#define MTL_STREAM_ARRAY_H

#include <MTL/SSE.h>
#include <MTL/AVX.h>
#include <MTL/Array.h>
#include <MTL/OpenMP.h>

namespace MTL
{

//
// Regular single-threaded operations.
//
template <class T> MTL_INLINE void
UnaryMinus_Sequential(T* pDst, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst = -*pDst;
}
template <class T> MTL_INLINE void
ScalarAddition_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst += scalar;
}
template <class T> MTL_INLINE void
ScalarSubtraction_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst -= scalar;
}
template <class T> MTL_INLINE void
ScalarMultiplication_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst *= scalar;
}
template <class T> MTL_INLINE void
ScalarDivision_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst /= scalar;
}
template <class DstT, class SrcT> MTL_INLINE void
Addition_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst + *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE void
Subtraction_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst - *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE void
Multiplication_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst * *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE void
Division_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst / *pSrc);
}
template <class T> MTL_INLINE void
Min_Sequential(T* pDst, const T* pSrc, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = Min(*pDst, *pSrc);
}
template <class T> MTL_INLINE void
Max_Sequential(T* pDst, const T* pSrc, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = Max(*pDst, *pSrc);
}
template <class DstT, class SrcT> void
Convert_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pSrc);
}

template <class T> MTL_INLINE T Sum_Sequential(const T* p, const T* pEnd)
{
  if (p >= pEnd)
    return T();

  T sum = *p++;
  for (; p < pEnd; p++)
    sum += *p;

  return sum;
}
template <class T> MTL_INLINE T SumOfAbsolutes_Sequential(const T* p, const T* pEnd)
{
  if (p >= pEnd)
    return T();

  T sum = Abs(*p++);
  for (; p < pEnd; p++)
    sum += Abs(*p);

  return sum;
}
template <class T> MTL_INLINE T SumOfSquares_Sequential(const T* p, const T* pEnd)
{
  if (p >= pEnd)
    return T();

  T sum = Square(*p++);
  for (; p < pEnd; p++)
    sum = MultiplyAndAdd(*p, *p, sum);

  return sum;
}
template <class T> MTL_INLINE T SumOfSquaredDifferences_Sequential(const T* p1, const T* p2,
                                                                   const T* pEnd1)
{
  if (p1 >= pEnd1)
    return T();

  T sum = Square(*p1++ - *p2++);
  for (; p1 < pEnd1; p1++, p2++)
  {
    T difference = *p1 - *p2;
    sum = MultiplyAndAdd(difference, difference, sum);
  }

  return sum;
}

template <class T> MTL_INLINE T Min_Sequential(const T* p, const T* pEnd)
{
  if (p < pEnd)
  {
    T min = *p++;
    for (; p < pEnd; p++)
      min = Min(min, *p);

    return min;
  }
  else
    return T();
}
template <class T> MTL_INLINE T Max_Sequential(const T* p, const T* pEnd)
{
  if (p < pEnd)
  {
    T max = *p++;
    for (; p < pEnd; p++)
      max = Max(max, *p);

    return max;
  }
  else
    return T();
}
template <class T> MTL_INLINE T MinOfAbsolutes_Sequential(const T* p, const T* pEnd)
{
  if (p < pEnd)
  {
    T min = Abs(*p++);
    for (; p < pEnd; p++)
      min = Min(min, Abs(*p));

    return min;
  }
  else
    return T();
}
template <class T> MTL_INLINE T MaxOfAbsolutes_Sequential(const T* p, const T* pEnd)
{
  if (p < pEnd)
  {
    T max = Abs(*p++);
    for (; p < pEnd; p++)
      max = Max(max, Abs(*p));

    return max;
  }
  else
    return T();
}

template <class T> MTL_INLINE
T DotProduct_Sequential(const T* p1, const T* p2, const T* pEnd1)
{
  if (p1 >= pEnd1)
    return T();

  T sum = *p1++ * *p2++;
  for (; p1 < pEnd1; p1++, p2++)
    sum = MultiplyAndAdd(*p1, *p2, sum);

  return sum;
}

template <class T> MTL_INLINE
void AdditionScaled_Sequential(T* pDst, const T* pSrc, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = MultiplyAndAdd(*pSrc, scalar, *pDst);
}


#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
//
// Streamed single-threaded operations.
//
template <class T> MTL_INLINE void
UnaryMinus_StreamAligned_Sequential(T* pDst, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedAligned(pDst);
    xDst = -xDst;
    xDst.StorePackedAligned(pDst);
  }
  UnaryMinus_Sequential(pDst, pDstEnd);
}
template <class T> MTL_INLINE void
UnaryMinus_StreamUnaligned_Sequential(T* pDst, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst(pDst);
    xDst = -xDst;
    xDst.StorePackedUnaligned(pDst);
  }
  UnaryMinus_Sequential(pDst, pDstEnd);
}

template <class T> MTL_INLINE
void Addition_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst += xSrc;
    xDst.StorePackedAligned(pDst);
  }
  Addition_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Addition_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst += xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Addition_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
Subtraction_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst -= xSrc;
    xDst.StorePackedAligned(pDst);
  }
  Subtraction_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Subtraction_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst -= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Subtraction_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
Multiplication_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst *= xSrc;
    xDst.StorePackedAligned(pDst);
  }
  Multiplication_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Multiplication_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst *= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Multiplication_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
Division_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst /= xSrc;
    xDst.StorePackedAligned(pDst);
  }
  Division_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Division_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst /= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Division_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
Min_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst = Min(xDst, xSrc);
    xDst.StorePackedAligned(pDst);
  }
  Min_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Min_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst = Min(xDst, xSrc);
    xDst.StorePackedUnaligned(pDst);
  }
  Min_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
Max_StreamAligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst = Max(xDst, xSrc);
    xDst.StorePackedAligned(pDst);
  }
  Max_Sequential(pDst, pSrc, pDstEnd);
}
template <class T> MTL_INLINE void
Max_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst = Max(xDst, xSrc);
    xDst.StorePackedUnaligned(pDst);
  }
  Max_Sequential(pDst, pSrc, pDstEnd);
}
template <class T, class SrcT> void
Convert_StreamAligned_Sequential(T* pDst, const SrcT* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    Array<XX<T>::Increment,T>::Set(&xDst[0], pSrc);
    xDst.StorePackedAligned(pDst);
  }
  Convert_Sequential(pDst, pSrc, pDstEnd);
}
template <class T, class SrcT> void
Convert_StreamUnaligned_Sequential(T* pDst, const SrcT* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    Array<XX<T>::Increment,T>::Set(&xDst[0], pSrc);
    xDst.StorePackedUnaligned(pDst);
  }
  Convert_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE void
ScalarAddition_StreamAligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedAligned(pDst);
    xDst += xScalar;
    xDst.StorePackedAligned(pDst);
  }
  ScalarAddition_Sequential(pDst, scalar, pDstEnd);
}
template <class T> MTL_INLINE void
ScalarAddition_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst(pDst);
    xDst += xScalar;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarAddition_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE void
ScalarSubtraction_StreamAligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedAligned(pDst);
    xDst -= xScalar;
    xDst.StorePackedAligned(pDst);
  }
  ScalarSubtraction_Sequential(pDst, scalar, pDstEnd);
}
template <class T> MTL_INLINE void
ScalarSubtraction_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst(pDst);
    xDst -= xScalar;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarSubtraction_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE void
ScalarMultiplication_StreamAligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedAligned(pDst);
    xDst *= xScalar;
    xDst.StorePackedAligned(pDst);
  }
  ScalarMultiplication_Sequential(pDst, scalar, pDstEnd);
}
template <class T> MTL_INLINE void
ScalarMultiplication_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst(pDst);
    xDst *= xScalar;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarMultiplication_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE void
ScalarDivision_StreamAligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedAligned(pDst);
    xDst /= xScalar;
    xDst.StorePackedAligned(pDst);
  }
  ScalarDivision_Sequential(pDst, scalar, pDstEnd);
}
template <class T> MTL_INLINE void
ScalarDivision_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst(pDst);
    xDst /= xScalar;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarDivision_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE T
Sum_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx;
    xx.LoadPackedAligned(p);
    xSum += xx;
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + Sum_Sequential(p, pEnd);
}
template <class T> MTL_INLINE void
Sum_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx(p);
    xSum += xx;
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + Sum_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
SumOfAbsolutes_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx;
    xx.LoadPackedAligned(p);
    xSum += Abs(xx);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfAbsolutes_Sequential(p, pEnd);
}
template <class T> MTL_INLINE void
SumOfAbsolutes_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx(p);
    xSum += Abs(xx);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfAbsolutes_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
SumOfSquares_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx;
    xx.LoadPackedAligned(p);
    xSum = MultiplyAndAdd(xx, xx, xSum);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfSquares_Sequential(p, pEnd);
}
template <class T> MTL_INLINE T
SumOfSquares_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xSum(T(0));
  FOR_STREAM(p, size)
  {
    XX<T> xx(p);
    xSum = MultiplyAndAdd(xx, xx, xSum);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfSquares_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
SumOfSquaredDifferences_StreamAligned_Sequential(const T* p1, const T* p2, SizeType size)
{
  const T* pEnd1 = p1 + size;

  XX<T> xSum(T(0));
  FOR_STREAM2(p1, p2, size)
  {
    XX<T> xx1, xx2;
    xx1.LoadPackedAligned(p1);
    xx2.LoadPackedAligned(p2);
    xx1 -= xx2;
    xSum = MultiplyAndAdd(xx1, xx1, xSum);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfSquaredDifferences_Sequential(p1, p2, pEnd1);
}
template <class T> MTL_INLINE T
SumOfSquaredDifferences_StreamUnaligned_Sequential(const T* p1, const T* p2, SizeType size)
{
  const T* pEnd1 = p1 + size;

  XX<T> xSum(T(0));
  FOR_STREAM2(p1, p2, size)
  {
    XX<T> xx1(p1);
    XX<T> xx2(p2);
    xx1 -= xx2;
    xSum = MultiplyAndAdd(xx1, xx1, xSum);
  }
  return Sum< XX<T>::Increment >(xSum.pData()) + SumOfSquaredDifferences_Sequential(p1, p2, pEnd1);
}

template <class T> MTL_INLINE T
Min_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMin;
    xMin.LoadPackedAligned(p);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx;
      xx.LoadPackedAligned(p);
      xMin = Min(xMin, xx);
    }

    if (p == pEnd)
      return Minimum< XX<T>::Increment >(xMin.pData());
    else
      return Min(Minimum< XX<T>::Increment >(xMin.pData()), Min_Sequential(p, pEnd));
  }
  else
    return Min_Sequential(p, pEnd);
}
template <class T> MTL_INLINE T
Min_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMin(p);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx(p);
      xMin = Min(xMin, xx);
    }

    if (p == pEnd)
      return Minimum< XX<T>::Increment >(xMin.pData());
    else
      return Min(Minimum< XX<T>::Increment >(xMin.pData()), Min_Sequential(p, pEnd));
  }
  else
    return Min_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
Max_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMax;
    xMax.LoadPackedAligned(p);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx;
      xx.LoadPackedAligned(p);
      xMax = Max(xMax, xx);
    }

    if (p == pEnd)
      return Maximum< XX<T>::Increment >(xMax.pData());
    else
      return Max(Maximum< XX<T>::Increment >(xMax.pData()), Max_Sequential(p, pEnd));
  }
  else
    return Max_Sequential(p, pEnd);
}
template <class T> MTL_INLINE T
Max_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMax(p);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx(p);
      xMax = Max(xMax, xx);
    }

    if (p == pEnd)
      return Maximum< XX<T>::Increment >(xMax.pData());
    else
      return Max(Maximum< XX<T>::Increment >(xMax.pData()), Max_Sequential(p, pEnd));
  }
  else
    return Max_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
MinOfAbsolutes_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMin;
    xMin.LoadPackedAligned(p);
    xMin = Abs(xMin);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx;
      xx.LoadPackedAligned(p);
      xMin = Min(xMin, Abs(xx));
    }

    if (p == pEnd)
      return Minimum< XX<T>::Increment >(xMin.pData());
    else
      return Min(Minimum< XX<T>::Increment >(xMin.pData()), MinOfAbsolutes_Sequential(p, pEnd));
  }
  else
    return MinOfAbsolutes_Sequential(p, pEnd);
}
template <class T> MTL_INLINE T
MinOfAbsolutes_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMin(p);
    xMin = Abs(xMin);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx(p);
      xMin = Min(xMin, Abs(xx));
    }

    if (p == pEnd)
      return Minimum< XX<T>::Increment >(xMin.pData());
    else
      return Min(Minimum< XX<T>::Increment >(xMin.pData()), MinOfAbsolutes_Sequential(p, pEnd));
  }
  else
    return MinOfAbsolutes_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
MaxOfAbsolutes_StreamAligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMax;
    xMax.LoadPackedAligned(p);
    xMax = Abs(xMax);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx;
      xx.LoadPackedAligned(p);
      xMax = Max(xMax, Abs(xx));
    }

    if (p == pEnd)
      return Maximum< XX<T>::Increment >(xMax.pData());
    else
      return Max(Maximum< XX<T>::Increment >(xMax.pData()), MaxOfAbsolutes_Sequential(p, pEnd));
  }
  else
    return MaxOfAbsolutes_Sequential(p, pEnd);
}
template <class T> MTL_INLINE T
MaxOfAbsolutes_StreamUnaligned_Sequential(const T* p, SizeType size)
{
  const T* pEnd = p + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xMax(p);
    xMax = Abs(xMax);
    p += XX<T>::Increment;

    FOR_STREAM(p, size - XX<T>::Increment)
    {
      XX<T> xx(p);
      xMax = Max(xMax, Abs(xx));
    }

    if (p == pEnd)
      return Maximum< XX<T>::Increment >(xMax.pData());
    else
      return Max(Maximum< XX<T>::Increment >(xMax.pData()), MaxOfAbsolutes_Sequential(p, pEnd));
  }
  else
    return MaxOfAbsolutes_Sequential(p, pEnd);
}

template <class T> MTL_INLINE T
DotProduct_StreamAligned_Sequential(const T* p1, const T* p2, SizeType size)
{
  const T* pEnd1 = p1 + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xDot, x1, x2;
    xDot.LoadPackedAligned(p1);
    x2.LoadPackedAligned(p2);
    xDot *= x2;
    p1 += XX<T>::Increment;
    p2 += XX<T>::Increment;

    FOR_STREAM2(p1, p2, size - XX<T>::Increment)
    {
      x1.LoadPackedAligned(p1);
      x2.LoadPackedAligned(p2);
      xDot = MultiplyAndAdd(x1, x2, xDot);
    }

    return Sum< XX<T>::Increment >(xDot.pData()) + DotProduct_Sequential(p1, p2, pEnd1);
  }
  else
    return DotProduct_Sequential(p1, p2, pEnd1);
}
template <class T> MTL_INLINE T
DotProduct_StreamUnaligned_Sequential(const T* p1, const T* p2, SizeType size)
{
  const T* pEnd1 = p1 + size;

  if (size >= XX<T>::Increment)
  {
    XX<T> xDot(p1);
    xDot *= XX<T>(p2);
    p1 += XX<T>::Increment;
    p2 += XX<T>::Increment;

    FOR_STREAM2(p1, p2, size - XX<T>::Increment)
      xDot = MultiplyAndAdd(XX<T>(p1), XX<T>(p2), xDot);

    return Sum< XX<T>::Increment >(xDot.pData()) + DotProduct_Sequential(p1, p2, pEnd1);
  }
  else
    return DotProduct_Sequential(p1, p2, pEnd1);
}

template <class T> MTL_INLINE
void AdditionScaled_StreamAligned_Sequential(T* pDst, const T* pSrc, const T& scalar,
                                             SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedAligned(pSrc);
    xDst.LoadPackedAligned(pDst);
    xDst = MultiplyAndAdd(xSrc, xScalar, xDst);
    xDst.StorePackedAligned(pDst);
  }
  AdditionScaled_Sequential(pDst, pSrc, scalar, pDstEnd);
}
template <class T> MTL_INLINE
void AdditionScaled_StreamUnaligned_Sequential(T* pDst, const T* pSrc, const T& scalar,
                                               SizeType size)
{
  const T* pDstEnd = pDst + size;

  XX<T> xScalar(scalar);

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc(pSrc);
    XX<T> xDst(pDst);
    xDst = MultiplyAndAdd(xSrc, xScalar, xDst);
    xDst.StorePackedUnaligned(pDst);
  }
  AdditionScaled_Sequential(pDst, pSrc, scalar, pDstEnd);
}

//
// Multi-threaded operations.
//
template <class T> MTL_INLINE void
UnaryMinus_StreamAligned_Parallel(T* pDst, SizeType size)
{
  Parallel_1Dst< T, UnaryMinus_StreamAligned_Sequential<T> >(pDst, size);
}
template <class T> MTL_INLINE void
UnaryMinus_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, UnaryMinus_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Addition_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Addition_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Addition_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Addition_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Subtraction_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Subtraction_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Subtraction_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Subtraction_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Multiplication_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Multiplication_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Multiplication_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Multiplication_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Division_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Division_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Division_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Division_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Min_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Min_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Min_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Min_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
Max_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Max_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE void
Max_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Max_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE void
ScalarAddition_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarAddition_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE void
ScalarAddition_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarAddition_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE void
ScalarSubtraction_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarSubtraction_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE void
ScalarSubtraction_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarSubtraction_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE void
ScalarMultiplication_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarMultiplication_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE void
ScalarMultiplication_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarMultiplication_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE void
ScalarDivision_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarDivision_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE void
ScalarDivision_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarDivision_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE T
Sum_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Sum_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
Sum_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Sum_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResults.Size);
}

template <class T> MTL_INLINE T
SumOfAbsolutes_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, SumOfAbsolutes_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
SumOfAbsolutes_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, SumOfAbsolutes_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
SumOfSquares_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, SumOfSquares_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
SumOfSquares_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, SumOfSquares_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
SumOfSquaredDifferences_StreamAligned_Parallel(const T* p1, const T* p2, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_2Src<T, T, SumOfSquaredDifferences_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p1, p2, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
SumOfSquaredDifferences_StreamUnaligned_Parallel(const T* p1, const T* p2, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_2Src<T, T, SumOfSquaredDifferences_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p1, p2, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
Min_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Min_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Min_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
Min_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Min_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Min_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
Max_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Max_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Max_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
Max_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, Max_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Max_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
MinOfAbsolutes_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, MinOfAbsolutes_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Min_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
MinOfAbsolutes_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, MinOfAbsolutes_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Min_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
MaxOfAbsolutes_StreamAligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
    ParallelReduction_1Src<T, T, MaxOfAbsolutes_StreamAligned_Sequential<T> >
      (subResults, subResultsSize, p, size);

  return Max_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
MaxOfAbsolutes_StreamUnaligned_Parallel(const T* p, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_1Src<T, T, MaxOfAbsolutes_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p, size);

  return Max_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE T
DotProduct_StreamAligned_Parallel(const T* p1, const T* p2, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_2Src<T, T, DotProduct_StreamAligned_Sequential<T> >
    (subResults, subResultsSize, p1, p2, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}
template <class T> MTL_INLINE T
DotProduct_StreamUnaligned_Parallel(const T* p1, const T* p2, SizeType size)
{
  T subResults[MTL_MAX_THREADS];
  SizeType subResultsSize;
  ParallelReduction_2Src<T, T, DotProduct_StreamUnaligned_Sequential<T> >
    (subResults, subResultsSize, p1, p2, size);

  return Sum_StreamAligned_Sequential(subResults, subResultsSize);
}

template <class T> MTL_INLINE void
AdditionScaled_StreamAligned_Parallel(T* pDst, const T* pSrc, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Src_1Val< T, AdditionScaled_StreamAligned_Sequential<T> >(pDst, pSrc,
                                                                           scalar, size);
}
template <class T> MTL_INLINE void
AdditionScaled_StreamUnaligned_Parallel(T* pDst, const T* pSrc, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Src_1Val< T, AdditionScaled_StreamUnaligned_Sequential<T> >(pDst, pSrc,
                                                                             scalar, size);
}
#endif  // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX

}  // namespace MTL

#endif  // MTL_STREAM_ARRAY_H

#endif  // #if MTL_ENABLE_SSE || MTL_ENABLE_AVX
