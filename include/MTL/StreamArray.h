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


#ifndef MTL_STREAM_ARRAY_H
#define MTL_STREAM_ARRAY_H

#include "SSE.h"
#include "AVX.h"

namespace MTL
{

template <class T> MTL_INLINE static void
UnaryMinus_Sequential(T* pDst, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst = -*pDst;
}
template <class T> MTL_INLINE static void
ScalarAddition_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst += scalar;
}
template <class T> MTL_INLINE static void
ScalarSubtraction_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst -= scalar;
}
template <class T> MTL_INLINE static void
ScalarMultiplication_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst *= scalar;
}
template <class T> MTL_INLINE static void
ScalarDivision_Sequential(T* pDst, const T& scalar, const T* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++)
    *pDst /= scalar;
}
template <class DstT, class SrcT> MTL_INLINE static void
Addition_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst + *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE static void
Subtraction_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst - *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE static void
Multiplication_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst * *pSrc);
}
template <class DstT, class SrcT> MTL_INLINE static void
Division_Sequential(DstT* pDst, const SrcT* pSrc, const DstT* pDstEnd)
{
  for (; pDst < pDstEnd; pDst++, pSrc++)
    *pDst = DstT(*pDst / *pSrc);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
UnaryMinus_StreamUnaligned_Sequential(T* pDst, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM(pDst, size)
  {
    XX<T> xDst;
    xDst.LoadPackedUnaligned(pDst);
    xDst = -xDst;
    xDst.StorePackedUnaligned(pDst);
  }
  UnaryMinus_Sequential(pDst, pDstEnd);
}

template <class T> MTL_INLINE static
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
template <class T> MTL_INLINE static void
Addition_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedUnaligned(pSrc);
    xDst.LoadPackedUnaligned(pDst);
    xDst += xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Addition_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
Subtraction_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedUnaligned(pSrc);
    xDst.LoadPackedUnaligned(pDst);
    xDst -= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Subtraction_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
Multiplication_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedUnaligned(pSrc);
    xDst.LoadPackedUnaligned(pDst);
    xDst *= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Multiplication_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
Division_StreamUnaligned_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xSrc, xDst;
    xSrc.LoadPackedUnaligned(pSrc);
    xDst.LoadPackedUnaligned(pDst);
    xDst /= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  Division_Sequential(pDst, pSrc, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
ScalarAddition_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    xDst.LoadPackedUnaligned(pDst);
    xDst += xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarAddition_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
ScalarSubtraction_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    xDst.LoadPackedUnaligned(pDst);
    xDst -= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarSubtraction_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
ScalarMultiplication_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    xDst.LoadPackedUnaligned(pDst);
    xDst *= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarMultiplication_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE static void
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
template <class T> MTL_INLINE static void
ScalarDivision_StreamUnaligned_Sequential(T* pDst, const T& scalar, SizeType size)
{
  const T* pDstEnd = pDst + size;

  FOR_STREAM2(pDst, pSrc, size)
  {
    XX<T> xDst;
    xDst.LoadPackedUnaligned(pDst);
    xDst /= xSrc;
    xDst.StorePackedUnaligned(pDst);
  }
  ScalarDivision_Sequential(pDst, scalar, pDstEnd);
}

template <class T> MTL_INLINE static void
UnaryMinus_StreamAligned_Parallel(T* pDst, SizeType size)
{
  Parallel_1Dst< T, UnaryMinus_StreamAligned_Sequential<T> >(pDst, size);
}
template <class T> MTL_INLINE static void
UnaryMinus_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, UnaryMinus_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE static void
Addition_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Addition_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE static void
Addition_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Addition_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE static void
Subtraction_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Subtraction_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE static void
Subtraction_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Subtraction_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE static void
Multiplication_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Multiplication_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE static void
Multiplication_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Multiplication_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE static void
Division_StreamAligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Division_StreamAligned_Sequential<T> >(pDst, pSrc, size);
}
template <class T> MTL_INLINE static void
Division_StreamUnaligned_Parallel(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, Division_StreamUnaligned_Sequential<T> >(pDst, pSrc, size);
}

template <class T> MTL_INLINE static void
ScalarAddition_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarAddition_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE static void
ScalarAddition_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarAddition_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE static void
ScalarSubtraction_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarSubtraction_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE static void
ScalarSubtraction_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarSubtraction_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE static void
ScalarMultiplication_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarMultiplication_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE static void
ScalarMultiplication_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarMultiplication_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

template <class T> MTL_INLINE static void
ScalarDivision_StreamAligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarDivision_StreamAligned_Sequential<T> >(pDst, scalar, size);
}
template <class T> MTL_INLINE static void
ScalarDivision_StreamUnaligned_Parallel(T* pDst, const T& scalar, SizeType size)
{
  Parallel_1Dst_1Val< T, ScalarDivision_StreamUnaligned_Sequential<T> >(pDst, scalar, size);
}

}  // namespace MTL

#endif  // MTL_STREAM_ARRAY_H
