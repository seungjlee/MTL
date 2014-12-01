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

#ifndef MTL_OPENMP_H
#define MTL_OPENMP_H

#include "CPU.h"

// Note that I have not really tested the code with OpenMP disabled.
#ifndef MTL_ENABLE_OPENMP
  #define MTL_ENABLE_OPENMP 1
#endif

#include <omp.h>

#ifndef MIN_OPENMP_DATA_SIZE
  #define MIN_OPENMP_DATA_SIZE 4096  // Minimum data size per thread for OpenMP (in bytes).
#endif

namespace MTL
{

template <class T> class DynamicVector;

// Returns true if OpenMP should be used for the data size and number of threads.
template <class T>
MTL_INLINE static bool DoOpenMP(SizeType size, SizeType numberOfThreads)
{
  return numberOfThreads > 1 && size >= numberOfThreads * MIN_OPENMP_DATA_SIZE / sizeof(T);
}

//
// Helper classes for multithreading.
//
template <class T, void (*Func)(T*, SizeType)>
MTL_INLINE static void Parallel_1Dst(T* p, SizeType size)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(size, numberOfThreads))
  {
    DynamicVector<SizeType> subSizes, offsets;
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(p + offsets[i], subSizes[i]);
  }
  else
#endif
    Func(p, size);
}

template <class T, void (*Func)(T*, const T*, SizeType)>
MTL_INLINE static void Parallel_1Dst_1Src(T* pDst, const T* pSrc, SizeType size)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(size, numberOfThreads))
  {
    DynamicVector<SizeType> subSizes, offsets;
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(pDst + offsets[i], pSrc + offsets[i], subSizes[i]);
  }
  else
#endif
    Func(pDst, pSrc, size);
}

template <class T, void (*Func)(T*, const T&, SizeType)>
MTL_INLINE static void Parallel_1Dst_1Val(T* p, const T& val, SizeType size)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(size, numberOfThreads))
  {
    DynamicVector<SizeType> subSizes, offsets;
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(p + offsets[i], val, subSizes[i]);
  }
  else
#endif
    Func(p, val, size);
}

template <class ReductionT, class T, ReductionT (*Func)(const T*, SizeType)>
MTL_INLINE static DynamicVector<ReductionT> ParallelReduction_1Src(const T* pSrc, SizeType size)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(size, numberOfThreads))
  {
    DynamicVector<ReductionT> subResults(numberOfThreads);

    DynamicVector<SizeType> subSizes, offsets;
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      subResults[i] = Func(pSrc + offsets[i], subSizes[i]);

    return subResults;
  }
  else
#endif
    return DynamicVector<ReductionT>(1, Func(pSrc, size));
}

template <class ReductionT, class T, ReductionT (*Func)(const T*, const T*, SizeType)>
MTL_INLINE static DynamicVector<ReductionT> ParallelReduction_2Src(const T* pSrc1, const T* pSrc2,
                                                                   SizeType size)
{
#if MTL_ENABLE_OPENMP
  I64 numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (DoOpenMP<T>(size, numberOfThreads))
  {
    DynamicVector<ReductionT> subResults(numberOfThreads);

    DynamicVector<SizeType> subSizes, offsets;
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (long i = 0; i < numberOfThreads; i++)
      subResults[i] = Func(pSrc1 + offsets[i], pSrc2 + offsets[i], subSizes[i]);

    return subResults;
  }
  else
#endif
    return DynamicVector<ReductionT>(1, Func(pSrc1, pSrc2, size));
}

}  // namespace MTL


#endif  // MTL_OPENMP_H
