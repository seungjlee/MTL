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

#ifndef MTL_OPENMP_H
#define MTL_OPENMP_H

#include <MTL/CPU.h>
#include <MTL/SSE.h>
#include <MTL/AVX.h>
#include <assert.h>

//
// Microsoft Visual Studio's OpenMP implementation by default sets the environment variable
// OMP_WAIT_POLICY to ACTIVE. I think they changed this behaviour from VS 2005 to VS 2008.
// They have claimed that this is better for performance but through the years but in most regular
// systems it causes a penalty by wasting CPU cycles waiting for the next OpenMP threaded
// task. I do recommend setting/creating the environment variable OMP_WAIT_POLICY with value
// PASSIVE for those using VS 2008 and newer unless you have a system with relatively large number
// of cores. Anyway, it should be tested fully before you decide which setting of OMP_WAIT_POLICY
// is better for your system.
//

#ifndef MTL_MAX_THREADS
  #define MTL_MAX_THREADS 1024
#endif

#include <omp.h>

#ifndef MIN_OPENMP_DATA_SIZE
  #define MIN_OPENMP_DATA_SIZE 4096  // Minimum data size per thread for OpenMP (in bytes).
#endif

namespace MTL
{

// Returns true if OpenMP should be used for the data size and number of threads.
template <class T>
MTL_INLINE static bool DoOpenMP(SizeType size, SizeType numberOfThreads)
{
  return numberOfThreads > 1 && size >= numberOfThreads * MIN_OPENMP_DATA_SIZE / sizeof(T);
}

//
// Helper fuctions for multithreading.
//
MTL_INLINE static SizeType ComputeParallelSubSizesBlockSize(SizeType size, U64 numberOfThreads)
{
  SizeType blockSize = size / SizeType(numberOfThreads);

  if (blockSize * numberOfThreads < size)
    blockSize++;

  return blockSize;
}

// Computes sizes and offsets for parallel processing.
template <class T> MTL_INLINE static void ComputeParallelSubSizes
(SizeType* subSizes, SizeType* offsets,
 SizeType totalSize, U64 numberOfThreads)
{
  assert(numberOfThreads <= MTL_MAX_THREADS);

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
  SizeType chunkSize = MTL::XX<T>::StreamSize(totalSize / SizeType(numberOfThreads));
#else
  SizeType chunkSize = totalSize / SizeType(numberOfThreads);
#endif

  for (U64 i = 0; i < numberOfThreads; i++)
  {
    offsets[i] = 0;
    subSizes[i] = chunkSize;
  }

  // Doing the laziest thing for now. Add the leftover to the last chunk.
  SizeType remainder = totalSize - chunkSize * SizeType(numberOfThreads);
  subSizes[numberOfThreads - 1] += remainder;

  // Compute offsets.
  for (U64 i = 1; i < numberOfThreads; i++)
    offsets[i] = offsets[i-1] + subSizes[i-1];
}
MTL_INLINE static void ComputeParallelSubHeights
(SizeType* subHeights, SizeType* offsets,
 SizeType height, U64 numberOfThreads, SizeType overlap = 0)
{
  assert(numberOfThreads <= MTL_MAX_THREADS);

  SizeType effectiveHeight = height - overlap;
  SizeType chunkHeight = effectiveHeight / SizeType(numberOfThreads);

  for (U64 i = 0; i < numberOfThreads; i++)
  {
    offsets[i] = 0;
    subHeights[i] = chunkHeight + overlap;
  }

  // Distribute the remainder.
  SizeType remainder = effectiveHeight - chunkHeight * SizeType(numberOfThreads);
  for (SizeType i = 0; i < remainder; i++)
    subHeights[i]++;

  // Compute offsets.
  for (U64 i = 1; i < numberOfThreads; i++)
    offsets[i] = offsets[i-1] + subHeights[i-1] - overlap;
}

template <class T, void (*Func)(T*, SizeType)>
MTL_INLINE static void Parallel_1Dst(T* p, SizeType size, I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<T>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
    {
      SizeType index = i;
      Func(p + offsets[index], subSizes[index]);
    }
  }
  else
#endif
    Func(p, size);
}

template <class DstT, class SrcT, void (*Func)(DstT*, const SrcT*, SizeType)>
MTL_INLINE static void Parallel_1Dst_1Src(DstT* pDst, const SrcT* pSrc, SizeType size,
                                          I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<DstT>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<DstT>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(pDst + offsets[i], pSrc + offsets[i], subSizes[i]);
  }
  else
#endif
    Func(pDst, pSrc, size);
}

template <class T, void (*Func)(T*, const T&, SizeType)>
MTL_INLINE static void Parallel_1Dst_1Val(T* p, const T& val, SizeType size, I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<T>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(p + offsets[i], val, subSizes[i]);
  }
  else
#endif
    Func(p, val, size);
}

template <class T, void (*Func)(T*, const T*, const T&, SizeType)>
MTL_INLINE static void Parallel_1Dst_1Src_1Val(T* pDst, const T* pSrc, const T& val, SizeType size,
                                               I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<T>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      Func(pDst + offsets[i], pSrc + offsets[i], val, subSizes[i]);
  }
  else
#endif
    Func(pDst, pSrc, val, size);
}

template <class ReductionT, class T, ReductionT (*Func)(const T*, SizeType)>
MTL_INLINE static void ParallelReduction_1Src(ReductionT* subResults, SizeType& subResultsSize,
                                              const T* pSrc, SizeType size, I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<T>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      subResults[i] = Func(pSrc + offsets[i], subSizes[i]);

    subResultsSize = SizeType(numberOfThreads);
  }
  else
#endif
  {
    subResults[0] = Func(pSrc, size);
    subResultsSize = 1;
  }
}

template <class ReductionT, class T, ReductionT (*Func)(const T*, const T*, SizeType)>
MTL_INLINE static void ParallelReduction_2Src(ReductionT* subResults, SizeType& subResultsSize,
                                              const T* pSrc1, const T* pSrc2, SizeType size,
                                              I64 numberOfThreads)
{
#if MTL_ENABLE_OPENMP
  if (numberOfThreads < 1)
    numberOfThreads = MTL::CPU::Instance().NumberOfThreads();
  if (numberOfThreads > MTL_MAX_THREADS)
    numberOfThreads = MTL_MAX_THREADS;

  if (DoOpenMP<T>(size, SizeType(numberOfThreads)))
  {
    SizeType subSizes[MTL_MAX_THREADS], offsets[MTL_MAX_THREADS];
    ComputeParallelSubSizes<T>(subSizes, offsets, size, numberOfThreads);

    #pragma omp parallel for
    for (I32 i = 0; i < numberOfThreads; i++)
      subResults[i] = Func(pSrc1 + offsets[i], pSrc2 + offsets[i], subSizes[i]);

    subResultsSize = SizeType(numberOfThreads);
  }
  else
#endif
  {
    subResults[0] = Func(pSrc1, pSrc2, size);
    subResultsSize = 1;
  }
}

}  // namespace MTL


#if MTL_ENABLE_OPENMP
#ifdef WIN32
#define MTL_PARALLEL_FOR_BLOCKS(Size)                                                       \
  int Number_Of_Threads = (int)MTL::CPU::Instance().NumberOfThreads();                      \
  int Block_Size = (int)MTL::ComputeParallelSubSizesBlockSize(Size, Number_Of_Threads);     \
  __pragma(omp parallel for num_threads(Number_Of_Threads) schedule(dynamic, Block_Size))
#define MTL_PARALLEL_FOR_BLOCKS_THREADS(Size, NumberOfThreads)                              \
  int Block_Size = (int)MTL::ComputeParallelSubSizesBlockSize(Size, NumberOfThreads);       \
  __pragma(omp parallel for num_threads(NumberOfThreads) schedule(dynamic, Block_Size))
#else
#define MTL_PRAGMA(x) _Pragma(#x)
#define MTL_PARALLEL_FOR_BLOCKS(Size)                                                       \
  int Number_Of_Threads = (int)MTL::CPU::Instance().NumberOfThreads();                      \
  int Block_Size = (int)MTL::ComputeParallelSubSizesBlockSize(Size, Number_Of_Threads);     \
  MTL_PRAGMA(omp parallel for num_threads(Number_Of_Threads) schedule(dynamic, Block_Size))
#define MTL_PARALLEL_FOR_BLOCKS_THREADS(Size, NumberOfThreads)                              \
  int Number_Of_Threads = (int)NumberOfThreads;                                             \
  int Block_Size = (int)MTL::ComputeParallelSubSizesBlockSize(Size, NumberOfThreads);       \
  MTL_PRAGMA(omp parallel for num_threads(Number_Of_Threads) schedule(dynamic, Block_Size))
#endif
#else  // MTL_ENABLE_OPENMP
#define MTL_PARALLEL_FOR_BLOCKS(Size)
#define MTL_PARALLEL_FOR_BLOCKS_THREADS(Size, NumberOfThreads)
#endif // MTL_ENABLE_OPENMP

#endif  // MTL_OPENMP_H
