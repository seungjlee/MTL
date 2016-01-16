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


#ifndef MTL_DYNAMIC_VECTOR_H
#define MTL_DYNAMIC_VECTOR_H

#include <assert.h>
#include <omp.h>
#include <MTL/Matrix.h>
#include <MTL/StreamMath.h>

//
// Some user overridable macros/constants.
//

// Initial allocation block.
#ifndef MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK
  #define MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK 1
#endif

// A limit to the maximum additional block allocated to avoid bad bugs.
#ifndef MTL_MAX_ALLOCATION_BLOCK
  #define MTL_MAX_ALLOCATION_BLOCK 512 * 1024 * 1024
#endif

// Block size for optimizing AssignAll() operations. 
#ifndef MTL_ASSIGN_ALL_BLOCK_SIZE
  #define MTL_ASSIGN_ALL_BLOCK_SIZE 1024
#endif

// Loop helpers.
#define FOR_EACH(Index, Initial, Size)  for(MTL::SizeType Index = Initial; Index < Size; Index++)
#define FOR_EACH_(V, Index)             FOR_EACH(Index, 0, V.Size())
#define FOR_EACH_INDEX(V)               FOR_EACH_(V, V##Index)


namespace MTL
{

template <class T>
class DynamicVector
{
public:
  // Default constructor.
  MTL_INLINE DynamicVector()
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK) {}

  // Constructs array of specified size. Elements might be uninitialized.
  MTL_INLINE DynamicVector(SizeType size)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    Resize(size);
  }

  // Constructs array of specified size. Elements are initialized with the initialValue.
  MTL_INLINE DynamicVector(SizeType size, const T& initialValue)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    Resize(size, initialValue);
  }

  // Copy constructor.
  MTL_INLINE DynamicVector(const DynamicVector& rhs)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    *this = rhs;
  }

  // Constructors that initialize from array.
  MTL_INLINE DynamicVector(const T* pSrc, SizeType size)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    Assign(pSrc, size);
  }
  MTL_INLINE DynamicVector(const T* pSrcBegin, const T* pSrcEnd)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    Assign(pSrcBegin, pSrcEnd);
  }

  // Constructor that casts types.
  template <class T2>
  MTL_INLINE explicit DynamicVector(const DynamicVector<T2>& rhs)
    : Buffer_(0), First_(0), Size_(0), BufferSize_(0),
      AllocBlock_(MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK)
  {
    Resize(rhs.Size());

    SizeType size = Size();
    T* pDst = First_;
    const T2* pSrc = rhs.Begin();

#if MTL_ENABLE_OPENMP
    if (DoOpenMP<T>(size, MTL::CPU::Instance().NumberOfThreads()))
    {
      #pragma omp parallel for
      for (long k = 0; k < (long)size; k++)
        pDst[k] = T(pSrc[k]);
    }
    else
#endif
    {
      const T* pDstEnd = pDst + size;

      for (; pDst < pDstEnd; pDst++, pSrc++)
        *pDst = T(*pSrc);
    }
  }

  ~DynamicVector()  { DestroyAndDelete(); }

  // Resizes vector to the new number of elements.
  MTL_INLINE void Resize(SizeType newSize) 
  {
    Reserve(newSize);
    Size_ = newSize;
  }

  // Resizes vector to the new number of elements and assigns all elements to the new value.
  MTL_INLINE void Resize(SizeType newSize, const T& value) 
  {
    if (newSize > BufferSize_)
    {
      SizeType newBufferSize;
      T* pBuffer;
      T* pFirst = AllocateAlignedMemory(newBufferSize, &pBuffer, newSize);
      assert((reinterpret_cast<U64>(pFirst) & (MTL_STREAM_BYTES-1)) == 0);

      ConstructElements(pFirst, pFirst + newBufferSize, value);

      // Destroy old elements.
      DestroyAndDelete();

      Buffer_ = pBuffer;
      First_ = pFirst;
      BufferSize_ = newBufferSize;

      AllocBlock_ = Min(BufferSize_, (SizeType)MTL_MAX_ALLOCATION_BLOCK);
      Size_ = newSize;
    }
    else
    {
      Size_ = newSize;
      AssignAll(First_, First_ + Size_, value);
    }
  }

  MTL_INLINE void Reserve(SizeType newSize)
  {
    if (newSize > BufferSize_)
    {
      SizeType newBufferSize;
      T* pBuffer;
      T* pFirst = AllocateAlignedMemory(newBufferSize, &pBuffer, newSize);
      assert((reinterpret_cast<SizeType>(pFirst) & (MTL_STREAM_BYTES-1)) == 0);

      // Call default constructors for new elements. This can be disabled with the macro
      // MTL_DYNAMIC_VECTOR_NO_CONSTRUCTOR_DESTRUCTOR.
      ConstructElements(pFirst, pFirst + newBufferSize);

      // User can enable this feature with the macro MTL_DYNAMIC_VECTOR_ZERO_INIT.
      FastZeroInit(pFirst, newBufferSize);

      OptimizedCopy(pFirst, First_, Size_);  // Assuming no overlap in buffers.

      // Destroy old elements.
      DestroyAndDelete();

      Buffer_ = pBuffer;
      First_ = pFirst;
      BufferSize_ = newBufferSize;

      AllocBlock_ = Min(BufferSize_, (SizeType)MTL_MAX_ALLOCATION_BLOCK);
    }
  }

  MTL_INLINE void PushBack(const T& newElement) 
  {
    if (Size_ == BufferSize_)
    {
      T temp = newElement;  // Copy in case the new element is in this vector.
      Resize(Size_ + 1);
      First_[Size_ - 1] = temp;
    }
    else
    {
      Size_++;
      First_[Size_ - 1] = newElement;
    }
  }

  MTL_INLINE void AddBack(const DynamicVector& newElements)
  {
    AddBack(newElements.Begin(), newElements.Size());
  }
  MTL_INLINE void AddBack(const T* pSrc, SizeType sourceSize)
  {
    Insert(End(), pSrc, sourceSize);
  }

  MTL_INLINE void Insert(T* pDst, const T* pSrc, SizeType sourceSize)
  {
    assert(pDst >= First_ && pDst <= End());
    assert(pSrc >= End() || pSrc + sourceSize < First_);

    SizeType moveSize = End() - pDst;
    SizeType insertOffset = pDst - First_;
    SizeType newTotalSize = Size_ + sourceSize;
    Reserve(newTotalSize);
    T* pNewDst = First_ + insertOffset;
    CopyBackwards(pNewDst + sourceSize, pNewDst, moveSize);
    OptimizedCopy(pNewDst, pSrc, sourceSize);
    Size_ = newTotalSize;
  }
  MTL_INLINE void Insert(T* pDst, const T* pSrc, const T* pSrcEnd)
  {
    assert(pSrcEnd >= pSrc);
    Insert(pDst, pSrc, pSrcEnd - pSrc);
  }
  MTL_INLINE void Insert(T* pDst, const T& element)
  {
    Insert(pDst, &element, 1);
  }

  MTL_INLINE static void CopyBackwards(T* pDst, const T* pSrc, SizeType size)
  {
    const T* pDstBegin = pDst;
    for (pDst += size-1, pSrc += size-1; pDst >= pDstBegin; pDst--, pSrc--)
      *pDst = *pSrc;
  }

  MTL_INLINE SizeType Size() const      { return Size_;       }
  MTL_INLINE SizeType Capacity() const  { return BufferSize_; }

  MTL_INLINE DynamicVector& operator=(const DynamicVector& rhs) 
  {
    if (this != &rhs)
    {
      Resize(rhs.Size_);
      OptimizedCopy(First_, rhs.First_, Size_);  // Assuming no overlap in buffers.
    }

    return *this;
  }

  MTL_INLINE bool operator==(const DynamicVector& rhs) const
  {
    return Size() == rhs.Size() && !memcmp(Begin(), rhs.Begin(), Size()*sizeof(T));
  }

  MTL_INLINE T& operator[](SizeType i)
  {
    assert(i >= 0 && i < Size_);
    return First_[i];
  }
  MTL_INLINE const T& operator[](SizeType i) const
  {
    assert(i >= 0 && i < Size_);
    return First_[i];
  }

  // Sets all vector values to zero.
  MTL_INLINE void Zeros();

  // Sets all values to input val.
  MTL_INLINE void SetAll(const T& val)
  {
    AssignAll(Begin(), End(), val);
  }

  MTL_INLINE T& Back()               { return *(First_ + Size_ - 1); }
  MTL_INLINE const T& Back() const   { return *(First_ + Size_ - 1); }

  MTL_INLINE T* Begin()              { return First_;                }
  MTL_INLINE const T* Begin() const  { return First_;                }
  MTL_INLINE T* End()                { return First_ + Size_;        }
  MTL_INLINE const T* End() const    { return First_ + Size_;        }

  MTL_INLINE T* begin()              { return First_;                }
  MTL_INLINE const T* begin() const  { return First_;                }
  MTL_INLINE T* end()                { return First_ + Size_;        }
  MTL_INLINE const T* end() const    { return First_ + Size_;        }

  // Note that Clear() does not free resources. Call Release() if you really want to free up
  // resources.
  MTL_INLINE void Clear()  { Size_ = 0; }
  MTL_INLINE void Release()
  {
    DestroyAndDelete();
    Buffer_ = 0;
    First_ = 0;
    BufferSize_ = 0;
    AllocBlock_ = MTL_DEFAULT_INITIAL_ALLOCATION_BLOCK;
    Size_ = 0;
  }

  MTL_INLINE void Assign(const T* pSrc, SizeType size)
  {
    assert(pSrc + size <= Begin() || pSrc >= End());
    Resize(size);
    OptimizedCopy(First_, pSrc, Size_);
  }
  MTL_INLINE void Assign(const T* pSrcBegin, const T* pSrcEnd)
  {
    Assign(pSrcBegin, pSrcEnd - pSrcBegin);
  }

protected:
  MTL_INLINE static void OptimizedCopy(T* pDst, const T* pSrc, SizeType size)
  {
    Copy(pDst, pSrc, size);
  }

  MTL_INLINE static void Copy(T* pDst, const T* pSrc, SizeType size)
  {
    for (SizeType i = 0; i < size; i++)
      pDst[i] = pSrc[i];
  }

  MTL_INLINE static void AssignAll(T* p, const T* pEnd, const T& val)
  {
#if MTL_ENABLE_OPENMP
    SizeType size = SizeType(pEnd - p);
    if (DoOpenMP<T>(size, MTL::CPU::Instance().NumberOfThreads()))
    {
      #pragma omp parallel for
      for (long k = 0; k < (long)size; k++)
        p[k] = val;
    }
    else
#endif
    {
      for (; p < pEnd; p++)
        *p = val;
    }
  }

private:
  T *Buffer_;
  T *First_;
  SizeType Size_;        // Vector size in number of elements.
  SizeType BufferSize_;  // Actual buffer size for this vector.
  SizeType AllocBlock_;  // Allocation block size.

  MTL_INLINE T* AllocateAlignedMemory(SizeType& newBufferSize, T** pBuffer, SizeType newSize)
  {
    // Might need to align the pointer to first element here?
    newBufferSize = (newSize / AllocBlock_ + 1) * AllocBlock_;
    *pBuffer = (T*) ::operator new (sizeof(T) * newBufferSize + MTL_STREAM_BYTES-1);
    SizeType address = reinterpret_cast<SizeType>(*pBuffer);
    SizeType offset = address & (MTL_STREAM_BYTES-1);

	T* pFirst = *pBuffer;
	if (offset != 0)
    {
      SizeType complement = MTL_STREAM_BYTES - offset;
      pFirst = reinterpret_cast<T*>(address + complement);
    }
    return pFirst;
  }

  MTL_INLINE void DestroyAndDelete()
  {
    DestroyElements(First_, First_ + BufferSize_);
    ::operator delete(Buffer_);
  }

  MTL_INLINE static void ConstructElements(T *p, const T* pEnd)
  {
    for (; p < pEnd; p++)
      ::new (p) T;
  }
  MTL_INLINE static void ConstructElements(T *p, const T* pEnd, const T& value)
  {
    for (; p < pEnd; p++)
      ::new (p) T(value);
  }

  MTL_INLINE static void DestroyElements(T *p, const T* pEnd)
  {
    for (; p < pEnd; p++)
      p->~T();
  }

  // User can override this to actually zero elements during initialization by defining the macro
  // MTL_DYNAMIC_VECTOR_ZERO_INIT.
  MTL_INLINE static void FastZeroInit(T* p, SizeType size) {}
};

}  // namespace MTL

#include <MTL/OpenMP.h>
#include <MTL/StreamArray.h>

namespace MTL
{

#ifdef WIN32
template <class T> MTL_INLINE static void OptimizedCopy(T* pDst, const T* pSrc, SizeType size);
template <class T> MTL_INLINE static void OptimizedZeros(T* p, SizeType size);
#else
template <class T> inline void OptimizedCopy(T* pDst, const T* pSrc, SizeType size);
template <class T> inline void OptimizedZeros(T* p, SizeType size);
#endif

template <class T> MTL_INLINE void MTL::DynamicVector<T>::Zeros()
{
  OptimizedZeros(First_, Size_);
}

template <class T> MTL_INLINE static void OptimizedZeros_Sequential(T* p, SizeType size)
{
  memset(p, 0, size * sizeof(T));
}
template <class T> MTL_INLINE static void OptimizedZeros(T* p, SizeType size)
{
  Parallel_1Dst< T, OptimizedZeros_Sequential >(p, size);
}

template <class T>
MTL_INLINE static void OptimizedCopy_Sequential(T* pDst, const T* pSrc, SizeType size)
{
  memcpy(pDst, pSrc, size * sizeof(T));
}
template <class T> MTL_INLINE static void OptimizedCopy(T* pDst, const T* pSrc, SizeType size)
{
  Parallel_1Dst_1Src< T, OptimizedCopy_Sequential >(pDst, pSrc, size);
}

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
template <class T>
MTL_INLINE static void AssignAll_Stream_Unaligned_Sequential(T* p, const T& val, SizeType size)
{
  const T* pEnd = p + size;

  XX<T> xVal = XX_SetPacked(val);

  FOR_STREAM(p, size)
    xVal.StorePackedUnaligned(p);

  for (; p < pEnd; p++)
    *p = val;
}

template <class T>
MTL_INLINE static void AssignAll_Stream(T* p, const T& val, SizeType size)
{
  SizeType byteSize = size * sizeof(T);

  if (byteSize >= MTL_ASSIGN_ALL_BLOCK_SIZE)
  {
    XX<T> xVal(val);

    SizeType blockCopySize = byteSize & ~(MTL_ASSIGN_ALL_BLOCK_SIZE-1);

    U8* pDst = (U8*)p;
    const U8* pBlocksEnd = pDst + blockCopySize;
    const U8* pFirstBlockEnd = pDst + MTL_ASSIGN_ALL_BLOCK_SIZE;

    for (; pDst < pFirstBlockEnd; pDst += 16)
      memcpy(pDst, &xVal, sizeof(__m128));

    for (; pDst < pBlocksEnd; pDst += MTL_ASSIGN_ALL_BLOCK_SIZE)
      memcpy(pDst, pDst - MTL_ASSIGN_ALL_BLOCK_SIZE, MTL_ASSIGN_ALL_BLOCK_SIZE);

    memcpy(pDst, pDst - MTL_ASSIGN_ALL_BLOCK_SIZE, byteSize - blockCopySize);
  }
  else
  {
    AssignAll_Stream_Unaligned_Sequential(p, val, size);
  }
}
#endif

template <class T> MTL_INLINE static void OptimizedAssignAll(T* p, const T& val, SizeType size)
{
#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
  if (MTL::CPU::Instance().NumberOfThreads() > 1)
    Parallel_1Dst_1Val< T, AssignAll_Stream<T> >(p, val, size);
  else
    AssignAll_Stream<T>(p, val, size);
#else
  #if MTL_ENABLE_OPENMP
    if (DoOpenMP<T>(size, MTL::CPU::Instance().NumberOfThreads()))
    {
      #pragma omp parallel for
      for (long k = 0; k < (long)size; k++)
        p[k] = val;
    }
    else
  #endif
    {
      const T* pEnd = p + size;
      for (; p < pEnd; p++)
        *p = val;
    }
#endif
}

}  // namespace MTL

#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ASSIGN_ALL(T)                                                \
namespace MTL                                                                                     \
{                                                                                                 \
template <> MTL_INLINE void MTL::DynamicVector<T>::AssignAll(T* p, const T* pEnd, const T& val)   \
{                                                                                                 \
  assert(p <= pEnd);                                                                              \
  MTL::SizeType size = MTL::SizeType(pEnd - p);                                                   \
  MTL::OptimizedAssignAll(p, val, size);                                                          \
}                                                                                                 \
}

#if MTL_ENABLE_SSE || MTL_ENABLE_AVX
#ifdef WIN32
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL(T)                                      \
namespace MTL                                                                                     \
{                                                                                                 \
template <> MTL_INLINE static void OptimizedZeros_Sequential(T* p, SizeType size)                 \
{                                                                                                 \
  AssignAll_Stream<T>(p, T(0), size);                                                             \
}                                                                                                 \
template <> MTL_INLINE static void OptimizedZeros(T* p, SizeType size)                            \
{                                                                                                 \
  Parallel_1Dst< T, OptimizedZeros_Sequential<T> >(p, size);                                      \
}                                                                                                 \
}
#else
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL(T)                                      \
namespace MTL                                                                                     \
{                                                                                                 \
template <> inline void OptimizedZeros_Sequential(T* p, SizeType size)                            \
{                                                                                                 \
  AssignAll_Stream<T>(p, T(0), size);                                                             \
}                                                                                                 \
template <> inline void OptimizedZeros(T* p, SizeType size)                                       \
{                                                                                                 \
  Parallel_1Dst< T, OptimizedZeros_Sequential<T> >(p, size);                                      \
}                                                                                                 \
}
#endif
#else
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL(T)
#endif

#ifdef WIN32
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL_CAST(T)                                 \
namespace MTL                                                                                     \
{                                                                                                 \
template <> MTL_INLINE static void OptimizedZeros(T* p, SizeType size)                            \
{                                                                                                 \
  Parallel_1Dst< I8, OptimizedZeros_Sequential<I8> >((I8*)p, size * sizeof(T));                   \
}                                                                                                 \
}
#else
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL_CAST(T)                                 \
namespace MTL                                                                                     \
{                                                                                                 \
template <> inline void OptimizedZeros(T* p, SizeType size)                                       \
{                                                                                                 \
  Parallel_1Dst< I8, OptimizedZeros_Sequential<I8> >((I8*)p, size * sizeof(T));                   \
}                                                                                                 \
}
#endif

// Initializes all elements in the vector buffer to zero if this is defined for the element.
#define MTL_DYNAMIC_VECTOR_ZERO_INIT(T)                                                           \
namespace MTL                                                                                     \
{                                                                                                 \
MTL_INLINE void FastDynamicVector<T>::FastZeroInit(T *p, SizeType size)                           \
{                                                                                                 \
  MTL::OptimizedZeros((I8*)p, size * sizeof(T));                                                  \
}                                                                                                 \
}

// Faster copy using intrinsic memcpy.
#define MTL_DYNAMIC_VECTOR_OPTIMIZED_COPY(T)                                                      \
namespace MTL                                                                                     \
{                                                                                                 \
template <>                                                                                       \
MTL_INLINE void DynamicVector<T>::OptimizedCopy(T* pDst, const T* pSrc, SizeType size)            \
{                                                                                                 \
  MTL::OptimizedCopy(pDst, pSrc, size);                                                           \
}                                                                                                 \
}

#define MTL_DYNAMIC_VECTOR_OPTIMIZED_COPY_CAST(T)                                                 \
namespace MTL                                                                                     \
{                                                                                                 \
template <>                                                                                       \
MTL_INLINE void DynamicVector<T>::OptimizedCopy(T* pDst, const T* pSrc, SizeType size)            \
{                                                                                                 \
  MTL::OptimizedCopy((I8*)pDst, (I8*)pSrc, size * sizeof(T));                                     \
}                                                                                                 \
}

// No initialization required for some cases.
#define MTL_DYNAMIC_VECTOR_NO_CONSTRUCTOR_DESTRUCTOR(T)                                           \
namespace MTL                                                                                     \
{                                                                                                 \
template <> MTL_INLINE void DynamicVector<T>::ConstructElements(T*, const T*) {}                  \
template <> MTL_INLINE void DynamicVector<T>::DestroyElements  (T*, const T*) {}                  \
}


#define MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(T)       \
MTL_DYNAMIC_VECTOR_OPTIMIZED_COPY(T)                  \
MTL_DYNAMIC_VECTOR_NO_CONSTRUCTOR_DESTRUCTOR(T)       \
MTL_DYNAMIC_VECTOR_OPTIMIZED_ASSIGN_ALL(T)            \
MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL(T)

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(I8)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(U8)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(I16)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(U16)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(I32)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(U32)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(I64)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(U64)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(F32)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS(F64)

#define MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(T)       \
MTL_DYNAMIC_VECTOR_OPTIMIZED_COPY_CAST(T)                  \
MTL_DYNAMIC_VECTOR_NO_CONSTRUCTOR_DESTRUCTOR(T)            \
MTL_DYNAMIC_VECTOR_OPTIMIZED_ZEROS_USE_ASSIGN_ALL_CAST(T)

MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector1D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector2D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector3D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector4D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector5D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(ColumnVector6D)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix1x1)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix2x2)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix3x3)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix4x4)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix5x5)
MTL_DYNAMIC_VECTOR_ALL_OPTIMIZATIONS_CAST(SquareMatrix6x6)

#endif  // MTL_DYNAMIC_VECTOR_H
