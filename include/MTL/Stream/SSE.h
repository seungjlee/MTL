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


#ifndef MTL_SSE_H
#define MTL_SSE_H

#include <MTL/Math.h>
#include <MTL/Stream/Stream.h>

#if MTL_ENABLE_SSE

#include <emmintrin.h>

namespace MTL
{

// Set all values in the XMM register with the same single input value.
MTL_INLINE static __m128  X128_SetPacked(F32 val)    { return _mm_set1_ps(val);                    }
MTL_INLINE static __m128d X128_SetPacked(F64 val)    { return _mm_set1_pd(val);                    }
MTL_INLINE static __m128i X128_SetPacked(I8  val)    { return _mm_set1_epi8(val);                  }
MTL_INLINE static __m128i X128_SetPacked(I16 val)    { return _mm_set1_epi16(val);                 }
MTL_INLINE static __m128i X128_SetPacked(I32 val)    { return _mm_set1_epi32(val);                 }
MTL_INLINE static __m128i X128_SetPacked(U8  val)    { return _mm_set1_epi8(val);                  }
MTL_INLINE static __m128i X128_SetPacked(U16 val)    { return _mm_set1_epi16(val);                 }
MTL_INLINE static __m128i X128_SetPacked(U32 val)    { return _mm_set1_epi32(val);                 }
MTL_INLINE static __m128i X128_SetPacked(I64 val)
{
  __m128d xx = X128_SetPacked((F64&)val);
  return (__m128i&)xx;
}
MTL_INLINE static __m128i X128_SetPacked(U64 val)
{
  __m128d xx = X128_SetPacked((F64&)val);
  return (__m128i&)xx;
}

// Common constants.
static const __m128  kX128_OnesF32       = X128_SetPacked(1.0f);
static const __m128  kX128_ThreeHalves32 = X128_SetPacked(1.5f);
static const __m128  kX128_TwosF32       = X128_SetPacked(2.0f);
static const __m128  kX128_ZerosF32      = _mm_setzero_ps();
static const __m128  kX128_HalvesF32     = X128_SetPacked((F32)kHalf);
static const __m128d kX128_OnesF64       = X128_SetPacked(1.0 );
static const __m128d kX128_TwosF64       = X128_SetPacked(2.0 );
static const __m128d kX128_ZerosF64      = _mm_setzero_pd();
static const __m128d kX128_HalvesF64     = X128_SetPacked(kHalf);
static const __m128  kX128_SignF32       = X128_SetPacked(( F32&)kSign32  );
static const __m128  kX128_NoSignF32     = X128_SetPacked(( F32&)kNoSign32);
static const __m128d kX128_SignF64       = X128_SetPacked((double&)kSign64  );
static const __m128d kX128_NoSignF64     = X128_SetPacked((double&)kNoSign64);
static const __m128i kX128_ZerosI        = _mm_setzero_si128();
static const __m128i kX128_NoZerosI      = X128_SetPacked(I32(-1));

template <class T> struct X128_Type {};
template <> struct X128_Type<F64>  { typedef __m128d Type; };
template <> struct X128_Type<F32>  { typedef __m128  Type; };
template <> struct X128_Type<I8>   { typedef __m128i Type; };
template <> struct X128_Type<I16>  { typedef __m128i Type; };
template <> struct X128_Type<I32>  { typedef __m128i Type; };
template <> struct X128_Type<I64>  { typedef __m128i Type; };
template <> struct X128_Type<U8>   { typedef __m128i Type; };
template <> struct X128_Type<U16>  { typedef __m128i Type; };
template <> struct X128_Type<U32>  { typedef __m128i Type; };
template <> struct X128_Type<U64>  { typedef __m128i Type; };


// Base class.
template <class T> class X128_Base : public X128_Type<T>
{
public:
  typedef typename X128_Base<T>::Type DataType;

  MTL_INLINE X128_Base() {}
  MTL_INLINE X128_Base(const DataType& data) : Data_(data) {}

  enum
  {
    Increment = sizeof(DataType) / sizeof(T)
  };

  MTL_INLINE static SizeType StreamSize(SizeType size)  { return size & ~(Increment-1); }

  MTL_INLINE static bool DoSSE(SizeType size)           { return size >= 2 * Increment; }

  MTL_INLINE const DataType& Data() const               { return Data_;                 }
  MTL_INLINE const T* pData() const                     { return (T*)&Data_;            }
  MTL_INLINE const T& operator[](SizeType i) const      { return pData()[i];            }
  MTL_INLINE       T& operator[](SizeType i)            { return ((T*)pData())[i];      }

protected:
  DataType Data_;
};

template <class T> class X128;

template<> class X128<I8> : public X128_Base<I8>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<I8>(data) {}
  MTL_INLINE X128(I8 val)          { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const I8 *ptr)   { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I8* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I8* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I8* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I8* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I8* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I8* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }
};

template<> class X128<U8> : public X128_Base<U8>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<U8>(data) {}
  MTL_INLINE X128(U8 val)          { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const U8 *ptr)   { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U8* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U8* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U8* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U8* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U8* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U8* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }
};

template<> class X128<I32> : public X128_Base<I32>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<I32>(data) {}
  MTL_INLINE X128(I32 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const I32 *ptr)  { LoadPackedUnaligned(ptr); }
  MTL_INLINE X128(I32 val3, I32 val2, I32 val1, I32 val0)  { Set(val3, val2, val1, val0); }

  MTL_INLINE void Set(I32 val3, I32 val2, I32 val1, I32 val0)
  {
    Data_ = _mm_set_epi32(val3, val2, val1, val0);
  }

  MTL_INLINE void Load(const I32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I32* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I32* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I32* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I32* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator-() const                { return _mm_sub_epi32(kX128_ZerosI, Data_); }
  MTL_INLINE X128 operator==(const X128& y) const  { return _mm_cmpeq_epi32(Data_, y.Data_);    }
  MTL_INLINE X128 operator<(const X128& y) const   { return _mm_cmplt_epi32(Data_, y.Data_);    }
  MTL_INLINE X128 operator>(const X128& y) const   { return _mm_cmpgt_epi32(Data_, y.Data_);    }
  MTL_INLINE X128 operator+(const X128& y) const   { return _mm_add_epi32(Data_, y.Data_);      }
  MTL_INLINE X128 operator-(const X128& y) const   { return _mm_sub_epi32(Data_, y.Data_);      }
  MTL_INLINE X128 operator&(const X128& y) const   { return _mm_and_si128(Data_, y.Data_);      }
  MTL_INLINE X128 operator|(const X128& y) const   { return _mm_or_si128(Data_, y.Data_);       }
  MTL_INLINE X128 operator^(const X128& y) const   { return _mm_xor_si128(Data_, y.Data_);      }

  MTL_INLINE X128 operator>>(int shift) const      { return _mm_srai_epi32(Data_, shift);       }
  MTL_INLINE X128 operator<<(int shift) const      { return _mm_slli_epi32(Data_, shift);       }
  MTL_INLINE X128& operator>>=(int shift)          { return *this = *this >> shift;             }
  MTL_INLINE X128& operator<<=(int shift)          { return *this = *this << shift;             }

  MTL_STREAM_EXTRA_INTEGER_OPERATORS(128);
};

template<> class X128<U32> : public X128_Base<U32>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<U32>(data) {}
  MTL_INLINE X128(U32 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const U32 *ptr)  { LoadPackedUnaligned(ptr); }
  MTL_INLINE X128(U32 val3, U32 val2, U32 val1, U32 val0)  { Set(val3, val2, val1, val0); }

  MTL_INLINE void Set(U32 val3, U32 val2, U32 val1, U32 val0)
  {
    Data_ = _mm_set_epi32(val3, val2, val1, val0);
  }

  MTL_INLINE void Load(const U32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U32* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U32* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U32* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U32* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator+(const X128& y) const  { return _mm_add_epi32(Data_, y.Data_); }
  MTL_INLINE X128 operator-(const X128& y) const  { return _mm_sub_epi32(Data_, y.Data_); }
  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }

  MTL_INLINE X128 operator>>(int shift) const     { return _mm_srli_epi32(Data_, shift);  }
  MTL_INLINE X128 operator<<(int shift) const     { return _mm_slli_epi32(Data_, shift);  }
  MTL_INLINE X128& operator>>=(int shift)         { return *this = *this >> shift;        }
  MTL_INLINE X128& operator<<=(int shift)         { return *this = *this << shift;        }

  MTL_STREAM_EXTRA_INTEGER_OPERATORS(128);
};

template<> class X128<I16> : public X128_Base<I16>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<I16>(data) {}
  MTL_INLINE X128(I16 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const I16 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I16* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I16* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I16* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I16* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I16* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I16* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const    { return _mm_and_si128(Data_, y.Data_);    }
  MTL_INLINE X128 operator|(const X128& y) const    { return _mm_or_si128(Data_, y.Data_);     }
  MTL_INLINE X128 operator^(const X128& y) const    { return _mm_xor_si128(Data_, y.Data_);    }

  // Returns the low 16-bits of the 32-bit products.
  MTL_INLINE X128  operator* (const X128& y) const  { return _mm_mullo_epi16(Data_, y.Data_);  }
  MTL_INLINE X128& operator*=(const X128& y)        { return *this = *this * y;                }

  // Returns the high 16-bits of the 32-bit products.
  MTL_INLINE X128  MultiplyHi(const X128& y) const  { return _mm_mulhi_epi16(Data_, y.Data_);  }

  // Converts low/high 4 16-bit signed integers into 32-bit signed integers.
  MTL_INLINE X128<I32> loToI32() const
  { return _mm_unpacklo_epi16(Data_, _mm_cmplt_epi16(Data_, kX128_ZerosI)); }
  MTL_INLINE X128<I32> hiToI32() const
  { return _mm_unpackhi_epi16(Data_, _mm_cmplt_epi16(Data_, kX128_ZerosI)); }

  MTL_INLINE void packSaturate(const X128<I32>& lo, const X128<I32>& hi)
  { Data_ = _mm_packs_epi32(lo.Data(), hi.Data()); }
};

template<> class X128<U16> : public X128_Base<U16>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<U16>(data) {}
  MTL_INLINE X128(U16 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const U16 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U16* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U16* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U16* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U16* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U16* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U16* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }
};

template<> class X128<I64> : public X128_Base<I64>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<I64>(data) {}
  MTL_INLINE X128(I64 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const I64 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I64* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I64* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I64* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I64* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }
};

template<> class X128<U64> : public X128_Base<U64>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<U64>(data) {}
  MTL_INLINE X128(U64 val)         { Data_ = X128_SetPacked(val); }
  MTL_INLINE X128(const U64 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U64* pSrc)
  { Data_ = _mm_loadu_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U64* pDst) const
  { _mm_storeu_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U64* pSrc)
  { Data_ = _mm_load_si128((const __m128i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U64* pDst) const
  { _mm_store_si128((__m128i*)pDst, Data_);       }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosI;  }

  MTL_INLINE X128 operator&(const X128& y) const  { return _mm_and_si128(Data_, y.Data_); }
  MTL_INLINE X128 operator|(const X128& y) const  { return _mm_or_si128(Data_, y.Data_);  }
  MTL_INLINE X128 operator^(const X128& y) const  { return _mm_xor_si128(Data_, y.Data_); }
};

template<> class X128<F32> : public X128_Base<F32>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<F32>(data) {}
  MTL_INLINE X128(F32 val)                                 { Set(val);                    }
  MTL_INLINE X128(F32 val3, F32 val2, F32 val1, F32 val0)  { Set(val3, val2, val1, val0); }
  MTL_INLINE X128(const F32 *ptr)                          { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Set(F32 val)                             { Data_ = X128_SetPacked(val); }
  MTL_INLINE void Set(F32 val3, F32 val2, F32 val1, F32 val0)
  {
    Data_ = _mm_setr_ps(val3, val2, val1, val0);
  }

  MTL_INLINE explicit X128(const X128<I32>& x)    { Data_ = _mm_cvtepi32_ps(x.Data());         }
  MTL_INLINE X128<I32> RoundedIntegers() const    { return X128<I32>(_mm_cvtps_epi32(Data_));  }
  MTL_INLINE X128<I32> TruncatedIntegers() const  { return X128<I32>(_mm_cvttps_epi32(Data_)); }

  // Conversion helpers.
  MTL_INLINE explicit X128(F64 val)  { Set((F32)val); }
  MTL_INLINE explicit X128(F64 val3, F64 val2, F64 val1, F64 val0)
  { Set((F32)val3, (F32)val2, (F32)val1, (F32)val0); }

  MTL_INLINE static X128 Zeros()        { return kX128_ZerosF32 ;     }
  MTL_INLINE static X128 Ones()         { return kX128_OnesF32;       }
  MTL_INLINE static X128 Halves()       { return kX128_HalvesF32;     }
  MTL_INLINE static X128 ThreeHalves()  { return kX128_ThreeHalves32; }
  MTL_INLINE static X128 Twos()         { return kX128_TwosF32;       }

  MTL_INLINE void Load(const F32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(F32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const F32* pSrc)   { Data_ = _mm_loadu_ps(pSrc); }
  MTL_INLINE void StorePackedUnaligned(F32* pDst) const  { _mm_storeu_ps(pDst, Data_); }
  MTL_INLINE void LoadPackedAligned(const F32* pSrc)     { Data_ = _mm_load_ps(pSrc);  }
  MTL_INLINE void StorePackedAligned(F32* pDst) const    { _mm_store_ps(pDst, Data_);  }

  MTL_INLINE void LoadSingle(const F32* pSrc)   { Data_ = _mm_load_ss(pSrc);  }
  MTL_INLINE void StoreSingle(F32* pDst) const  { _mm_store_ss(pDst, Data_);  }

  MTL_INLINE X128 operator==(const X128& y) const  { return _mm_cmpeq_ps(Data_, y.Data_);       }
  MTL_INLINE X128 operator<(const X128& y)  const  { return _mm_cmplt_ps(Data_, y.Data_);       }
  MTL_INLINE X128 operator>(const X128& y)  const  { return _mm_cmpgt_ps(Data_, y.Data_);       }
  MTL_INLINE X128 operator<=(const X128& y) const  { return _mm_cmple_ps(Data_, y.Data_);       }
  MTL_INLINE X128 operator>=(const X128& y) const  { return _mm_cmpge_ps(Data_, y.Data_);       }
  MTL_INLINE X128 operator-()               const  { return _mm_xor_ps(Data_, kX128_SignF32);   }
  MTL_INLINE X128 operator+(const X128& y)  const  { return _mm_add_ps(Data_, y.Data_);         }
  MTL_INLINE X128 operator-(const X128& y)  const  { return _mm_sub_ps(Data_, y.Data_);         }
  MTL_INLINE X128 operator*(const X128& y)  const  { return _mm_mul_ps(Data_, y.Data_);         }
  MTL_INLINE X128 operator/(const X128& y)  const  { return _mm_div_ps(Data_, y.Data_);         }
  MTL_INLINE X128 operator&(const X128& y)  const  { return _mm_and_ps(Data_, y.Data_);         }
  MTL_INLINE X128 operator|(const X128& y)  const  { return _mm_or_ps(Data_, y.Data_);          }
  MTL_INLINE X128 operator^(const X128& y)  const  { return _mm_xor_ps(Data_, y.Data_);         }
  MTL_INLINE X128 SquareRoot()              const  { return _mm_sqrt_ps(Data_);                 }
  MTL_INLINE X128 SquareRootSingle()        const  { return _mm_sqrt_ss(Data_);                 }

  // Precision of about 11-bits.
  MTL_INLINE X128 Reciprocal()                  const  { return _mm_rcp_ps(Data_);   }
  MTL_INLINE X128 ReciprocalSquareRoot()        const  { return _mm_rsqrt_ps(Data_); }
  MTL_INLINE X128 ReciprocalSquareRootSingle()  const  { return _mm_rsqrt_ss(Data_); }

  // Precision of about 22-bits. Doing a Newton iteration to improve precision.
  MTL_INLINE X128 ReciprocalPrecise() const
  {
    X128 r = Reciprocal();
    return r * (Twos() - *this * r);
  }

  // Accuracy of about 22-bits. Doing a Newton iteration to improve precision.
  MTL_INLINE X128 ReciprocalSquareRootPrecise() const
  {
    X128 r = ReciprocalSquareRoot();
    return r * (ThreeHalves() - Halves() * *this * r * r);
  }
  MTL_INLINE X128 ReciprocalSquareRootSinglePrecise() const
  {
    X128 r = ReciprocalSquareRootSingle();
    return r * (ThreeHalves() - Halves() * *this * r * r);
  }

  MTL_STREAM_EXTRA_OPERATORS(128);
};

template<> class X128<F64> : public X128_Base<F64>
{
public:
  MTL_INLINE X128() : X128_Base() {}
  MTL_INLINE X128(const DataType& data) : X128_Base<double>(data) {}
  MTL_INLINE X128(F64 val)             { Set(val);                 }
  MTL_INLINE X128(F64 val1, F64 val0)  { Set(val1, val0);          }
  MTL_INLINE X128(const F64 *ptr)      { LoadPackedUnaligned(ptr); }

  MTL_INLINE void Set(F64 val)             { Data_ = X128_SetPacked(val);     }
  MTL_INLINE void Set(F64 val1, F64 val0)  { Data_ = _mm_setr_pd(val1, val0); }

  MTL_INLINE explicit X128(const X128<I32>& x)    { Data_ = _mm_cvtepi32_pd(x.Data());         }
  MTL_INLINE X128<I32> RoundedIntegers() const    { return X128<I32>(_mm_cvtpd_epi32(Data_));  }
  MTL_INLINE X128<I32> TruncatedIntegers() const  { return X128<I32>(_mm_cvttpd_epi32(Data_)); }

  MTL_INLINE static X128 Zeros()   { return kX128_ZerosF64;  }
  MTL_INLINE static X128 Ones()    { return kX128_OnesF64;   }
  MTL_INLINE static X128 Halves()  { return kX128_HalvesF64; }
  MTL_INLINE static X128 Twos()    { return kX128_TwosF64;   }

  MTL_INLINE void Load(const F64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(F64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const F64* pSrc)   { Data_ = _mm_loadu_pd(pSrc); }
  MTL_INLINE void StorePackedUnaligned(F64* pDst) const  { _mm_storeu_pd(pDst, Data_); }
  MTL_INLINE void LoadPackedAligned(const F64* pSrc)     { Data_ = _mm_load_pd(pSrc);  }
  MTL_INLINE void StorePackedAligned(F64* pDst) const    { _mm_store_pd(pDst, Data_);  }

  MTL_INLINE void LoadSingle(const F64* pSrc)   { Data_ = _mm_load_sd(pSrc);  }
  MTL_INLINE void StoreSingle(F64* pDst) const  { _mm_store_sd(pDst, Data_);  }

  MTL_INLINE X128 operator==(const X128& y) const  { return _mm_cmpeq_pd(Data_, y.Data_);       }
  MTL_INLINE X128 operator<(const X128& y)  const  { return _mm_cmplt_pd(Data_, y.Data_);       }
  MTL_INLINE X128 operator>(const X128& y)  const  { return _mm_cmpgt_pd(Data_, y.Data_);       }
  MTL_INLINE X128 operator-() const                { return _mm_xor_pd(Data_, kX128_SignF64);   }
  MTL_INLINE X128 operator+(const X128& y)  const  { return _mm_add_pd(Data_, y.Data_);         }
  MTL_INLINE X128 operator-(const X128& y)  const  { return _mm_sub_pd(Data_, y.Data_);         }
  MTL_INLINE X128 operator*(const X128& y)  const  { return _mm_mul_pd(Data_, y.Data_);         }
  MTL_INLINE X128 operator/(const X128& y)  const  { return _mm_div_pd(Data_, y.Data_);         }
  MTL_INLINE X128 operator&(const X128& y)  const  { return _mm_and_pd(Data_, y.Data_);         }
  MTL_INLINE X128 operator|(const X128& y)  const  { return _mm_or_pd(Data_, y.Data_);          }
  MTL_INLINE X128 operator^(const X128& y)  const  { return _mm_xor_pd(Data_, y.Data_);         }
  MTL_INLINE X128 SquareRoot()              const  { return _mm_sqrt_pd(Data_);                 }
  MTL_INLINE X128 SquareRootSingle()        const  { return _mm_sqrt_sd(kX128_ZerosF64, Data_); }

  MTL_STREAM_EXTRA_OPERATORS(128);
};

// Minimum and maximum.
template <> MTL_INLINE X128<F32> Min(const X128<F32>& a, const X128<F32>& b)
{
  return _mm_min_ps(a.Data(), b.Data());
}
template <> MTL_INLINE X128<F64> Min(const X128<F64>& a, const X128<F64>& b)
{
  return _mm_min_pd(a.Data(), b.Data());
}
template <> MTL_INLINE X128<U8> Min(const X128<U8>& a, const X128<U8>& b)
{
  return _mm_min_epu8(a.Data(), b.Data());
}
template <> MTL_INLINE X128<I16> Min(const X128<I16>& a, const X128<I16>& b)
{
  return _mm_min_epi16(a.Data(), b.Data());
}
template <> MTL_INLINE X128<F32> Max(const X128<F32>& a, const X128<F32>& b)
{
  return _mm_max_ps(a.Data(), b.Data());
}
template <> MTL_INLINE X128<F64> Max(const X128<F64>& a, const X128<F64>& b)
{
  return _mm_max_pd(a.Data(), b.Data());
}
template <> MTL_INLINE X128<U8> Max(const X128<U8>& a, const X128<U8>& b)
{
  return _mm_max_epu8(a.Data(), b.Data());
}
template <> MTL_INLINE X128<I16> Max(const X128<I16>& a, const X128<I16>& b)
{
  return _mm_max_epi16(a.Data(), b.Data());
}

// Absolute value.
MTL_INLINE X128<F64> Abs(const X128<F64>& a)
{
  return a & kX128_NoSignF64;
}
MTL_INLINE X128<F32> Abs(const X128<F32>& a)
{
  return a & kX128_NoSignF32;
}

template <> MTL_INLINE
X128<F64> MultiplyAndAdd(const X128<F64>& a, const X128<F64>& b, const X128<F64>& c)
{
  return a * b + c;
}
template <> MTL_INLINE
X128<F32> MultiplyAndAdd(const X128<F32>& a, const X128<F32>& b, const X128<F32>& c)
{
  return a * b + c;
}

// Conditional equivalent to "condition ? a : b" in C/C++.
MTL_INLINE X128<I8> Conditional(const X128<I8>& condition, const X128<I8>& a, const X128<I8>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<U8> Conditional(const X128<U8>& condition, const X128<U8>& a, const X128<U8>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<I16> Conditional(const X128<I16>& condition, const X128<I16>& a, const X128<I16>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<U16> Conditional(const X128<U16>& condition, const X128<U16>& a, const X128<U16>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<I32> Conditional(const X128<I32>& condition, const X128<I32>& a, const X128<I32>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<U32> Conditional(const X128<U32>& condition, const X128<U32>& a, const X128<U32>& b)
{
  return _mm_or_si128(_mm_and_si128(condition.Data(), a.Data()),
                      _mm_andnot_si128(condition.Data(), b.Data()));;
}
MTL_INLINE X128<F32> Conditional(const X128<F32>& condition,
                                 const X128<F32>& a, const X128<F32>& b)
{
  return _mm_or_ps(_mm_and_ps(condition.Data(), a.Data()),
                   _mm_andnot_ps(condition.Data(), b.Data()));;
}
MTL_INLINE X128<double> Conditional(const X128<F64>& condition,
                                    const X128<F64>& a, const X128<F64>& b)
{
  return _mm_or_pd(_mm_and_pd(condition.Data(), a.Data()),
                   _mm_andnot_pd(condition.Data(), b.Data()));;
}

// Permute values in registers.
template<int PermutationBits>
MTL_INLINE X128<I32> Shuffle(const X128<I32>& a)
{
  return _mm_shuffle_epi32(a.Data(), PermutationBits);
}
template<int PermutationBits>
MTL_INLINE X128<U32> Shuffle(const X128<U32>& a)
{
  return _mm_shuffle_epi32(a.Data(), PermutationBits);
}
template<int PermutationBits>
MTL_INLINE X128<F32> Shuffle(const X128<F32>& a, const X128<F32>& b)
{
  return _mm_shuffle_ps(a.Data(), b.Data(), PermutationBits);
}
template<int PermutationBits>
MTL_INLINE X128<F64> Shuffle(const X128<F64>& a, const X128<F64>& b)
{
  return _mm_shuffle_pd(a.Data(), b.Data(), PermutationBits);
}


template<int SHIFT>
MTL_INLINE X128<U32> ShiftRight(const X128<U32>& a)  { return _mm_srli_epi32(a.Data(), SHIFT); }
template<int SHIFT>
MTL_INLINE X128<U32> ShiftLeft (const X128<U32>& a)  { return _mm_slli_epi32(a.Data(), SHIFT); }
template<int SHIFT>
MTL_INLINE X128<U32> RotateLeft(const X128<U32>& a)  { return ShiftLeft<SHIFT>(a) | ShiftRight<32 - SHIFT>(a); }

}  // namespace MTL

#endif  // MTL_ENABLE_SSE

#endif  // MTL_SSE_H
