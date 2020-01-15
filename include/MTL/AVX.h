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


#ifndef MTL_AVX_H
#define MTL_AVX_H

#include "Math.h"
#include "Stream.h"

#if MTL_ENABLE_AVX

#include <immintrin.h>

namespace MTL
{
// Set all values in the XMM register with the same single input value.
MTL_INLINE static __m256  X256_SetPacked(F32 val)    { return _mm256_set1_ps(val);                 }
MTL_INLINE static __m256d X256_SetPacked(F64 val)    { return _mm256_set1_pd(val);                 }
MTL_INLINE static __m256i X256_SetPacked(I8  val)    { return _mm256_set1_epi8(val);               }
MTL_INLINE static __m256i X256_SetPacked(I16 val)    { return _mm256_set1_epi16(val);              }
MTL_INLINE static __m256i X256_SetPacked(I32 val)    { return _mm256_set1_epi32(val);              }
MTL_INLINE static __m256i X256_SetPacked(U8  val)    { return _mm256_set1_epi8(val);               }
MTL_INLINE static __m256i X256_SetPacked(U16 val)    { return _mm256_set1_epi16(val);              }
MTL_INLINE static __m256i X256_SetPacked(U32 val)    { return _mm256_set1_epi32(val);              }
MTL_INLINE static __m256i X256_SetPacked(I64 val)    { return (__m256i&)X256_SetPacked((F64&)val); }
MTL_INLINE static __m256i X256_SetPacked(U64 val)    { return (__m256i&)X256_SetPacked((F64&)val); }

// Common constants.
static const __m256  kX256_OnesF32       = X256_SetPacked(1.0f);
static const __m256  kX256_ThreeHalves32 = X256_SetPacked(1.5f);
static const __m256  kX256_TwosF32       = X256_SetPacked(2.0f);
static const __m256  kX256_ZerosF32      = _mm256_setzero_ps();
static const __m256  kX256_HalvesF32     = X256_SetPacked((F32)kHalf);
static const __m256d kX256_OnesF64       = X256_SetPacked(1.0 );
static const __m256d kX256_TwosF64       = X256_SetPacked(2.0 );
static const __m256d kX256_ZerosF64      = _mm256_setzero_pd();
static const __m256d kX256_HalvesF64     = X256_SetPacked(kHalf);
static const __m256  kX256_SignF32       = X256_SetPacked(*((F32*)&kSign32)  );
static const __m256  kX256_NoSignF32     = X256_SetPacked(*((F32*)&kNoSign32));
static const __m256d kX256_SignF64       = X256_SetPacked(*((F64*)&kSign64  ));
static const __m256d kX256_NoSignF64     = X256_SetPacked(*((F64*)&kNoSign64));
static const __m256i kX256_ZerosI        = _mm256_setzero_si256();
static const __m256i kX256_NoZerosI      = X256_SetPacked(I32(-1));

template <class T> struct X256_Type {};
template <> struct X256_Type<F64>  { typedef __m256d Type; };
template <> struct X256_Type<F32>  { typedef __m256  Type; };
template <> struct X256_Type<I8>   { typedef __m256i Type; };
template <> struct X256_Type<I16>  { typedef __m256i Type; };
template <> struct X256_Type<I32>  { typedef __m256i Type; };
template <> struct X256_Type<I64>  { typedef __m256i Type; };
template <> struct X256_Type<U8>   { typedef __m256i Type; };
template <> struct X256_Type<U16>  { typedef __m256i Type; };
template <> struct X256_Type<U32>  { typedef __m256i Type; };
template <> struct X256_Type<U64>  { typedef __m256i Type; };


// Base class.
template <class T> class X256_Base : public X256_Type<T>
{
public:
  typedef typename X256_Base<T>::Type DataType;

  MTL_INLINE X256_Base() {}
  MTL_INLINE X256_Base(const DataType& data) : Data_(data) {}

  enum
  {
    Increment = sizeof(Type) / sizeof(T)
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

template <class T> class X256;

template<> class X256<I8> : public X256_Base<I8>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<I8>(data) {}
  MTL_INLINE X256(I8 val)          { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const I8 *ptr)   { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I8* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I8* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I8* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I8* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I8* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I8* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<U8> : public X256_Base<U8>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<U8>(data) {}
  MTL_INLINE X256(U8 val)          { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const U8 *ptr)   { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U8* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U8* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U8* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U8* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U8* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U8* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<I32> : public X256_Base<I32>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<I32>(data) {}
  MTL_INLINE X256(I32 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const I32 *ptr)  { LoadPackedUnaligned(ptr); }

  MTL_INLINE void Load(const I32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I32* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I32* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I32* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I32* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator-() const                { return _mm256_sub_epi32(kX256_ZerosI, Data_); }
  MTL_INLINE X256 operator==(const X256& y) const  { return _mm256_cmpeq_epi32(Data_, y.Data_);    }
  MTL_INLINE X256 operator<(const X256& y) const   { return _mm256_cmpgt_epi32(y.Data_, Data_);    }
  MTL_INLINE X256 operator>(const X256& y) const   { return _mm256_cmpgt_epi32(Data_, y.Data_);    }
  MTL_INLINE X256 operator+(const X256& y) const   { return _mm256_add_epi32(Data_, y.Data_);      }
  MTL_INLINE X256 operator-(const X256& y) const   { return _mm256_sub_epi32(Data_, y.Data_);      }
  MTL_INLINE X256 operator&(const X256& y) const   { return _mm256_and_si256(Data_, y.Data_);      }
  MTL_INLINE X256 operator|(const X256& y) const   { return _mm256_or_si256(Data_, y.Data_);       }
  MTL_INLINE X256 operator^(const X256& y) const   { return _mm256_xor_si256(Data_, y.Data_);      }

  MTL_INLINE X256 operator>>(int shift) const      { return _mm256_srai_epi32(Data_, shift);       }
  MTL_INLINE X256 operator<<(int shift) const      { return _mm256_slli_epi32(Data_, shift);       }
  MTL_INLINE X256& operator>>=(int shift)          { return *this = *this >> shift;                }
  MTL_INLINE X256& operator<<=(int shift)          { return *this = *this << shift;                }

  MTL_STREAM_EXTRA_INTEGER_OPERATORS(256);
};

template<> class X256<U32> : public X256_Base<U32>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<U32>(data) {}
  MTL_INLINE X256(U32 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const U32 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U32* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U32* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U32* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U32* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<I16> : public X256_Base<I16>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<I16>(data) {}
  MTL_INLINE X256(I16 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const I16 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I16* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I16* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I16* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I16* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I16* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I16* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const    { return _mm256_and_si256(Data_, y.Data_);   }
  MTL_INLINE X256 operator|(const X256& y) const    { return _mm256_or_si256(Data_, y.Data_);    }
  MTL_INLINE X256 operator^(const X256& y) const    { return _mm256_xor_si256(Data_, y.Data_);   }

  // Returns the low 16-bits of the 32-bit products.
  MTL_INLINE X256  operator* (const X256& y) const  { return _mm256_mullo_epi16(Data_, y.Data_); }
  MTL_INLINE X256& operator*=(const X256& y)        { return *this = *this * y;                  }

  // Returns the high 16-bits of the 32-bit products.
  MTL_INLINE X256  MultiplyHi(const X256& y) const  { return _mm256_mulhi_epi16(Data_, y.Data_); }

  // Converts low/high 4 16-bit signed integers into 32-bit signed integers.
  MTL_INLINE X256<I32> loToI32() const
  { return _mm256_unpacklo_epi16(Data_, _mm256_cmpgt_epi16(kX256_ZerosI, Data_)); }
  MTL_INLINE X256<I32> hiToI32() const
  { return _mm256_unpackhi_epi16(Data_, _mm256_cmpgt_epi16(kX256_ZerosI, Data_)); }

  MTL_INLINE void packSaturate(const X256<I32>& lo, const X256<I32>& hi)
  { Data_ = _mm256_packs_epi32(lo.Data(), hi.Data()); }
};

template<> class X256<U16> : public X256_Base<U16>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<U16>(data) {}
  MTL_INLINE X256(U16 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const U16 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U16* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U16* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U16* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U16* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U16* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U16* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<I64> : public X256_Base<I64>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<I64>(data) {}
  MTL_INLINE X256(I64 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const I64 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const I64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(I64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const I64* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(I64* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const I64* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(I64* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<U64> : public X256_Base<U64>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<U64>(data) {}
  MTL_INLINE X256(U64 val)         { Data_ = X256_SetPacked(val); }
  MTL_INLINE X256(const U64 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Load(const U64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(U64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const U64* pSrc)
  { Data_ = _mm256_loadu_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedUnaligned(U64* pDst) const
  { _mm256_storeu_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE void LoadPackedAligned(const U64* pSrc)
  { Data_ = _mm256_load_si256((const __m256i*)pSrc); }
  MTL_INLINE void StorePackedAligned(U64* pDst) const
  { _mm256_store_si256((__m256i*)pDst, Data_);       }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosI;  }

  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_si256(Data_, y.Data_); }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_si256(Data_, y.Data_);  }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_si256(Data_, y.Data_); }
};

template<> class X256<F32> : public X256_Base<F32>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<F32>(data) {}
  MTL_INLINE X256(F32 val)         { Set(val);                    }
  MTL_INLINE X256(F32 val0, F32 val1, F32 val2, F32 val3,
                  F32 val4, F32 val5, F32 val6, F32 val7)
  { Set(val0, val1, val2, val3, val4, val5, val6, val7); }
  MTL_INLINE X256(const F32 *ptr)  { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Set(F32 val)    { Data_ = X256_SetPacked(val); }
  MTL_INLINE void Set(F32 val0, F32 val1, F32 val2, F32 val3,
                      F32 val4, F32 val5, F32 val6, F32 val7)
  { Data_ = _mm256_setr_ps(val0, val1, val2, val3, val4, val5, val6, val7); }

  MTL_INLINE explicit X256(const X256<I32>& x)    { Data_ = _mm256_cvtepi32_ps(x.Data());         }
  MTL_INLINE X256<I32> RoundedIntegers() const    { return X256<I32>(_mm256_cvtps_epi32(Data_));  }
  MTL_INLINE X256<I32> TruncatedIntegers() const  { return X256<I32>(_mm256_cvttps_epi32(Data_)); }

  // Conversion helpers.
  MTL_INLINE explicit X256(F64 val)  { Set((F32)val); }
  MTL_INLINE explicit X256(F64 val0, F64 val1, F64 val2, F64 val3,
                           F64 val4, F64 val5, F64 val6, F64 val7)
  { Set((F32)val0, (F32)val1, (F32)val2, (F32)val3, (F32)val4, (F32)val5, (F32)val6, (F32)val7); }

  MTL_INLINE static X256 Zeros()        { return kX256_ZerosF32;      }
  MTL_INLINE static X256 Ones()         { return kX256_OnesF32;       }
  MTL_INLINE static X256 Halves()       { return kX256_HalvesF32;     }
  MTL_INLINE static X256 ThreeHalves()  { return kX256_ThreeHalves32; }
  MTL_INLINE static X256 Twos()         { return kX256_TwosF32;       }

  MTL_INLINE void Load(const F32* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(F32* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const F32* pSrc)   { Data_ = _mm256_loadu_ps(pSrc); }
  MTL_INLINE void StorePackedUnaligned(F32* pDst) const  { _mm256_storeu_ps(pDst, Data_); }
  MTL_INLINE void LoadPackedAligned(const F32* pSrc)     { Data_ = _mm256_load_ps(pSrc);  }
  MTL_INLINE void StorePackedAligned(F32* pDst) const    { _mm256_store_ps(pDst, Data_);  }

  MTL_INLINE X256 operator==(const X256& y) const  { return _mm256_cmp_ps(Data_, y.Data_, _CMP_EQ_OQ); }
  MTL_INLINE X256 operator< (const X256& y) const  { return _mm256_cmp_ps(Data_, y.Data_, _CMP_LT_OQ); }
  MTL_INLINE X256 operator> (const X256& y) const  { return _mm256_cmp_ps(Data_, y.Data_, _CMP_GT_OQ); }
  MTL_INLINE X256 operator<=(const X256& y) const  { return _mm256_cmp_ps(Data_, y.Data_, _CMP_LE_OQ); }
  MTL_INLINE X256 operator>=(const X256& y) const  { return _mm256_cmp_ps(Data_, y.Data_, _CMP_GE_OQ); }

  MTL_INLINE X256 operator-() const               { return _mm256_xor_ps(Data_, kX256_SignF32);   }
  MTL_INLINE X256 operator+(const X256& y) const  { return _mm256_add_ps(Data_, y.Data_);         }
  MTL_INLINE X256 operator-(const X256& y) const  { return _mm256_sub_ps(Data_, y.Data_);         }
  MTL_INLINE X256 operator*(const X256& y) const  { return _mm256_mul_ps(Data_, y.Data_);         }
  MTL_INLINE X256 operator/(const X256& y) const  { return _mm256_div_ps(Data_, y.Data_);         }
  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_ps(Data_, y.Data_);         }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_ps(Data_, y.Data_);          }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_ps(Data_, y.Data_);         }
  MTL_INLINE X256 SquareRoot()             const  { return _mm256_sqrt_ps(Data_);                 }

  // Precision of about 11-bits.
  MTL_INLINE X256 Reciprocal()                  const  { return _mm256_rcp_ps(Data_);   }
  MTL_INLINE X256 ReciprocalSquareRoot()        const  { return _mm256_rsqrt_ps(Data_); }

  // Precision of about 22-bits. Doing a Newton iteration to improve precision.
  MTL_INLINE X256 ReciprocalPrecise() const
  {
    X256 r = Reciprocal();
    return r * (Twos() - *this * r);
  }

  // Accuracy of about 22-bits. Doing a Newton iteration to improve precision.
  MTL_INLINE X256 ReciprocalSquareRootPrecise() const
  {
    X256 r = ReciprocalSquareRoot();
    return r * (ThreeHalves() - Halves() * *this * r * r);
  }

  MTL_STREAM_EXTRA_OPERATORS(256);
};

template<> class X256<F64> : public X256_Base<F64>
{
public:
  MTL_INLINE X256() : X256_Base() {}
  MTL_INLINE X256(const DataType& data) : X256_Base<double>(data) {}
  MTL_INLINE X256(F64 val)                                 { Set(val);                    }
  MTL_INLINE X256(F64 val0, F64 val1, F64 val2, F64 val3)  { Set(val0, val1, val2, val3); }
  MTL_INLINE X256(const F64 *ptr)                          { LoadPackedUnaligned(ptr);    }

  MTL_INLINE void Set(F64 val)             { Data_ = X256_SetPacked(val);     }
  MTL_INLINE void Set(F64 val0, F64 val1, F64 val2, F64 val3)
  { Data_ = _mm256_setr_pd(val0, val1, val2, val3); }

  MTL_INLINE explicit X256(const X128<I32>& x)    { Data_ = _mm256_cvtepi32_pd(x.Data());         }
  MTL_INLINE X128<I32> RoundedIntegers() const    { return X128<I32>(_mm256_cvtpd_epi32(Data_));  }
  MTL_INLINE X128<I32> TruncatedIntegers() const  { return X128<I32>(_mm256_cvttpd_epi32(Data_)); }

  MTL_INLINE static X256 Zeros()   { return kX256_ZerosF64;  }
  MTL_INLINE static X256 Ones()    { return kX256_OnesF64;   }
  MTL_INLINE static X256 Halves()  { return kX256_HalvesF64; }
  MTL_INLINE static X256 Twos()    { return kX256_TwosF64;   }

  MTL_INLINE void Load(const F64* pSrc)   { LoadPackedUnaligned(pSrc);  }
  MTL_INLINE void Store(F64* pDst) const  { StorePackedUnaligned(pDst); }
  MTL_INLINE void LoadPackedUnaligned(const F64* pSrc)   { Data_ = _mm256_loadu_pd(pSrc); }
  MTL_INLINE void StorePackedUnaligned(F64* pDst) const  { _mm256_storeu_pd(pDst, Data_); }
  MTL_INLINE void LoadPackedAligned(const F64* pSrc)     { Data_ = _mm256_load_pd(pSrc);  }
  MTL_INLINE void StorePackedAligned(F64* pDst) const    { _mm256_store_pd(pDst, Data_);  }

  MTL_INLINE X256 operator==(const X256& y) const 
  { return _mm256_cmp_pd(Data_, y.Data_, _CMP_EQ_OQ); }
  MTL_INLINE X256 operator<(const X256& y) const 
  { return _mm256_cmp_pd(Data_, y.Data_, _CMP_LT_OQ); }
  MTL_INLINE X256 operator>(const X256& y) const
  { return _mm256_cmp_pd(Data_, y.Data_, _CMP_GT_OQ); }
  MTL_INLINE X256 operator-() const               { return _mm256_xor_pd(Data_, kX256_SignF64);   }
  MTL_INLINE X256 operator+(const X256& y) const  { return _mm256_add_pd(Data_, y.Data_);         }
  MTL_INLINE X256 operator-(const X256& y) const  { return _mm256_sub_pd(Data_, y.Data_);         }
  MTL_INLINE X256 operator*(const X256& y) const  { return _mm256_mul_pd(Data_, y.Data_);         }
  MTL_INLINE X256 operator/(const X256& y) const  { return _mm256_div_pd(Data_, y.Data_);         }
  MTL_INLINE X256 operator&(const X256& y) const  { return _mm256_and_pd(Data_, y.Data_);         }
  MTL_INLINE X256 operator|(const X256& y) const  { return _mm256_or_pd(Data_, y.Data_);          }
  MTL_INLINE X256 operator^(const X256& y) const  { return _mm256_xor_pd(Data_, y.Data_);         }
  MTL_INLINE X256 SquareRoot()             const  { return _mm256_sqrt_pd(Data_);                 }

  MTL_STREAM_EXTRA_OPERATORS(256);
};

// Minimum and maximum.
template <> MTL_INLINE X256<F32> Min(const X256<F32>& a, const X256<F32>& b)
{
  return _mm256_min_ps(a.Data(), b.Data());
}
template <> MTL_INLINE X256<F64> Min(const X256<F64>& a, const X256<F64>& b)
{
  return _mm256_min_pd(a.Data(), b.Data());
}
template <> MTL_INLINE X256<U8> Min(const X256<U8>& a, const X256<U8>& b)
{
  return _mm256_min_epu8(a.Data(), b.Data());
}
template <> MTL_INLINE X256<I16> Min(const X256<I16>& a, const X256<I16>& b)
{
  return _mm256_min_epi16(a.Data(), b.Data());
}
template <> MTL_INLINE X256<F32> Max(const X256<F32>& a, const X256<F32>& b)
{
  return _mm256_max_ps(a.Data(), b.Data());
}
template <> MTL_INLINE X256<F64> Max(const X256<F64>& a, const X256<F64>& b)
{
  return _mm256_max_pd(a.Data(), b.Data());
}
template <> MTL_INLINE X256<U8> Max(const X256<U8>& a, const X256<U8>& b)
{
  return _mm256_max_epu8(a.Data(), b.Data());
}
template <> MTL_INLINE X256<I16> Max(const X256<I16>& a, const X256<I16>& b)
{
  return _mm256_max_epi16(a.Data(), b.Data());
}

// Absolute value.
MTL_INLINE X256<F64> Abs(const X256<F64>& a)
{
  return a & kX256_NoSignF64;
}
MTL_INLINE X256<F32> Abs(const X256<F32>& a)
{
  return a & kX256_NoSignF32;
}

template <> MTL_INLINE
X256<F64> MultiplyAndAdd(const X256<F64>& a, const X256<F64>& b, const X256<F64>& c)
{
  return _mm256_fmadd_pd(a.Data(), b.Data(), c.Data());
}
template <> MTL_INLINE
X256<F32> MultiplyAndAdd(const X256<F32>& a, const X256<F32>& b, const X256<F32>& c)
{
  return _mm256_fmadd_ps(a.Data(), b.Data(), c.Data());
}

// Conditional equivalent to "condition ? a : b" in C/C++.
MTL_INLINE X256<I8> Conditional(const X256<I8>& condition, const X256<I8>& a, const X256<I8>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<U8> Conditional(const X256<U8>& condition, const X256<U8>& a, const X256<U8>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<I16> Conditional(const X256<I16>& condition, const X256<I16>& a, const X256<I16>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<U16> Conditional(const X256<U16>& condition, const X256<U16>& a, const X256<U16>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<I32> Conditional(const X256<I32>& condition, const X256<I32>& a, const X256<I32>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<U32> Conditional(const X256<U32>& condition, const X256<U32>& a, const X256<U32>& b)
{
  return _mm256_or_si256(_mm256_and_si256(condition.Data(), a.Data()),
                         _mm256_andnot_si256(condition.Data(), b.Data()));;
}
MTL_INLINE X256<F32> Conditional(const X256<F32>& condition,
                                 const X256<F32>& a, const X256<F32>& b)
{
  return _mm256_or_ps(_mm256_and_ps(condition.Data(), a.Data()),
                      _mm256_andnot_ps(condition.Data(), b.Data()));;
}
MTL_INLINE X256<double> Conditional(const X256<F64>& condition,
                                    const X256<F64>& a, const X256<F64>& b)
{
  return _mm256_or_pd(_mm256_and_pd(condition.Data(), a.Data()),
                      _mm256_andnot_pd(condition.Data(), b.Data()));;
}

// Permute values in registers.
template<int PermutationBits>
MTL_INLINE X256<F32> Shuffle(const X256<F32>& a, const X256<F32>& b)
{
  return _mm256_shuffle_ps(a.Data(), b.Data(), PermutationBits);
}
template<int PermutationBits>
MTL_INLINE X256<F64> Shuffle(const X256<F64>& a, const X256<F64>& b)
{
  return _mm256_shuffle_pd(a.Data(), b.Data(), PermutationBits);
}

}  // namespace MTL

#endif  // MTL_ENABLE_AVX

#endif  // MTL_AVX_H
