//
// Math Template Library
//
// Copyright (c) 2026: Seung Jae Lee, https://github.com/seungjlee/MTL
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

// Targets cold paths in StreamArray, AVX intrinsic wrappers, the small-array
// fallback in Array.h, and special-case branches in Rotation3D that the wider
// test suite happens not to exercise.

#include <MTL/Array.h>
#include <MTL/Math/DynamicVector.h>
#include <MTL/Math/DynamicVectorOperators.h>
#include <MTL/Math/Rotation3D.h>
#include <MTL/Tools/Test.h>
#include <MTL/Stream/StreamArray.h>

#if MTL_ENABLE_AVX
#include <MTL/Stream/AVX.h>
#endif

using namespace MTL;

// Drive the full template surface of MTL::Array<N,T> at compile-and-call time:
// And/Or, Set with conversion, MaxNorm, etc. on a small fixed-size buffer.
TEST(StreamCoverage_ArraySmall)
{
  U32 buf[5] = { 0xF0F0F0F0u, 0x0F0F0F0Fu, 0xFFFF0000u, 0x0000FFFFu, 0xAAAA5555u };
  U32 a = Array<5, U32>::And(buf);
  U32 o = Array<5, U32>::Or(buf);
  MTL_EQUAL(a, U32(0));
  MTL_EQUAL(o, U32(0xFFFFFFFFu));

  // Set(T*, const T&) writes a scalar to every slot.
  F32 dst[4] = { 1, 2, 3, 4 };
  Array<4, F32>::Set(dst, 7.5f);
  for (int i = 0; i < 4; i++)
    MTL_EQUAL_FLOAT(dst[i], 7.5f, 0.0);

  // Set with type conversion.
  I32 src[4] = { 10, 20, 30, 40 };
  F64 dstD[4] = {};
  Array<4, F64>::Set(dstD, src);
  for (int i = 0; i < 4; i++)
    MTL_EQUAL_FLOAT(dstD[i], (F64)src[i], 0.0);

  F64 vals[4] = { -3.0, 1.5, -0.25, 2.0 };
  F64 maxNorm = Array<4, F64>::MaxNorm(vals);
  MTL_EQUAL_FLOAT(maxNorm, 3.0, 0.0);
}

// Drives RotationToXAxis / RotationToYAxis / RotationToZAxis through their
// axis-aligned special cases (the conditional with assert(v[i] != 0)).
TEST(StreamCoverage_Rotation3D_AxisAlignedSpecialCases)
{
  static const F64 kTol = 1e-14;

  // Each public RotationTo*Axis dispatches through the private 4-arg helper,
  // so feeding in axis-aligned inputs still drives the special-case branches.
  {
    Vector3D<F64> v(0.0, 1.5, 0.0);
    Rotation3D<F64> Rx = Rotation3D<F64>::RotationToXAxis(v);
    Point3D<F64> p = Rx * v;
    MTL_EQUAL_FLOAT(p.y(), 0.0, kTol);
    MTL_EQUAL_FLOAT(p.z(), 0.0, kTol);
  }
  {
    Vector3D<F64> v(2.0, 0.0, 0.0);
    Rotation3D<F64> Ry = Rotation3D<F64>::RotationToYAxis(v);
    Point3D<F64> p = Ry * v;
    MTL_EQUAL_FLOAT(p.x(), 0.0, kTol);
    MTL_EQUAL_FLOAT(p.z(), 0.0, kTol);
  }
  {
    Vector3D<F64> v(0.0, 0.5, 0.0);
    Rotation3D<F64> Rz = Rotation3D<F64>::RotationToZAxis(v);
    Point3D<F64> p = Rz * v;
    MTL_EQUAL_FLOAT(p.x(), 0.0, kTol);
    MTL_EQUAL_FLOAT(p.y(), 0.0, kTol);
  }
}

// Walk the StreamArray Min/Max/MinOfAbsolutes/MaxOfAbsolutes paths over a
// vector wider than one SIMD lane and not a clean multiple of XX::Increment,
// so the remaining-tail branch executes alongside the SIMD body.
template <class T>
static void DriveStreamReductions()
{
  enum { Size = 17 };  // > Increment for any X128/X256/X512 over F32/F64.

  DynamicVector<T> v(Size);
  for (SizeType i = 0; i < v.Size(); i++)
    v[i] = T(I32(i) - 8);                   // negative + zero + positive

  T expectedMin = T(-8);
  T expectedMax = T(Size - 1 - 8);
  T expectedMinAbs = T(0);
  T expectedMaxAbs = T(8);

  MTL_EQUAL_FLOAT(F64(Min(v)), F64(expectedMin), 0.0);
  MTL_EQUAL_FLOAT(F64(Max(v)), F64(expectedMax), 0.0);
  MTL_EQUAL_FLOAT(F64(MinOfAbsolutes(v)), F64(expectedMinAbs), 0.0);
  MTL_EQUAL_FLOAT(F64(MaxOfAbsolutes(v)), F64(expectedMaxAbs), 0.0);

  // Difference-norm path: SumOfSquaredDifferences (StreamArray line ~152-155).
  DynamicVector<T> u(Size);
  for (SizeType i = 0; i < u.Size(); i++)
    u[i] = v[i] + T(1);
  T ssd = SumOfSquaredDifferences(v, u);
  MTL_EQUAL_FLOAT(F64(ssd), F64(Size), 1e-6);

  // ScalarAddition aligned-stream path (StreamArray line ~570-574).
  DynamicVector<T> w(Size, T(2));
  w += 3.0;
  for (SizeType i = 0; i < w.Size(); i++)
    MTL_EQUAL_FLOAT(F64(w[i]), 5.0, 0.0);
}

TEST(StreamCoverage_ReductionsAndScalarOps_F32) { DriveStreamReductions<F32>(); }
TEST(StreamCoverage_ReductionsAndScalarOps_F64) { DriveStreamReductions<F64>(); }

#if MTL_ENABLE_AVX
// X256 wrappers around _mm256 intrinsics for F32/F64 that the wider suite
// happens not to call (Min/Max/FMA/aligned-load/U8 broadcast).
TEST(StreamCoverage_AVX_X256_FloatIntrinsics)
{
  alignas(32) F32 a32[8] = { 1, 5, 3, 7, 2, 6, 4, 8 };
  alignas(32) F32 b32[8] = { 4, 2, 6, 1, 8, 3, 5, 7 };
  X256<F32> a, b;
  a.LoadPackedAligned(a32);
  b.LoadPackedAligned(b32);
  X256<F32> mn = Min(a, b);
  X256<F32> mx = Max(a, b);
  for (int i = 0; i < 8; i++)
  {
    MTL_EQUAL_FLOAT(mn.pData()[i], MTL::Min(a32[i], b32[i]), 0.0);
    MTL_EQUAL_FLOAT(mx.pData()[i], MTL::Max(a32[i], b32[i]), 0.0);
  }

  X256<F32> c(0.5f);
  X256<F32> fma = MultiplyAndAdd(a, b, c);
  for (int i = 0; i < 8; i++)
    MTL_EQUAL_FLOAT(fma.pData()[i], a32[i] * b32[i] + 0.5f, 1e-5);

  alignas(32) F64 a64[4] = { 1, 5, 3, 7 };
  alignas(32) F64 b64[4] = { 4, 2, 6, 1 };
  X256<F64> ad, bd;
  ad.LoadPackedAligned(a64);
  bd.LoadPackedAligned(b64);
  X256<F64> mnd = Min(ad, bd);
  X256<F64> mxd = Max(ad, bd);
  for (int i = 0; i < 4; i++)
  {
    MTL_EQUAL_FLOAT(mnd.pData()[i], MTL::Min(a64[i], b64[i]), 0.0);
    MTL_EQUAL_FLOAT(mxd.pData()[i], MTL::Max(a64[i], b64[i]), 0.0);
  }

  // X256(U8) broadcast constructor.
  X256<U8> u(U8(0xAB));
  for (int i = 0; i < X256<U8>::Increment; i++)
    MTL_EQUAL(u.pData()[i], U8(0xAB));
}
#endif
