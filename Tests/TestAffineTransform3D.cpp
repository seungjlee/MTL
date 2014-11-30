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

#include <MTL/Test.h>
#include <MTL/AxisAngle.h>
#include <MTL/RotationTranslation3D.h>

using namespace MTL;

static const double kTol = 1e-12;

TEST(TestRotation3D)
{
  // Some simple tests.
  Rotation3D<F64> Rx = Rotation3D<F64>::RotationX(90 * kDegreesToRadians);
  Rotation3D<F64> Ry = Rotation3D<F64>::RotationY(90 * kDegreesToRadians);
  Rotation3D<F64> Rz = Rotation3D<F64>::RotationZ(90 * kDegreesToRadians);

  Point3D<F64> v1, v2, v3;

  v1 = Rx * Point3D<F64>(1,0,0);
  v2 = Ry * Point3D<F64>(1,0,0);
  v3 = Rz * Point3D<F64>(1,0,0);
  MTL_EQUAL_FLOAT(v1.x(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v2.z(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v3.y(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);

  v1 = Rx * Point3D<F64>(0,1,0);
  v2 = Ry * Point3D<F64>(0,1,0);
  v3 = Rz * Point3D<F64>(0,1,0);
  MTL_EQUAL_FLOAT(v1.z(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v2.y(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v3.x(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);

  v1 = Rx * Point3D<F64>(0,0,1);
  v2 = Ry * Point3D<F64>(0,0,1);
  v3 = Rz * Point3D<F64>(0,0,1);
  MTL_EQUAL_FLOAT(v1.y(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v2.x(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v3.z(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);

  AxisAngle<F64> Ax(Rx);
  AxisAngle<F64> Ay(Ry);
  AxisAngle<F64> Az(Rz);

  v1 = Ax * Point3D<F64>(1,0,0);
  v2 = Ay * Point3D<F64>(1,0,0);
  v3 = Az * Point3D<F64>(1,0,0);
  MTL_EQUAL_FLOAT(v1.x(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v2.z(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v3.y(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);

  v1 = Ax * Point3D<F64>(0,1,0);
  v2 = Ay * Point3D<F64>(0,1,0);
  v3 = Az * Point3D<F64>(0,1,0);
  MTL_EQUAL_FLOAT(v1.z(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v2.y(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v3.x(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);

  v1 = Ax * Point3D<F64>(0,0,1);
  v2 = Ay * Point3D<F64>(0,0,1);
  v3 = Az * Point3D<F64>(0,0,1);
  MTL_EQUAL_FLOAT(v1.y(),   -1.0, kTol);
  MTL_EQUAL_FLOAT(v1.Sum(), -1.0, kTol);
  MTL_EQUAL_FLOAT(v2.x(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v2.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v3.z(),    1.0, kTol);
  MTL_EQUAL_FLOAT(v3.Sum(),  1.0, kTol);
  MTL_EQUAL_FLOAT(v1.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v2.SumOfSquares(), 1.0, kTol);
  MTL_EQUAL_FLOAT(v3.SumOfSquares(), 1.0, kTol);
}

TEST(TestRotationTranslation3D)
{
  RotationTranslation3D<F64> T;
}
