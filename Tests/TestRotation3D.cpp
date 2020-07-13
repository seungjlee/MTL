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

#include <MTL/Tools/Test.h>
#include <MTL/AxisAngle.h>
#include <MTL/Random.h>

using namespace MTL;

TEST(TestRotation3D)
{
  static const double kTol = 1e-12;

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

TEST(TestAxisAngle)
{
  static const double kPointTol     = 1e-15;
  static const double kAxisAngleTol = 1e-11;
  static const double kIdentityTol  = 1e-14;

  enum
  {
    N = 100000
  };

  Random random;

  for (int i = 0; i < N; i++)
  {
    Rotation3D<F64> R1(random.GetNext(0.0, kTwoPi),
                       random.GetNext(0.0, kTwoPi),
                       random.GetNext(0.0, kTwoPi));

    Rotation3D<F64> R2(random.GetNext(0.0, kTwoPi),
                       random.GetNext(0.0, kTwoPi),
                       random.GetNext(0.0, kTwoPi));

    AxisAngle<F64> V1(R1);
    AxisAngle<F64> V2(R2);

    Point3D<F64>
      point(random.GetNext(-1.0, 1.0), random.GetNext(-1.0, 1.0), random.GetNext(-1.0, 1.0));

    Point3D<F64> pointR = R1 * R2 * point;
    Point3D<F64> pointV = R1 * R2 * point;

    MTL_EQUAL_FLOAT(pointV.x(), pointR.x(), kPointTol);
    MTL_EQUAL_FLOAT(pointV.y(), pointR.y(), kPointTol);
    MTL_EQUAL_FLOAT(pointV.z(), pointR.z(), kPointTol);
  }

  Timer timeR;
  Timer timeV;

  for (int i = 0; i < N; i++)
  {
    Rotation3D<F64> R(random.GetNext(0.0, kTwoPi),
                      random.GetNext(0.0, kTwoPi),
                      random.GetNext(0.0, kTwoPi));

    AxisAngle<F64> V(R);

    timeR.Start();
    Rotation3D<F64> RR = V.GetRotationMatrix();
    timeR.Stop();

    Rotation3D<F64> I(RR.MultiplyByTranspose());

    for (I32 row = 0; row < I.Rows(); row++)
    {
      for (I32 col = 0; col < I.Cols(); col++)
      {
        if (row == col)
          MTL_EQUAL_FLOAT(I[row][col], 1.0, kIdentityTol);
        else
          MTL_EQUAL_FLOAT(I[row][col], 0.0, kIdentityTol);
      }
    }

    timeV.Start();
    AxisAngle<F64> VV(RR);
    timeV.Stop();

    MTL_EQUAL_FLOAT(VV.Vector().x(), V.Vector().x(), kAxisAngleTol);
    MTL_EQUAL_FLOAT(VV.Vector().y(), V.Vector().y(), kAxisAngleTol);
    MTL_EQUAL_FLOAT(VV.Vector().z(), V.Vector().z(), kAxisAngleTol);
  }

  wprintf(L"  GetRotationMatrix time:     %9.3f msecs (%d times)\n", timeR.Milliseconds(), N);
  wprintf(L"  Convert to axis angle time: %9.3f msecs (%d times)\n", timeV.Milliseconds(), N);
}

TEST(RotateToAxis)
{
  static const double kTol = 1e-14;

  enum
  {
    N = 100000
  };

  Random random;

  for (int i = 0; i < N; i++)
  {
    Vector3D<F64> v(random.GetNext(-10.0, 10.0),
                    random.GetNext(-10.0, 10.0),
                    random.GetNext(-10.0, 10.0));

    Rotation3D<F64> Rx = Rotation3D<F64>::RotationToXAxis(v);
    Rotation3D<F64> Ry = Rotation3D<F64>::RotationToYAxis(v);
    Rotation3D<F64> Rz = Rotation3D<F64>::RotationToZAxis(v);

    Point3D<F64> ptX = Rx * v;
    MTL_EQUAL_FLOAT(ptX.y(), 0.0, kTol);
    MTL_EQUAL_FLOAT(ptX.z(), 0.0, kTol);

    Point3D<F64> ptY = Ry * v;
    MTL_EQUAL_FLOAT(ptY.x(), 0.0, kTol);
    MTL_EQUAL_FLOAT(ptY.z(), 0.0, kTol);

    Point3D<F64> ptZ = Rz * v;
    MTL_EQUAL_FLOAT(ptZ.x(), 0.0, kTol);
    MTL_EQUAL_FLOAT(ptZ.y(), 0.0, kTol);

    Rotation3D<F64> Ix = Rx * Rotation3D<F64>::RotationFromXAxis(v);
    Rotation3D<F64> Iy = Ry * Rotation3D<F64>::RotationFromYAxis(v);
    Rotation3D<F64> Iz = Rz * Rotation3D<F64>::RotationFromZAxis(v);

    for (I32 row = 0; row < Ix.Rows(); row++)
    {
      for (I32 col = 0; col < Ix.Cols(); col++)
      {
        if (row == col)
        {
          MTL_EQUAL_FLOAT(Ix[row][col], 1.0, kTol);
          MTL_EQUAL_FLOAT(Iy[row][col], 1.0, kTol);
          MTL_EQUAL_FLOAT(Iz[row][col], 1.0, kTol);
        }
        else
        {
          MTL_EQUAL_FLOAT(Ix[row][col], 0.0, kTol);
          MTL_EQUAL_FLOAT(Iy[row][col], 0.0, kTol);
          MTL_EQUAL_FLOAT(Iz[row][col], 0.0, kTol);
        }
      }
    }
  }
}

static double Normalize(double angle)
{
  while (angle < -kPi)
    angle += 2*kPi;
  while (angle >= kPi)
    angle -= 2*kPi;

  return angle;
}

TEST(ComputeAngles)
{
  static const double kTol = 1e-14;

  enum
  {
    N = 100000
  };

  Random random;

  for (int i = 0; i < N; i++)
  {
    double angleX = random.GetNext(-1.0, 1.0) * 1.57;
    double angleY = random.GetNext(-1.0, 1.0) * 1.57;
    double angleZ = random.GetNext(-1.0, 1.0) * 1.57;

    Rotation3D<double> R(angleX, angleY, angleZ);

    Vector3D<double> angles = R.ComputeAngles();
    
    MTL_EQUAL_FLOAT(Normalize(angles.x()), Normalize(angleX), kTol);
    MTL_EQUAL_FLOAT(Normalize(angles.y()), Normalize(angleY), kTol);
    MTL_EQUAL_FLOAT(Normalize(angles.z()), Normalize(angleZ), kTol);
  }

  for (int i = 0; i < N; i++)
  {
    double angleX = random.GetNext(-1.0, 1.0) * 10;
    double angleY = random.GetNext(-1.0, 1.0) * 1.57;
    double angleZ = random.GetNext(-1.0, 1.0) * 10;

    Rotation3D<double> R(angleX, angleY, angleZ);

    Vector3D<double> angles = R.ComputeAngles();
    
    MTL_EQUAL_FLOAT(Normalize(angles.x()), Normalize(angleX), kTol);
    MTL_EQUAL_FLOAT(Normalize(angles.y()), Normalize(angleY), kTol);
    MTL_EQUAL_FLOAT(Normalize(angles.z()), Normalize(angleZ), kTol);
  }
}