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

#include <MTL/Test.h>
#include <MTL/LinearAlgebra.h>
#include <MTL/LDLt.h>
#include <MTL/QR.h>
#include <MTL/SVD.h>
#include <MTL/Random.h>
#include <MTL/ProjectionToImageTransform.h>

using namespace MTL;

static const double kTol = 1e-12;

TEST(TestMatrix)
{
  double a[][3] = {{ 4, -3,  2},
                   { 2,  4, -5},
                   {-2,  1,  6}};
  double b[] = {1, -4,  3};
  double solution[] = {-0.211267605633803, -0.295774647887324, 0.478873239436620};

  SquareMatrix3x3 A(a);
  ColumnVector3D x(b);

  A.SolveLUP(x);

  for (int i = 0; i < x.Size(); i++)
    MTL_EQUAL_FLOAT(x[i], solution[i], kTol);

  ColumnVector3D y = A.Inverse() * ColumnVector3D(b);

  for (int i = 0; i < x.Size(); i++)
    MTL_EQUAL_FLOAT(y[i], solution[i], kTol);
}

TEST(TestMatrixConcatenation)
{
  double a[][2] = {{ 1, 2 },
                   { 4, 5 }};

  double col[] = { 3, 6 };
  double row[] = { 7, 8, 9 };

  SquareMatrix2x2 A(a);
  ColumnVector2D C(col);
  RowVector3D R(row);

  SquareMatrix3x3 M = (A || C) && R;

  for (int row = 0; row < M.Rows(); row++)
    for (int col = 0; col < M.Cols(); col++)
      MTL_EQUAL_FLOAT(M[row][col], row * 3.0 + col + 1.0, kTol);
}

TEST(TestNullSpace)
{
  double v[] = { 0, 0, 1, 1 };

  RowVector4D V(v);
  Matrix<4,3> Null = NullSpace(V);

  for (int row = 0; row < 4; row++)
    wprintf(L"  %18.15f %18.15f %18.15f\n", Null[row][0], Null[row][1], Null[row][2]);

  RowVector3D zeros = V * Null;

  for (int i = 0; i < zeros.Cols(); i++)
    MTL_EQUAL_FLOAT(zeros[i], 0.0, kTol);
}

TEST(TestTranspose)
{
  double a[][3] = {{ 4, -3,  2},
                   { 2,  4, -5},
                   {-2,  1,  6}};

  SquareMatrix3x3 At(a);
  At.Transpose();

  for (int row = 0; row < At.Rows(); row++)
    for (int col = 0; col < At.Cols(); col++)
      MTL_EQUAL_FLOAT(At[row][col], a[col][row], kTol);

  SquareMatrix3x3 A(a);

  SquareMatrix3x3 AAt = A*At;
  SquareMatrix3x3 AtA = At*A;

  SquareMatrix3x3 M1 = A.MultiplyByTranspose();
  SquareMatrix3x3 M2 = A.MultiplyTransposeByThis();

  for (int row = 0; row < M1.Rows(); row++)
    for (int col = 0; col < M1.Cols(); col++)
      MTL_EQUAL_FLOAT(M1[row][col], AAt[row][col], kTol);

  for (int row = 0; row < M2.Rows(); row++)
    for (int col = 0; col < M2.Cols(); col++)
      MTL_EQUAL_FLOAT(M2[row][col], AtA[row][col], kTol);
}

TEST(TestInverseAndDeterminant)
{
  double a[][3] = {{ 4, -3,  2},
                   { 2,  4, -5},
                   {-2,  1,  6}};

  SquareMatrix3x3 A(a);
  SquareMatrix3x3 InvA1 = A.Inverse();
  SquareMatrix3x3 InvA2 = InverseRecursive(A);

  for (int row = 0; row < InvA1.Rows(); row++)
    for (int col = 0; col < InvA1.Cols(); col++)
      MTL_EQUAL_FLOAT(InvA1[row][col], InvA2[row][col], kTol);

  MTL_EQUAL_FLOAT(A.Determinant(), DeterminantRecursive(A), kTol);
  MTL_EQUAL_FLOAT(A.Determinant(), A.StableDeterminant(), kTol);

  I32 rank;
  F64 conditionNumber;
  ColumnVector3D xSVD;
  SquareMatrix3x3 temp(A);
  SolveJacobiSVD<3>(temp, xSVD, rank, conditionNumber);
  wprintf(L"SVD Condition number of A = %.15f\n", conditionNumber);
  wprintf(L"LUP Condition number of A = %.15f\n", A.ConditionNumber());

  ProjectionToImageTransform<F64> P1;
  ProjectionToImageTransform<F64> P2;

#ifdef WIN32
  SquareMatrix3x3 F = ComputeFundamentalMatrix(P1,P2);
#endif
}

TEST(TestPseudoinverse)
{
  static const double kTol = 1e-14;

  enum
  {
    kRepeats = 1000
  };

  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    Matrix<5,3> A = random.Matrix<5,3,F64>(-1, 1);

    for (I32 col = 0; col < A.Cols(); col++)
      A[col][col] += 1.0;
    
    Matrix<3,5> pinvA = ComputePseudoinverseJacobiSVD(A);
    Matrix<3,3> I = pinvA * A;

    for (I32 row = 0; row < I.Rows(); row++)
    {
      for (I32 col = 0; col < I.Cols(); col++)
      {
        if (row == col)
          MTL_EQUAL_FLOAT(I[row][col], 1.0, kTol);
        else
          MTL_EQUAL_FLOAT(I[row][col], 0.0, kTol);
      }
    }
  }
}

template<I32 N, class T>
void TestSolvers(const T& tol)
{
  enum
  {
    kRepeats = 10000
  };

  Timer time_SVD;
  Timer time_QR;
  Timer time_LDLt;
  Timer time_LUP;
  Timer time_Inverse;
  T maxRMS_SVD = 0;
  T maxRMS_QR = 0;
  T maxRMS_LUP = 0;
  T maxRMS_LDLt = 0;
  T maxRMS_Inverse = 0;

  T maxConditionNumber = 0;

  ColumnVector<N,T> residuals;

  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    SquareMatrix<N,T> A = random.Matrix<N,N,T>(-1, 1);
    ColumnVector<N,T> b = random.Matrix<N,1,T>(-1, 1);;

    // Create a positive definite symmetric matrix that is well conditioned.
    A = A.MultiplyByTranspose();
    A += SquareMatrix<N,T>(SquareMatrix<N,T>::eIdentity) * 10;

    I32 rank;
    T conditionNumber;

    ColumnVector<N,T> xSVD = b;
    time_SVD.Start();
    SquareMatrix<N,T> temp(A);
    SolveJacobiSVD<N>(temp, xSVD, rank, conditionNumber, tol);
    time_SVD.Stop();
    residuals = A * xSVD - b;
    maxRMS_SVD = Max(maxRMS_SVD, residuals.RMS());

    maxConditionNumber = Max(maxConditionNumber, conditionNumber);

    ColumnVector<N,T> xQR = b;
    time_QR.Start();
    temp = A;
    SolveHouseholderQR<N>(xQR, temp);
    time_QR.Stop();
    residuals = A * xQR - b;
    maxRMS_QR = Max(maxRMS_QR, residuals.RMS());

    ColumnVector<N,T> xLDLt = b;
    time_LDLt.Start();
    temp = A;
    SolveLDLt<N>(xLDLt, temp, tol);
    time_LDLt.Stop();
    residuals = A * xLDLt - b;
    maxRMS_LDLt = Max(maxRMS_LDLt, residuals.RMS());

    ColumnVector<N,T> xLUP = b;
    time_LUP.Start();
    A.SolveLUP(xLUP);
    time_LUP.Stop();
    residuals = A * xLUP - b;
    maxRMS_LUP = Max(maxRMS_LUP, residuals.RMS());

    time_Inverse.Start();
    ColumnVector<N,T> xInverse = A.Inverse() * b;
    time_Inverse.Stop();
    residuals = A * xInverse - b;
    maxRMS_Inverse = Max(maxRMS_Inverse, residuals.RMS());
  }

  wprintf(L"  Maximum condition number was: %f\n", maxConditionNumber);
  wprintf(L"  Every solver ran %d times (%dx%d matrices).\n", kRepeats, N, N);

  wprintf(L"  Solve SVD:     %9.3f msecs, Max RMS = %e\n",
          time_SVD.Milliseconds(), maxRMS_SVD);
  wprintf(L"  Solve QR:      %9.3f msecs, Max RMS = %e\n",
          time_QR.Milliseconds(), maxRMS_QR);
  wprintf(L"  Solve LDLt:    %9.3f msecs, Max RMS = %e\n",
          time_LDLt.Milliseconds(), maxRMS_LDLt);
  wprintf(L"  Solve LUP:     %9.3f msecs, Max RMS = %e\n",
          time_LUP.Milliseconds(), maxRMS_LUP);
  wprintf(L"  Solve Inverse: %9.3f msecs, Max RMS = %e\n\n",
          time_Inverse.Milliseconds(), maxRMS_Inverse);

  MTL_LESS_THAN(maxRMS_SVD, tol);
  MTL_LESS_THAN(maxRMS_LDLt, tol);
  MTL_LESS_THAN(maxRMS_LUP, tol);
  MTL_LESS_THAN(maxRMS_Inverse, tol);
}

TEST(TestSolversF32)
{
  static const float kTol = 1e-6f;

  TestSolvers< 3,F32>(kTol);
  TestSolvers< 4,F32>(kTol);
  TestSolvers< 5,F32>(kTol);
  TestSolvers< 6,F32>(kTol);
  TestSolvers< 9,F32>(kTol);
  //TestSolvers<11,F32>(kTol);
}

TEST(TestSolversF64)
{
  static const double kTol = 1e-14;

  TestSolvers< 3,F64>(kTol);
  TestSolvers< 4,F64>(kTol);
  TestSolvers< 5,F64>(kTol);
  TestSolvers< 6,F64>(kTol);
  TestSolvers< 9,F64>(kTol);
  //TestSolvers<11,F64>(kTol);
}

TEST(Test_SubColumn)
{
  static const double kTol = 1e-14;

  ColumnVector<9> v;
  
  for (int i = 0; i < v.Size(); i++)
    v[i] = i;

  SquareMatrix3x3 M = v.SubColumn<0,3>() || v.SubColumn<3,3>() || v.SubColumn<6,3>();

  double expected[3][3] = {{ 0, 3, 6 },
                           { 1, 4, 7 },
                           { 2, 5, 8 }};

  for (int row = 0; row < M.Rows(); row++)
    for (int col = 0; col < M.Cols(); col++)
      MTL_EQUAL_FLOAT(M[row][col], expected[row][col], kTol);
}

TEST(Test_ColumnVector_MaxNorm)
{
  double v1[] = { -1, 5, 3, -2 };
  double v2[] = { -7, 5, 3, -2 };
  double v3[] = { -7, 5, 8, -2 };
  double v4[] = { -7, 5, 8, -9 };
  ColumnVector4D V1(v1);
  ColumnVector4D V2(v2);
  ColumnVector4D V3(v3);
  ColumnVector4D V4(v4);
  MTL_EQUAL(V1.MaxNorm(), 5);
  MTL_EQUAL(V2.MaxNorm(), 7);
  MTL_EQUAL(V3.MaxNorm(), 8);
  MTL_EQUAL(V4.MaxNorm(), 9);
}

TEST(Test_SetAll)
{
  static const double kTol = 1e-14;

  enum
  {
    kRepeats = 10000
  };

  SquareMatrix3x3 M;
  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    double val = random.GetNext(-11.0, 11.0);
    M.SetAll(val);

    for (int row = 0; row < M.Rows(); row++)
      for (int col = 0; col < M.Cols(); col++)
        MTL_EQUAL_FLOAT(M[row][col], val, kTol);
  }
}

TEST(Test_SetAll_SetColumn)
{
  static const double kTol = 1e-14;

  enum
  {
    kRepeats = 10000
  };

  SquareMatrix3x3 M;
  Random random;

  for (I32 i = 0; i < kRepeats; i++)
  {
    double val = random.GetNext(-11.0, 11.0);
    M.SetAll(val);

    for (int row = 0; row < M.Rows(); row++)
      for (int col = 0; col < M.Cols(); col++)
        MTL_EQUAL_FLOAT(M[row][col], val, kTol);

    ColumnVector3D cols[3];
    cols[0] = ColumnVector3D(random.Matrix<3,1,F64>(-11.0, 11.0));
    cols[1] = ColumnVector3D(random.Matrix<3,1,F64>(-11.0, 11.0));
    cols[2] = ColumnVector3D(random.Matrix<3,1,F64>(-11.0, 11.0));

    for (int col = 0; col < M.Cols(); col++)
      M.SetColumn(col, cols[col]);

    for (int row = 0; row < M.Rows(); row++)
      for (int col = 0; col < M.Cols(); col++)
        MTL_EQUAL_FLOAT(M[row][col], cols[col][row], kTol);
  }
}
