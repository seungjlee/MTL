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
#include <MTL/Math/SVD.h>
#include <MTL/Math/Random.h>
#include <MTL/Timer.h>

using namespace MTL;

static const double kTol = 1e-14;

TEST(TestSVD_FixedMatrix)
{
  double a[4][4] = {{1, 5, 0, 0},
                    {0, 2, 6, 0},
                    {0, 0, 3, 7},
                    {0, 0, 0, 4}};

  SquareMatrix4x4 A(a);
  SquareMatrix4x4 V;
  double S[4];

  JacobiSVD<4,4>(A, S, V);

  wprintf(L"  %20.15f  %20.15f  %20.15f  %20.15f\n", S[0], S[1], S[2], S[3]);
  MTL_EQUAL_FLOAT(S[0], 8.895008771746831, kTol);
  MTL_EQUAL_FLOAT(S[1], 6.357395271827689, kTol);
  MTL_EQUAL_FLOAT(S[2], 4.522558768895992, kTol);
  MTL_EQUAL_FLOAT(S[3], 0.093842901552412, kTol);
}

TEST(TestSVD)
{
  double a[5][5] = {{ 1.0, 0.1, 0.9, 0.3, 0.6 },
                    { 0.1, 2.0, 0.7, 0.2, 0.4 },
                    { 0.9, 0.7, 3.0, 0.5, 0.8 },
                    { 0.3, 0.2, 0.5, 4.0, 0.3 },
                    { 0.6, 0.4, 0.8, 0.3, 5.0 }};

  DynamicMatrix<F64> A = Matrix<5,5,F64>(a);
  
  DynamicMatrix<F64> Ut1 = A;
  DynamicVector<F64> D1;
  DynamicMatrix<F64> Vt1;

  bool converged = JacobiSVDTransposed(Ut1, D1, Vt1);
  MTL_VERIFY(converged);
  
  for (int row = 0; row < Ut1.Rows(); row++)
    for (int col = 0; col < Ut1.Cols(); col++)
      MTL_EQUAL_FLOAT(Ut1[row][col], Vt1[row][col], kTol);

  for (int i = 0; i < (int)D1.Size() - 1; i++)
    MTL_LESS_THAN(D1[i+1], D1[i]);

  DynamicMatrix<F64> S(A.Rows(), A.Cols(), DynamicMatrix<F64>::eZeros);

  S.SetDiagonals(D1);
  DynamicMatrix<F64> A1 = Ut1.ComputeTranspose() * S * Vt1;
  for (int row = 0; row < A.Rows(); row++)
    for (int col = 0; col < A.Cols(); col++)
      MTL_EQUAL_FLOAT(A1[row][col], A[row][col], kTol);
  
  Ut1 = A;
  converged = JacobiSVDTransposed(Ut1, D1, Vt1, true);
  MTL_VERIFY(converged);
  
  for (int row = 0; row < Ut1.Rows(); row++)
    for (int col = 0; col < Ut1.Cols(); col++)
      MTL_EQUAL_FLOAT(Ut1[row][col], Vt1[row][col], kTol);

  for (int i = 0; i < (int)D1.Size() - 1; i++)
    MTL_LESS_THAN(D1[i], D1[i+1]);

  S.SetDiagonals(D1);
  A1 = Vt1.ComputeTranspose() * S * Vt1;
  for (int row = 0; row < A.Rows(); row++)
    for (int col = 0; col < A.Cols(); col++)
      MTL_EQUAL_FLOAT(A1[row][col], A[row][col], kTol);

  DynamicMatrix<F64> Ut2 = A;
  DynamicVector<F64> D2;
  converged = JacobiEigen(Ut2, D2);

  for (int i = 0; i < (int)D1.Size(); i++)
    MTL_EQUAL_FLOAT(D2[i], D1[i], kTol);

  for (int row = 0; row < Ut1.Rows(); row++)
    for (int col = 0; col < Ut1.Cols(); col++)
      MTL_EQUAL_FLOAT(Ut2[row][col], Ut1[row][col], kTol);

  DynamicMatrix<F64> S2(A.Rows(), A.Cols(), DynamicMatrix<F64>::eZeros);

  S2.SetDiagonals(D2);
  DynamicMatrix<F64> A2 = Ut2.ComputeTranspose() * S2 * Ut2;
  for (int row = 0; row < A.Rows(); row++)
    for (int col = 0; col < A.Cols(); col++)
      MTL_EQUAL_FLOAT(A2[row][col], A[row][col], kTol);
}

TEST(TestHouseholderJacobiSVD_FixedMatrix)
{
  static const double kHJTol = 1e-12;

  double a[4][4] = {{1, 5, 0, 0},
                    {0, 2, 6, 0},
                    {0, 0, 3, 7},
                    {0, 0, 0, 4}};

  // Reference singular values via plain JacobiSVD.
  SquareMatrix4x4 Aref(a);
  SquareMatrix4x4 Vref;
  double Sref[4];
  JacobiSVD<4,4>(Aref, Sref, Vref);

  SquareMatrix4x4 A(a);
  SquareMatrix4x4 V;
  double S[4];

  bool converged = HouseholderJacobiSVD<4,4>(A, S, V);
  MTL_VERIFY(converged);

  wprintf(L"  %20.15f  %20.15f  %20.15f  %20.15f\n", S[0], S[1], S[2], S[3]);
  for (int i = 0; i < 4; i++)
    MTL_EQUAL_FLOAT(S[i], Sref[i], kHJTol);

  // Reconstruct A: U * diag(S) * V^T should equal original.
  SquareMatrix4x4 Aorig(a);
  for (int row = 0; row < 4; row++)
  {
    for (int col = 0; col < 4; col++)
    {
      double sum = 0;
      for (int k = 0; k < 4; k++)
        sum += A[row][k] * S[k] * V[col][k];
      MTL_EQUAL_FLOAT(sum, Aorig[row][col], kHJTol);
    }
  }
}

TEST(TestHouseholderJacobiSVD_Tall)
{
  static const double kHJTol = 1e-11;

  // 6x4 tall matrix.
  double a[6][4] = {{ 1.0,  2.0, -1.0,  3.0},
                    {-2.0,  0.5,  4.0,  1.0},
                    { 0.7, -1.3,  2.2, -0.4},
                    { 3.1,  2.5,  0.8,  1.7},
                    {-0.6,  1.9, -2.1,  0.9},
                    { 1.4, -0.8,  3.3,  2.0}};

  Matrix<6,4,double> A0(a);

  // Reference using JacobiSVD.
  Matrix<6,4,double> Aref = A0;
  SquareMatrix<4,double> Vref;
  double Sref[4];
  JacobiSVD<6,4>(Aref, Sref, Vref);

  Matrix<6,4,double> A = A0;
  SquareMatrix<4,double> V;
  double S[4];

  bool converged = HouseholderJacobiSVD<6,4>(A, S, V);
  MTL_VERIFY(converged);

  for (int i = 0; i < 4; i++)
    MTL_EQUAL_FLOAT(S[i], Sref[i], kHJTol);

  // Verify orthonormality of U columns.
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      double dot = 0;
      for (int r = 0; r < 6; r++)
        dot += A[r][i] * A[r][j];
      MTL_EQUAL_FLOAT(dot, (i == j ? 1.0 : 0.0), kHJTol);
    }
  }

  // Verify orthonormality of V columns.
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      double dot = 0;
      for (int r = 0; r < 4; r++)
        dot += V[r][i] * V[r][j];
      MTL_EQUAL_FLOAT(dot, (i == j ? 1.0 : 0.0), kHJTol);
    }
  }

  // Verify reconstruction: U * diag(S) * V^T == A0.
  for (int row = 0; row < 6; row++)
  {
    for (int col = 0; col < 4; col++)
    {
      double sum = 0;
      for (int k = 0; k < 4; k++)
        sum += A[row][k] * S[k] * V[col][k];
      MTL_EQUAL_FLOAT(sum, A0[row][col], kHJTol);
    }
  }
}

TEST(TestSolveHouseholderJacobiSVD)
{
  static const double kHJTol = 1e-10;

  // Overdetermined system that has an exact solution.
  double a[5][3] = {{ 2.0,  1.0,  0.0},
                    { 1.0,  3.0,  1.0},
                    { 0.0,  1.0,  4.0},
                    { 1.0,  0.0,  2.0},
                    { 2.0,  2.0,  1.0}};

  Matrix<5,3,double> A(a);
  ColumnVector<3,double> xTrue;
  xTrue[0] = 1.5;
  xTrue[1] = -2.0;
  xTrue[2] = 0.75;

  ColumnVector<5,double> b;
  for (int i = 0; i < 5; i++)
  {
    b[i] = 0;
    for (int j = 0; j < 3; j++)
      b[i] += A[i][j] * xTrue[j];
  }

  ColumnVector<3,double> x;
  I32 rank;
  double conditionNumber;
  bool converged = SolveHouseholderJacobiSVD<5,3>(A, x, rank, conditionNumber, b);
  MTL_VERIFY(converged);
  MTL_EQUAL(rank, 3);

  for (int i = 0; i < 3; i++)
    MTL_EQUAL_FLOAT(x[i], xTrue[i], kHJTol);
}

TEST(TestComputePseudoinverseHouseholderJacobiSVD)
{
  static const double kHJTol = 1e-11;

  double a[5][3] = {{ 2.0,  1.0,  0.0},
                    { 1.0,  3.0,  1.0},
                    { 0.0,  1.0,  4.0},
                    { 1.0,  0.0,  2.0},
                    { 2.0,  2.0,  1.0}};

  Matrix<5,3,double> A(a);

  Matrix<3,5,double> pinvRef = ComputePseudoinverseJacobiSVD(A);
  Matrix<3,5,double> pinv    = ComputePseudoinverseHouseholderJacobiSVD(A);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 5; j++)
      MTL_EQUAL_FLOAT(pinv[i][j], pinvRef[i][j], kHJTol);
}

TEST(TestHouseholderJacobiSVD_Dynamic)
{
  static const double kHJTol = 1e-10;

  // Tall A: 12 x 4 (M=12, N=4 columns).
  Random random;
  DynamicMatrix<F64> A0 = random.DynamicMatrix<F64>(12, 4, -1.0, 1.0);

  // Reference solve via SolveJacobiSVDTransposed.
  DynamicVector<F64> b = random.DynamicVector<F64>(12, -1.0, 1.0);
  DynamicMatrix<F64> At0 = A0.ComputeTranspose();

  DynamicVector<F64> xRef;
  I32 rankRef;
  F64 condRef;
  SolveJacobiSVDTransposed(xRef, rankRef, condRef, At0, b);

  DynamicVector<F64> xHJ;
  I32 rankHJ;
  F64 condHJ;
  bool converged = SolveHouseholderJacobiSVDTransposed(xHJ, rankHJ, condHJ, At0, b);
  MTL_VERIFY(converged);
  MTL_EQUAL(rankHJ, rankRef);

  for (int i = 0; i < 4; i++)
    MTL_EQUAL_FLOAT(xHJ[i], xRef[i], kHJTol);

  // Decomposition correctness: U * diag(D) * V^T == A.
  DynamicMatrix<F64> A = A0;
  DynamicMatrix<F64> V;
  DynamicVector<F64> D;
  HouseholderJacobiSVD(A, D, V);

  for (int row = 0; row < 12; row++)
  {
    for (int col = 0; col < 4; col++)
    {
      double sum = 0;
      for (int k = 0; k < 4; k++)
        sum += A[row][k] * D[k] * V[col][k];
      MTL_EQUAL_FLOAT(sum, A0[row][col], kHJTol);
    }
  }
}

template <I32 M, I32 N>
static void RunSVDPerformance(I32 repeats)
{
  static const double kPerfTol = 1e-9;

  Random random;

  // Pre-generate matrices so the timed loops only measure SVD work.
  std::vector<Matrix<M,N,double>> matrices(repeats);
  for (I32 i = 0; i < repeats; i++)
    matrices[i] = random.Matrix<M,N,double>(-1.0, 1.0);

  Timer tJ;
  Timer tHJ;

  double maxErrJ  = 0;
  double maxErrHJ = 0;
  I32    convFailJ  = 0;
  I32    convFailHJ = 0;

  // Plain Jacobi SVD.
  for (I32 r = 0; r < repeats; r++)
  {
    Matrix<M,N,double> A = matrices[r];
    SquareMatrix<N,double> V;
    double S[N];

    tJ.Start();
    bool converged = JacobiSVD<M,N>(A, S, V);
    tJ.Stop();

    if (!converged) convFailJ++;

    // Reconstruction error.
    double err = 0;
    for (I32 i = 0; i < M; i++)
    {
      for (I32 j = 0; j < N; j++)
      {
        double sum = 0;
        for (I32 k = 0; k < N; k++)
          sum += A[i][k] * S[k] * V[j][k];
        err = Max(err, Abs(sum - matrices[r][i][j]));
      }
    }
    maxErrJ = Max(maxErrJ, err);
  }

  // Householder bidiagonalization + Jacobi SVD.
  for (I32 r = 0; r < repeats; r++)
  {
    Matrix<M,N,double> A = matrices[r];
    SquareMatrix<N,double> V;
    double S[N];

    tHJ.Start();
    bool converged = HouseholderJacobiSVD<M,N>(A, S, V);
    tHJ.Stop();

    if (!converged) convFailHJ++;

    double err = 0;
    for (I32 i = 0; i < M; i++)
    {
      for (I32 j = 0; j < N; j++)
      {
        double sum = 0;
        for (I32 k = 0; k < N; k++)
          sum += A[i][k] * S[k] * V[j][k];
        err = Max(err, Abs(sum - matrices[r][i][j]));
      }
    }
    maxErrHJ = Max(maxErrHJ, err);
  }

  double speedup = (tHJ.Milliseconds() > 0) ? tJ.Milliseconds() / tHJ.Milliseconds() : 0;

  wprintf(L"  %dx%d  x %d:\n", M, N, repeats);
  wprintf(L"    JacobiSVD            : %9.3f msecs, max recon err = %.3e, conv-fail = %d\n",
          tJ.Milliseconds(), maxErrJ, convFailJ);
  wprintf(L"    HouseholderJacobiSVD : %9.3f msecs, max recon err = %.3e, conv-fail = %d\n",
          tHJ.Milliseconds(), maxErrHJ, convFailHJ);
  wprintf(L"    speedup (J / HJ)     : %9.3fx\n", speedup);

  MTL_LESS_THAN(maxErrJ,  kPerfTol);
  MTL_LESS_THAN(maxErrHJ, kPerfTol);
}

TEST(TestJacobiSVD_vs_HouseholderJacobiSVD_Performance)
{
  // Square: HouseholderJacobiSVD pays bidiag overhead with no Jacobi-size win.
  RunSVDPerformance< 4, 4>(20000);
  RunSVDPerformance< 8, 8>(5000);

  // Mildly rectangular.
  RunSVDPerformance<16, 8>(2000);

  // Tall: HouseholderJacobiSVD should start to win as M >> N.
  RunSVDPerformance<32, 6>(1000);
  RunSVDPerformance<64, 8>(500);
  RunSVDPerformance<128,8>(200);
  RunSVDPerformance<256,8>(100);
}

