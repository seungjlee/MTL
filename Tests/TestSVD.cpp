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
