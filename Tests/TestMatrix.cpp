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
#include <MTL/Matrix.h>

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
}
