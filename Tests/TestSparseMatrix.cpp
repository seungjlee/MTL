//
// Math Template Library
//
// Copyright (c) 2015-2017: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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
#include <MTL/Random.h>
#include <MTL/OptimizerLevenbergMarquardt.h>

using namespace MTL;

static const double kTol = 1e-11;

class SparseOptimizer : public SparseOptimizerLevenbergMarquardt<F64>
{
public:
  SparseOptimizer()  {}

protected:
  virtual void CostFunction(DynamicVector<F64>& residuals, const Parameters& parameters)
  {
  }

  virtual void ComputeJacobian(CompressedSparseMatrix<F64>& J, const Parameters& parameters)
  {
  }

  virtual void ComputeSparsityMatrix(I32 numberOfParams)
  {
  }
};

TEST(TestSparseOptimizer)
{
  SparseOptimizer optimizer;
}

TEST(TestMultiplication)
{
  F64 a[4][3] = {{ 1, 5, 0 },
                 { 2, 6, 0 },
                 { 3, 0, 7 },
                 { 4, 0, 8 }};

  F64 b[4] = { 3, -2, -1, 1 };

  DynamicMatrix<F64> fullA = Matrix<4,3,F64>(a);
  DynamicVector<F64> B(b, 4);

  CompressedSparseMatrix<F64> A(fullA);

  MTL_EQUAL(A.NumberOfElements(), 8);

  DynamicVector<F64> C;
  A.MultiplyTransposed(C, B);
  DynamicVector<F64> D = fullA.ComputeTranspose() * B;

  MTL_EQUAL(D.Size(), C.Size());

  for (U32 i = 0; i < C.Size(); i++)
    MTL_EQUAL_FLOAT(D[i], C[i], kTol);

  DynamicMatrix<F64> M2;
  A.MultiplyTransposeByThis(M2);
  DynamicMatrix<F64> M1 = fullA.ComputeTranspose() * fullA;

  for (I32 row = 0; row < M1.Rows(); row++)
    for (I32 col = 0; col < M1.Cols(); col++)
      MTL_EQUAL_FLOAT(M2[row][col], M1[row][col], kTol);

  CompressedSparseMatrix<F64> M3;
  A.MultiplyTransposeByThis(M3);

  for (I32 row = 0; row < M1.Rows(); row++)
    for (I32 col = 0; col < M1.Cols(); col++)
      MTL_EQUAL_FLOAT(M3(row,col), M1[row][col], kTol);
}

TEST(TestLargeMultiplication)
{
  enum
  {
    kRepeats = 10,
    kRows = 20000,
    kCols = 200
  };

  Timer timeSparseMtxV;
  Timer timeFullMtxV;
  Timer timeSparseMtxM;
  Timer timeFullMtxM;
  Timer timeOptimizeMtxM;
  Timer timeOptimizedSequentialMtxM;
  Timer timeOptimizedParallelMtxM;

  Random random;
  CompressedSparseMatrix<F64> M3, M4, M5;


  for (I32 i = 0; i < kRepeats; i++)
  {
    DynamicVector<F64> B = random.DynamicVector<F64>(kRows, -1, 1);
    DynamicMatrix<F64> fullA = random.DynamicMatrix<F64>(kRows, kCols, -1, 1);

    for (I32 row = 0; row < fullA.Rows(); row++)
      for (I32 col = 0; col < fullA.Cols(); col++)
        if (Abs(fullA[row][col]) < 0.9)
          fullA[row][col] = 0.0;

    DynamicMatrix<F64> fullAt = fullA.ComputeTranspose();

    CompressedSparseMatrix<F64> A(fullA);

    DynamicVector<F64> C;

    timeSparseMtxV.Start();
    A.MultiplyTransposed(C, B);
    timeSparseMtxV.Stop();

    timeFullMtxV.Start();
    DynamicVector<F64> D = fullAt * B;
    timeFullMtxV.Stop();

    MTL_EQUAL(D.Size(), C.Size());

    for (U32 i = 0; i < C.Size(); i++)
      MTL_EQUAL_FLOAT(D[i], C[i], kTol);

    DynamicMatrix<F64> M2;
    A.MultiplyTransposeByThis(M2);

    timeFullMtxM.Start();
    DynamicMatrix<F64> M1 = fullAt.MultiplyByTranspose();
    timeFullMtxM.Stop();

    for (I32 row = 0; row < M1.Rows(); row++)
      for (I32 col = 0; col < M1.Cols(); col++)
        MTL_EQUAL_FLOAT(M2[row][col], M1[row][col], kTol);

    timeSparseMtxM.Start();
    A.MultiplyTransposeByThisParallel(M3);
    timeSparseMtxM.Stop();

    for (I32 row = 0; row < M1.Rows(); row++)
      for (I32 col = 0; col < M1.Cols(); col++)
        MTL_EQUAL_FLOAT(M3(row,col), M1[row][col], kTol);

    timeOptimizeMtxM.Start();
    A.OptimizeMultiplyTransposeByThis();
    timeOptimizeMtxM.Stop();

    timeOptimizedSequentialMtxM.Start();
    A.MultiplyTransposeByThis(M4);
    timeOptimizedSequentialMtxM.Stop();

    for (I32 row = 0; row < M1.Rows(); row++)
      for (I32 col = 0; col < M1.Cols(); col++)
        MTL_EQUAL_FLOAT(M4(row,col), M1[row][col], kTol);

    timeOptimizedParallelMtxM.Start();
    A.MultiplyTransposeByThisParallel(M5);
    timeOptimizedParallelMtxM.Stop();

    for (I32 row = 0; row < M1.Rows(); row++)
      for (I32 col = 0; col < M1.Cols(); col++)
        MTL_EQUAL_FLOAT(M5(row,col), M1[row][col], kTol);
  }

  printf("\n");
  printf("Dense  Mt x V time: %9.3f msecs.\n", timeFullMtxV.Milliseconds());
  printf("Sparse Mt x V time: %9.3f msecs.\n", timeSparseMtxV.Milliseconds());
  printf("Dense  Mt x M time: %9.3f msecs.\n", timeFullMtxM.Milliseconds());
  printf("Sparse Mt x M time: %9.3f msecs.\n\n", timeSparseMtxM.Milliseconds());
  printf("Prepare optimized sparse Mt x M time:      %9.3f msecs.\n",
         timeOptimizeMtxM.Milliseconds());
  printf("Optimized sparse Mt x M time (sequential): %9.3f msecs.\n",
         timeOptimizedSequentialMtxM.Milliseconds());
  printf("Optimized sparse Mt x M time (parallel):   %9.3f msecs.\n",
         timeOptimizedParallelMtxM.Milliseconds());
  printf("\n");
}
