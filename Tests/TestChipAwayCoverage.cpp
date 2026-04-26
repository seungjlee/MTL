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

// Coverage chip-aways for cold paths in DynamicVector, CompressedSparseMatrix,
// and WorkerThread.

#include <MTL/CPU.h>
#include <MTL/Math/AxisAngle.h>
#include <MTL/Math/DynamicVector.h>
#include <MTL/Math/Matrix.h>
#include <MTL/Math/SparseMatrix.h>
#include <MTL/Tools/Test.h>
#include <MTL/Tools/WorkerThread.h>

using namespace MTL;

// AssignAll has both a parallel-OpenMP path and a single-thread fallback.
// The default fixture has multiple threads, so force one explicitly to
// drive the sequential branch in OptimizedAssignAll/AssignAll.
TEST(ChipAway_DynamicVector_AssignAll_SingleThread)
{
  U64 saved = CPU::Instance().NumberOfThreads();
  CPU::Instance().NumberOfThreads(1);

  // Resize(newSize, value) -> AssignAll path on a non-trivial size.
  DynamicVector<F64> v;
  v.Resize(257, 7.5);
  for (SizeType i = 0; i < v.Size(); i++)
    MTL_EQUAL_FLOAT(v[i], 7.5, 0.0);

  // Larger size to drive the SIMD body in AssignAll_Stream.
  DynamicVector<F32> w;
  w.Resize(4096, -0.5f);
  for (SizeType i = 0; i < w.Size(); i++)
    MTL_EQUAL_FLOAT(w[i], -0.5f, 0.0);

  // Resize(newSize, value) on an already-large buffer drives the in-place
  // assign branch (DynamicVector.h lines 176-177).
  DynamicVector<F64> reused;
  reused.Resize(1024);
  reused.Resize(64, 9.25);
  for (SizeType i = 0; i < reused.Size(); i++)
    MTL_EQUAL_FLOAT(reused[i], 9.25, 0.0);

  CPU::Instance().NumberOfThreads(saved);
}

// SequentialConvert is the elementwise type-converting copy used by the
// constructor that takes a DynamicVector of a different element type.
TEST(ChipAway_DynamicVector_SequentialConvert)
{
  DynamicVector<I32> src(33);
  for (SizeType i = 0; i < src.Size(); i++)
    src[i] = (I32)i - 10;
  DynamicVector<F64> dst(src);
  for (SizeType i = 0; i < dst.Size(); i++)
    MTL_EQUAL_FLOAT(dst[i], (F64)src[i], 0.0);
}

// Construct two sparse columns whose row indices interleave in both
// directions ("ai[p] < ai[q]" and "ai[p] > ai[q]" branches), and one column
// missing its diagonal entry to drive MaxDiagonal's "found == false" case.
// Then call MultiplyTransposeByThis (the non-Parallel path that the LM
// tests don't reach) to drive the row-pair scan branches in SparseMatrix.h.
TEST(ChipAway_SparseMatrix_MultiplyTransposeByThis_NonParallel)
{
  CompressedSparseMatrix<F64> A;

  enum { Rows = 6, Cols = 3 };
  DynamicVector<DynamicVector<I32>> sparseColumns(Cols);
  DynamicVector<DynamicVector<F64>> sparseValues(Cols);

  // Column 0: rows {0, 2, 4}. Includes diagonal at row 0.
  sparseColumns[0].PushBack(0);  sparseValues[0].PushBack(2.0);
  sparseColumns[0].PushBack(2);  sparseValues[0].PushBack(3.0);
  sparseColumns[0].PushBack(4);  sparseValues[0].PushBack(5.0);

  // Column 1: rows {1, 2, 5}. Includes diagonal at row 1. Interleaves
  // with column 0 so the row-merge has both ai[p]<ai[q] and ai[p]>ai[q].
  sparseColumns[1].PushBack(1);  sparseValues[1].PushBack(7.0);
  sparseColumns[1].PushBack(2);  sparseValues[1].PushBack(11.0);
  sparseColumns[1].PushBack(5);  sparseValues[1].PushBack(13.0);

  // Column 2: rows {3, 4, 5}. NO diagonal entry at row 2 -> exercises
  // MaxDiagonal's "found == false" branch.
  sparseColumns[2].PushBack(3);  sparseValues[2].PushBack(17.0);
  sparseColumns[2].PushBack(4);  sparseValues[2].PushBack(19.0);
  sparseColumns[2].PushBack(5);  sparseValues[2].PushBack(23.0);

  // updateOptimizedMultiplyStructure=false skips the symbolic prepass so
  // MultiplyTransposeByThis takes the unoptimised row-merge path.
  A.Create(Rows, Cols, sparseColumns, sparseValues, /*update*/ false);

  DynamicMatrix<F64> P;
  A.MultiplyTransposeByThis(P);

  // Spot-check a few entries against the dense product A^T * A.
  // P[0][0] = 4 + 9 + 25 = 38
  MTL_EQUAL_FLOAT(P[0][0], 38.0, 1e-12);
  // P[1][1] = 49 + 121 + 169 = 339
  MTL_EQUAL_FLOAT(P[1][1], 339.0, 1e-12);
  // P[2][2] = 289 + 361 + 529 = 1179
  MTL_EQUAL_FLOAT(P[2][2], 1179.0, 1e-12);
  // P[0][1] = 3*11 = 33 (only row 2 overlaps).
  MTL_EQUAL_FLOAT(P[0][1], 33.0, 1e-12);
  MTL_EQUAL_FLOAT(P[1][0], 33.0, 1e-12);
  // P[0][2] = 5*19 = 95 (only row 4 overlaps).
  MTL_EQUAL_FLOAT(P[0][2], 95.0, 1e-12);
  // P[1][2] = 13*23 = 299 (only row 5 overlaps).
  MTL_EQUAL_FLOAT(P[1][2], 299.0, 1e-12);

  MTL_EQUAL_FLOAT(A.MaxDiagonal(), 7.0, 0.0);  // column 2 contributes 0.

  // Re-run with updateOptimizedMultiplyStructure=true so MultiplyTransposeByThis
  // takes the RowIndexPairs-driven branch (lines 609-645 in SparseMatrix.h).
  CompressedSparseMatrix<F64> Aopt;
  Aopt.Create(Rows, Cols, sparseColumns, sparseValues, /*update*/ true);
  DynamicMatrix<F64> Popt;
  Aopt.MultiplyTransposeByThis(Popt);
  MTL_EQUAL_FLOAT(Popt[0][0], 38.0, 1e-12);
  MTL_EQUAL_FLOAT(Popt[1][1], 339.0, 1e-12);
  MTL_EQUAL_FLOAT(Popt[2][2], 1179.0, 1e-12);
  MTL_EQUAL_FLOAT(Popt[0][1], 33.0, 1e-12);
  MTL_EQUAL_FLOAT(Popt[0][2], 95.0, 1e-12);
  MTL_EQUAL_FLOAT(Popt[1][2], 299.0, 1e-12);
}

// Drives WorkerThread's three HandleException overloads. The base class
// re-throws inside these handlers, which would terminate the worker thread
// and crash the test, so we override them to count and swallow.
template <class WorkerClass>
class ThrowingWorker : public WorkerClass
{
public:
  ThrowingWorker(int mode)
    : WorkerClass(std::wstring(L"throwing")), Mode_(mode), HandledMTL_(0),
      HandledStd_(0), HandledUnknown_(0)
  {
  }

  void ProcessWork(const std::vector<int>&) override
  {
    if (Mode_ == 0)
      throw MTL::Exception(L"intentional MTL::Exception");
    else if (Mode_ == 1)
      throw std::runtime_error("intentional std::exception");
    else
      throw 7;
  }

  void HandleException(const MTL::Exception&, const MTL::String&) override
  {
    HandledMTL_++;
  }
  void HandleException(const std::exception&, const MTL::String&) override
  {
    HandledStd_++;
  }
  void HandleException(const MTL::String&) override
  {
    HandledUnknown_++;
  }

  int Mode_;
  int HandledMTL_;
  int HandledStd_;
  int HandledUnknown_;
};

TEST(ChipAway_WorkerThread_HandleExceptionVariants)
{
  using Worker = WorkerThread<int, std::vector<int>, std::mutex>;
  for (int mode = 0; mode < 3; mode++)
  {
    ThrowingWorker<Worker> w(mode);
    w.Start();
    w.QueueWork(1);
    // Give the worker time to wake, dispatch ProcessWork, throw, and have the
    // ProcessThread catch handler call our HandleException override before we
    // race a Shutdown that would set Running_=false ahead of the wake-up.
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    w.Shutdown();

    if (mode == 0) MTL_EQUAL(w.HandledMTL_, 1);
    if (mode == 1) MTL_EQUAL(w.HandledStd_, 1);
    if (mode == 2) MTL_EQUAL(w.HandledUnknown_, 1);
  }
}

// AxisAngle with the identity rotation drives the zero-angle special case
// in ComputeCachedValues (lines 193-196 in AxisAngle.h).
TEST(ChipAway_AxisAngle_ZeroAngle)
{
  AxisAngle<F64> identity(Vector3D<F64>(0, 0, 0));
  Rotation3D<F64> R = identity.GetRotationMatrix();
  for (I32 i = 0; i < 3; i++)
    for (I32 j = 0; j < 3; j++)
      MTL_EQUAL_FLOAT(R[i][j], i == j ? 1.0 : 0.0, 1e-15);
}

// LUP with a row swap drives the sign-flip branch (Matrix.h line ~804) and
// Inverse on a singular matrix returns the NaN sentinel (line ~759).
TEST(ChipAway_Matrix_LUP_PivotSwapAndSingularInverse)
{
  // Row 0 has the small pivot; row 1 has a much larger leading magnitude, so
  // the first iteration of LUP swaps them and flips the sign.
  SquareMatrix<3, F64> A;
  A[0][0] = 0.001; A[0][1] = 1.0; A[0][2] = 2.0;
  A[1][0] = 5.0;   A[1][1] = 4.0; A[1][2] = 3.0;
  A[2][0] = 1.0;   A[2][1] = 2.0; A[2][2] = 9.0;

  SquareMatrix<3, F64> LU;
  ColumnVector<3, I32> P;
  F64 sign = 1.0;
  bool ok = A.LUP(LU, P, sign);
  MTL_VERIFY(ok);
  MTL_VERIFY(sign != 1.0);  // a row swap occurred.

  // Singular matrix: rank-deficient -> Inverse returns NaN sentinel.
  SquareMatrix<3, F64> S;
  S.Zeros();
  SquareMatrix<3, F64> Sinv = S.Inverse();
  MTL_VERIFY(Sinv[0][0] != Sinv[0][0]);  // NaN compares unequal to itself.
}

// Drives WorkerThread's Pause/Continue branches and IsRunning getter. We use
// PeriodMilliseconds=10 so Pause/Continue are no-ops at the running-flag level
// (they only affect Running_ when PeriodMilliseconds==0, which would also
// cause the worker's main loop to exit irreversibly).
template <class WorkerClass>
class CountingWorker : public WorkerClass
{
public:
  CountingWorker() : WorkerClass(std::wstring(L"chipaway")), Total_(0) {}
  void ProcessWork(const std::vector<int>& data) override
  {
    for (int v : data) Total_ += v;
  }
  int Total_;
};

TEST(ChipAway_WorkerThread_PauseContinueIsRunning)
{
  using Worker = WorkerThread<int, std::vector<int>, std::mutex, /*PeriodMs*/ 10>;
  CountingWorker<Worker> w;
  w.Start();
  MTL_VERIFY(w.IsRunning());

  // Calls into Pause/Continue to drive the PeriodMilliseconds!=0 short-circuit.
  w.Pause();
  w.Continue();
  MTL_VERIFY(w.IsRunning());

  for (int i = 1; i <= 4; i++)
    w.QueueWork(i);

  // Wait briefly for the periodic worker to drain its queue.
  for (int spin = 0; spin < 500 && !w.IsIdle(); spin++)
    std::this_thread::sleep_for(std::chrono::milliseconds(2));

  w.Shutdown();
  MTL_EQUAL(w.Total_, 1 + 2 + 3 + 4);
}
