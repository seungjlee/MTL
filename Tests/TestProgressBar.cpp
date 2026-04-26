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

#include <MTL/Tools/Test.h>
#include <MTL/Tools/ProgressBar.h>

using namespace MTL;

// Drive a ProgressBar through a full cycle of updates so the rendering paths
// (full cells, partial sub-cells, empty trailing cells, both precisions, both
// the synchronous and asynchronous final-update flavors) all execute. The
// rendering writes to wcout; tests run with -DisableProgressBar in CI but here
// we explicitly construct the bar to exercise its code paths.
TEST(Test_ProgressBar_FullCycle)
{
  ProgressBar bar(/*finalUpdateIsSynchronous=*/true);

  // Drive a few partial-block fractions explicitly.
  bar.Update(0.0,    L"start", false);
  bar.Update(0.125,  L"1/8",   false);
  bar.Update(0.3125, L"5/16",  true);   // extra precision path
  bar.Update(0.5,    L"half",  false);
  bar.Update(0.875,  L"7/8",   false);
  bar.Update(1.0,    L"done",  false);  // synchronous final update
}

TEST(Test_ProgressBar_AsynchronousFinal)
{
  ProgressBar bar(/*finalUpdateIsSynchronous=*/false);
  bar.Update(0.0,  L"begin");
  bar.Update(0.5,  L"middle");
  bar.Update(1.0,  L"end");
}

TEST(Test_ProgressBar_DisableEnable)
{
  ProgressBar bar(false);
  bar.Disable();
  bar.Update(0.5, L"silent");  // should be ignored.
  bar.Enable();
  bar.Update(0.6, L"loud");
  bar.Update(1.0, L"end");
}

TEST(Test_ProgressBar_OutOfRangeAndZeroLength)
{
  ProgressBar bar(false);
  // Out-of-range percentages are clamped.
  bar.Update(-0.5, L"clamped low");
  bar.Update( 1.5, L"clamped high");
  // Zero bar length skips the bar-drawing block entirely.
  bar.Update(0.5, L"no bar", false, 0);
}
