# Math Template Library (MTL)

A header-only C++17 library for numerical and prototyping work. MTL gives you
matrix and vector types whose hot operations are already vectorized with SSE2,
AVX2/FMA, and optionally AVX-512, and parallelized with OpenMP — so a
research-style prototype runs at production-grade speed on day one.

## What MTL adds when you are prototyping

Most prototypes start with `std::vector` and a couple of nested loops, then
later get rewritten "for performance". MTL is designed so you do not have to
do that rewrite:

- **Drop-in replacements that are already fast.** `MTL::DynamicVector<T>` and
  `MTL::DynamicMatrix<T>` look like ordinary containers, but `+`, `-`, `*=`,
  dot products, sums, sums of squares, and matrix multiplies dispatch to
  SIMD- and OpenMP-accelerated kernels. You write straightforward linear
  algebra and get the parallel implementation for free.
- **A solver toolbox, not just a math kernel.** The `Math/` folder ships
  ready-to-use SVD, QR, LDLt, sparse matrices, Levenberg–Marquardt
  optimizers (dense and sparse), polynomial fitting, Givens/Householder
  rotations, sphere fitting, 2D/3D points and rotations, and an affine
  transform / projection-to-image stack. These are the pieces a typical
  computer-vision or estimation prototype needs, with consistent APIs.
- **Readable SIMD when you do drop down a level.** The `XX<T>` template
  (see below) wraps the intrinsics behind ordinary operators, so a
  hand-vectorized inner loop reads like normal C++ rather than `_mm256_*`
  soup. The same source compiles to SSE, AVX2, or AVX-512 by flipping a
  CMake flag.
- **A built-in test harness.** `MTL::Test`, `TEST(Name)`, `MTL_VERIFY`,
  `MTL_EQUAL_FLOAT`, and `MTL_APP()` let you turn a prototype into a
  self-checking executable in a few lines, with timing, colored output, a
  Unicode progress bar, and exception handling. Each test file in `Tests/`
  is auto-discovered by CMake and built into its own executable.
- **Concurrency utilities that match the style.** `Tools/` provides
  `WorkerThread`, `Pool`, `Event`, `SpinMutex`, and a small command-line
  options system (`Options.h`, `MTL_OPTION(...)`) so a prototype can grow
  into a small CLI tool without pulling in heavy frameworks.
- **Header-only and dependency-free.** Add `include/` to your include path
  and `#include <MTL/...>`. There is nothing to link, no transitive
  dependencies, and no build artifacts to manage.

The result is a short path from "I want to try this idea" to "I have a
fast, tested executable that I can compare against a baseline".

## Build Status

| Platform | Configuration | Status |
|----------|--------------|--------|
| Windows Server 2022 | Visual Studio 2022 | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Windows2022_VS2022.yml/badge.svg) |
| Ubuntu 22.04 | GCC | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22.yml/badge.svg) |
| Ubuntu 22.04 | GCC + AVX | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22_AVX.yml/badge.svg) |

## Getting started

```cpp
#include <MTL/Math/DynamicVector.h>
#include <MTL/Math/DynamicMatrix.h>
#include <MTL/Math/SVD.h>
#include <MTL/Tools/Test.h>

TEST(Solve_LeastSquares)
{
  MTL::DynamicMatrix<double> At = /* N-by-M, A transposed */;
  MTL::DynamicVector<double> b  = /* M */;

  // Vectorized + threaded under the hood; the call site stays simple.
  MTL::I32 rank;
  double conditionNumber;
  MTL::DynamicVector<double> x = b;
  MTL::SolveJacobiSVDTransposed(x, rank, conditionNumber, At, b);

  MTL_GREATER_THAN(rank, 0);
}
```

Build with CMake (3.14+):

```
cmake -B Build -DCMAKE_BUILD_TYPE=Release
cmake --build Build -j$(nproc)
python3 RunTests.py -ConsoleOut
```

CMake options: `MTL_ENABLE_SSE` (ON), `MTL_ENABLE_AVX` (ON),
`MTL_ENABLE_AVX512` (OFF). AVX-512 is off by default because not every CI
machine supports it.

## SIMD with the XX class

When a prototype needs a custom vectorized loop, the `XX<T>` template wraps
SSE / AVX / AVX-512 intrinsics behind a clean operator interface. It maps to
`X128<T>`, `X256<T>`, or `X512<T>` based on the enabled instruction set, so
the same source works at every SIMD width.

### Before (raw intrinsics)
```cpp
__m256d a = _mm256_loadu_pd(pSrc);
__m256d b = _mm256_set1_pd(2.0);
__m256d sum = _mm256_add_pd(a, b);
__m256d prod = _mm256_mul_pd(sum, sum);
__m256d result = _mm256_sqrt_pd(prod);
_mm256_storeu_pd(pDst, result);
```

### After (using XX)
```cpp
MTL::XX<MTL::F64> a(pSrc);
MTL::XX<MTL::F64> b(2.0);
MTL::XX<MTL::F64> sum = a + b;
MTL::XX<MTL::F64> result = (sum * sum).SquareRoot();
result.Store(pDst);
```

`XX<T>` supports the standard arithmetic (`+`, `-`, `*`, `/`), bitwise
(`&`, `|`, `^`), comparison (`==`, `<`, `>`, `<=`, `>=`), and compound
assignment (`+=`, `-=`, `*=`, `/=`) operators. Additional methods include
`SquareRoot()`, `Reciprocal()`, `Abs()`, `Min()`, `Max()`,
`MultiplyAndAdd()`, and `Conditional()`.

Supported element types: `F32`, `F64`, `I8`, `U8`, `I16`, `U16`, `I32`,
`U32`, `I64`, `U64`.

## Repository layout

- `include/MTL/Math/` — vectors, matrices, solvers, optimizers, transforms.
- `include/MTL/Stream/` — SIMD wrappers (`X128`/`X256`/`X512`/`XX`) and
  vectorized array kernels.
- `include/MTL/Tools/` — test framework, threading, progress bar, base64,
  worker pools, options parser.
- `Tests/` — one self-contained `Test*.cpp` per topic; CMake auto-discovers
  them.
- `ExampleApp/Example.cpp` — minimal `MTL_APP()` showing the options system.

