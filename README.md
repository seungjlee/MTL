# Math Template Library  
Matrix and vector templated classes with operations accelerated using SSE2, AVX2, FMA, and OpenMP. Contains efficient methods to solve linear and non-linear least square problems with great numerical stability. Note that AVX and FMA acceleration is disabled by default since not many systems support it currently. These features can be easily enabled for the included tests using the CMake GUI.  

## Build Status

| Platform | Configuration | Status |
|----------|--------------|--------|
| Windows Server 2022 | Visual Studio 2022 | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Windows2022_VS2022.yml/badge.svg) |
| Ubuntu 22.04 | GCC | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22.yml/badge.svg) |
| Ubuntu 22.04 | GCC + AVX | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22_AVX.yml/badge.svg) |

## SIMD with the XX Class

One of the main goals of this library is to make SIMD intrinsics easier to use and read. The `XX<T>` template class wraps SSE, AVX, and AVX512 intrinsics behind a clean C++ operator interface. It automatically maps to `X128<T>`, `X256<T>`, or `X512<T>` depending on which instruction set is enabled at compile time, so the same code works across all SIMD widths without changes.

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

The `XX<T>` class supports all standard arithmetic (`+`, `-`, `*`, `/`), bitwise (`&`, `|`, `^`), comparison (`==`, `<`, `>`, `<=`, `>=`), and compound assignment (`+=`, `-=`, `*=`, `/=`) operators. Additional methods include `SquareRoot()`, `Reciprocal()`, `Abs()`, `Min()`, `Max()`, `MultiplyAndAdd()`, and `Conditional()`.

Supported element types: `F32`, `F64`, `I8`, `U8`, `I16`, `U16`, `I32`, `U32`, `I64`, `U64`.
