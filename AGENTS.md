# Project Guidelines

## Architecture

Header-only C++ math template library in `include/MTL/`. Key areas:

- **Math/** — Matrix, vector, linear algebra, optimization (LM, sparse), geometric transforms.
- **Stream/** — SIMD wrappers (`X128`/`X256`/`X512` for SSE/AVX/AVX512), stream array operations.
- **Tools/** — Threading, synchronization, test framework, utilities.

`XX<T>` is the width-agnostic SIMD alias (maps to X128, X256, or X512 based on compile flags).

## Build and Test

CMake 3.14+, C++17. Build with:

```
cmake -B Build -G "Visual Studio 17 2022" -A x64   # Windows
cmake -B Build -DCMAKE_BUILD_TYPE=Release           # Linux
cmake --build Build --config Release
```

CMake options: `MTL_ENABLE_SSE` (ON), `MTL_ENABLE_AVX` (ON), `MTL_ENABLE_AVX512` (OFF).

Run all tests:

```
python RunTests.py -ConsoleOut          # Windows
python3 RunTests.py -ConsoleOut         # Linux
```

Each `Tests/Test*.cpp` compiles to a separate executable. Tests are auto-discovered by the CMake glob—just add a new file.

## Code Style

- **Namespace:** All library code in `MTL`.
- **Types:** Use `MTL::I8`/`I16`/`I32`/`I64`, `MTL::U8`/`U16`/`U32`/`U64`, `MTL::F32`/`F64` instead of primitive types.
- **Naming:** PascalCase for classes, methods, functions, and local variables. Constants use `k` prefix (`kTol`, `kInfinity32`). Member variables have trailing underscore (`Data_`, `Name_`). Macros use `MTL_` prefix with `SCREAMING_CASE`.
- **Indentation:** 2 spaces.
- **Comments:** Minimal. Use `//` only where logic is not self-evident.

## Conventions

- SIMD code uses `#if MTL_ENABLE_SSE` / `#if MTL_ENABLE_AVX` / `#if MTL_ENABLE_AVX512` guards.
- Prefer `XX<T>` for tests that work at any SIMD width. Use explicit `X128`/`X256`/`X512` only when testing width-specific features.
- Test macros: `TEST(Name)`, `MTL_VERIFY(expr)`, `MTL_EQUAL(actual, expected)`, `MTL_EQUAL_FLOAT(actual, expected, tol)`.
- Every source file starts with the BSD 2-Clause license header.
- Wide strings (`std::wstring`, `wchar_t`) are used throughout the library.
