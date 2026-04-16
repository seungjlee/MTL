# Math Template Library  
Matrix and vector templated classes with operations accelerated using SSE2, AVX2, FMA, and OpenMP. Contains efficient methods to solve linear and non-linear least square problems with great numerical stability. Note that AVX and FMA acceleration is disabled by default since not many systems support it currently. These features can be easily enabled for the included tests using the CMake GUI.  

## Build Status

| Platform | Configuration | Status |
|----------|--------------|--------|
| Windows Server 2022 | Visual Studio 2022 | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Windows2022_VS2022.yml/badge.svg) |
| Ubuntu 22.04 | GCC | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22.yml/badge.svg) |
| Ubuntu 22.04 | GCC + AVX | ![Build](https://github.com/seungjlee/MTL/actions/workflows/Ubuntu22_AVX.yml/badge.svg) |
