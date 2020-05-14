# Math Template Library  
Matrix and vector templated classes with operations accelerated using SSE2, AVX2, FMA, and OpenMP. Contains efficient methods to solve linear and non-linear least square problems with great numerical stability. Note that AVX and FMA acceleration is disabled by default since not many systems support it currently. These features can be easily enabled for the included tests using the CMake GUI.  
Currently, this software is getting tested with Microsoft Visual Studio 2015 and 2017. Also g++ 7.5 is being tested on Ubuntu.  
