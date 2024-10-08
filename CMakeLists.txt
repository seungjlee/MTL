#
# Math Template Library
#
# Copyright (c) 2014: Seung Jae Lee, https://github.com/seungjlee/MTL
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice, this list of
#      conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright notice, this list of
#      conditions and the following disclaimer in the documentation and/or other materials
#      provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PROJECT(MTL)

cmake_minimum_required(VERSION 3.14)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

include(MTL)

option (MTL_ENABLE_SSE    "Enables SSE instructions."          ON)
option (MTL_ENABLE_AVX    "Enables AVX and FMA instructions."  ON)
option (MTL_ENABLE_AVX512 "Enables AVX512 instructions."      OFF)
option (MTL_USE_FAST_FMA  "Use faster version of std::fma()."  ON)

message(STATUS "")

if(MSVC)

  SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /Zm1000 /EHsc /GR /openmp /FC /W4 /wd4127 /wd4505 /std:c++17")

  if(MTL_ENABLE_AVX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  endif(MTL_ENABLE_AVX)
  if(MTL_ENABLE_AVX512)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
  endif(MTL_ENABLE_AVX512)

  SET(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MDd /Zi /Ob0 /Od /RTC1")
  SET(CMAKE_CXX_FLAGS_RELEASE "/MD /D NDEBUG /Ox")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MD /Zi /D NDEBUG /Ox /Oy-")

  if(CMAKE_CL_64)
    ADD_DEFINITIONS(-DWIN64)
  endif(CMAKE_CL_64)

  SET(WINVER 0x0502)

  #ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS)
  ADD_DEFINITIONS(-DNOMINMAX -DUNICODE "-DWINVER=${WINVER}" "-D_WIN32_WINNT=${WINVER}")
  ADD_DEFINITIONS(-D_ENABLE_ATOMIC_ALIGNMENT_FIX -D_SILENCE_CXX17_CODECVT_HEADER_DEPRECATION_WARNING)

else(MSVC)

  if(NOT DEFINED CMAKE_BUILD_TYPE OR NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting CMAKE_BUILD_TYPE=Release")
    SET(CMAKE_BUILD_TYPE Release)
  endif()

  message(STATUS "")
  message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
  message(STATUS "")

  SET(CMAKE_CXX_FLAGS "-fopenmp -std=c++17 -Wall")

  if (MTL_ENABLE_AVX512)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mavx512dq -mavx512f -mfma")
  elseif (MTL_ENABLE_AVX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
  endif()

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-char-subscripts")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-strict-aliasing")

  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")
endif(MSVC)

if(MTL_ENABLE_SSE)
  ADD_DEFINITIONS(-DMTL_ENABLE_SSE=1)
else(MTL_ENABLE_SSE)
  ADD_DEFINITIONS(-DMTL_ENABLE_SSE=0)
endif(MTL_ENABLE_SSE)

if(MTL_ENABLE_AVX)
  ADD_DEFINITIONS(-DMTL_ENABLE_AVX=1)
else(MTL_ENABLE_AVX)
  ADD_DEFINITIONS(-DMTL_ENABLE_AVX=0)
endif(MTL_ENABLE_AVX)

if(MTL_ENABLE_AVX512)
  ADD_DEFINITIONS(-DMTL_ENABLE_AVX512=1)
else(MTL_ENABLE_AVX512)
  ADD_DEFINITIONS(-DMTL_ENABLE_AVX512=0)
endif(MTL_ENABLE_AVX512)

if(MTL_USE_FAST_FMA)
  ADD_DEFINITIONS(-DFP_FAST_FMA -DFP_FAST_FMAF -DFP_FAST_FMAL)
endif(MTL_USE_FAST_FMA)

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "")

message(STATUS "Compiler definitions:")
get_directory_property(DirDefs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMPILE_DEFINITIONS)
foreach(d ${DirDefs})
  message(STATUS "  ${d}" )
endforeach()
message(STATUS "")

include_directories(
                    ./include
                    ./Tests
                   )

add_subdirectory(Tests)
