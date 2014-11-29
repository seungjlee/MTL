#
# Math Template Library
#
# Copyright (c) 2014: Seung Jae Lee
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

cmake_minimum_required(VERSION 3.0)

if(MSVC11)

  SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /EHsc /GR /openmp /FC")
  SET(CMAKE_CXX_FLAGS_DEBUG "/D_DEBUG /MDd /Zi /Ob0 /Od /RTC1")

  SET(CMAKE_CXX_FLAGS_RELEASE "/MD /D NDEBUG /Ox")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MD /Zi /D NDEBUG /Ox /Oy-")

  if(CMAKE_CL_64)
    ADD_DEFINITIONS(-DWIN64)
  endif(CMAKE_CL_64)
  
endif(MSVC11)

SET(WINVER 0x0501)

#ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -DNOMINMAX -D_SCL_SECURE_NO_DEPRECATE)
ADD_DEFINITIONS(-DUNICODE -DWNT "-DWINVER=${WINVER}" "-D_WIN32_WINNT=${WINVER}")

include_directories(
                    ./include
                    ./Tests
                   )

add_subdirectory(Tests)
