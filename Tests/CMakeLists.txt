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

PROJECT(TestMTL)

FILE(GLOB headers ../include/MTL/*.h)
FILE(GLOB math_headers ../include/MTL/Math/*.h)
FILE(GLOB stream_headers ../include/MTL/Stream/*.h)
FILE(GLOB tools_headers ../include/MTL/Tools/*.h)
FILE(GLOB ALL_TESTS "Test*.cpp")

SOURCE_GROUP(MTL FILES ${headers})
SOURCE_GROUP(MTL\\Math FILES ${math_headers})
SOURCE_GROUP(MTL\\Stream FILES ${stream_headers})
SOURCE_GROUP(MTL\\Tools FILES ${tools_headers})

FOREACH(t ${ALL_TESTS})
  STRING(REPLACE ".cpp" "" TestNameWithPath ${t})
  STRING(REGEX MATCH "([0-9_a-zA-Z]+)$" TestName ${TestNameWithPath})
  ADD_EXECUTABLE( ${TestName} ${t} ${headers} ${math_headers} ${stream_headers} ${tools_headers})
ENDFOREACH(t ${ALL_TESTS})
