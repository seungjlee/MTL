//
// Math Template Library
//
// Copyright (c) 2016: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/StringHelpers.h>

using namespace MTL;

TEST(TestStringFunctions)
{
  ResetOutputStream();
  std::cout << ToLowerCase("AbC\n");
  std::cout << ToUpperCase("AbC\n");
  ResetOutputStream();

  std::wcout << ToLowerCase(L"AbC\n");
  std::wcout << ToUpperCase(L"AbC\n");

  std::wcout << ToUTF16(ToUTF8(ToLowerCase(L"AbC\n")));
  std::wcout << ToUTF16(ToUTF8(ToUpperCase(L"AbC\n")));

  MTL_EQUAL(ToUTF16(ToLowerCase("xYz")), ToUTF16(ToLowerCase("xyZ")));
  MTL_EQUAL(ToLowerCase(L"xYz"), ToLowerCase(L"xyZ"));
}
