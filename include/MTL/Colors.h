//
// Math Template Library
//
// Copyright (c) 2019: Seung Jae Lee, https://sourceforge.net/projects/mathtemplatelibrary/
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


#ifndef MTL_COLORS_H
#define MTL_COLORS_H

#include <MTL/Definitions.h>
#include <stdio.h>

// Color macros.
#define COLOR_RESET     L"\033[0m"

#define COLOR_BLACK     L"\033[30m"
#define COLOR_RED       L"\033[31m"
#define COLOR_GREEN     L"\033[32m"
#define COLOR_YELLOW    L"\033[33m"
#define COLOR_BLUE      L"\033[34m"
#define COLOR_MAGENTA   L"\033[35m"
#define COLOR_CYAN      L"\033[36m"
#define COLOR_WHITE     L"\033[37m"

#define COLOR_LBLACK    L"\033[90m"
#define COLOR_LRED      L"\033[91m"
#define COLOR_LGREEN    L"\033[92m"
#define COLOR_LYELLOW   L"\033[93m"
#define COLOR_LBLUE     L"\033[94m"
#define COLOR_LMAGENTA  L"\033[95m"
#define COLOR_LCYAN     L"\033[96m"
#define COLOR_LWHITE    L"\033[97m"

// RGB background and foreground colors.
#define COLOR_BG(R,G,B)  L"\033[48;2;" L ## #R L";" L ## #G L";" L ## #B L"m"
#define COLOR_FG(R,G,B)  L"\033[38;2;" L ## #R L";" L ## #G L";" L ## #B L"m"

namespace MTL
{

// Color helper class.
class ColorScope
{
public:
  ColorScope(const wchar_t* color)
  {
    wprintf(color);
  }
  ~ColorScope()
  {
    wprintf(COLOR_RESET);
  }
};

}

#endif  // MTL_COLORS_H
