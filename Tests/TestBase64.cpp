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
#include <MTL/Tools/Base64.h>

using namespace MTL;

TEST(Test_Encode_Decode)
{
  // From Wikipedia example.
  const char* testString = "any carnal pleasure.";
  Out() << ToUTF16(Base64::Encode(testString, 16)) << std::endl;
  Out() << ToUTF16(Base64::Encode(testString, 17)) << std::endl;
  Out() << ToUTF16(Base64::Encode(testString, 18)) << std::endl;
  Out() << ToUTF16(Base64::Encode(testString, 19)) << std::endl;
  Out() << ToUTF16(Base64::Encode(testString, 20)) << std::endl;

  MTL_VERIFY(Base64::Encode(testString, 16) == "YW55IGNhcm5hbCBwbGVhcw==");
  MTL_VERIFY(Base64::Encode(testString, 17) == "YW55IGNhcm5hbCBwbGVhc3U=");
  MTL_VERIFY(Base64::Encode(testString, 18) == "YW55IGNhcm5hbCBwbGVhc3Vy");
  MTL_VERIFY(Base64::Encode(testString, 19) == "YW55IGNhcm5hbCBwbGVhc3VyZQ==");
  MTL_VERIFY(Base64::Encode(testString, 20) == "YW55IGNhcm5hbCBwbGVhc3VyZS4=");

  for (int length = 11; length < 21; length++)
  {
    Out() << ToUTF16(ToString(Base64::Decode(Base64::Encode(testString, length)))) << std::endl;
    MTL_VERIFY(ToString(Base64::Decode(Base64::Encode(testString, length))) == std::string(testString, length));
  }
}
