//
// Math Template Library
//
// Copyright (c) 2019: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/BinaryStream.h>
#include <MTL/File.h>
#include <MTL/Math.h>

using namespace MTL;

TEST(Test_Binary_Stream)
{
  const std::string tempFile = "_Temp_BinaryStream.bin";

  int8_t   a0 = -7;
  int16_t  a1 = -333;
  int32_t  a2 = -77777;
  int64_t  a3 = Pow<5>(-77777);
  uint8_t  a4 = 155;
  uint16_t a5 = 999;
  uint32_t a6 = 123456789;
  uint64_t a7 = Pow<6>(-77777);
  float    a8 =  0.111f;
  double   a9 = -1.234;

  BinaryStream stream;
  stream << a0 << a1 << a2 << a3 << a4 << a5 << a6 << a7 << a8 << a9;

  int8_t   b0;
  int16_t  b1;
  int32_t  b2;
  int64_t  b3;
  uint8_t  b4;
  uint16_t b5;
  uint32_t b6;
  uint64_t b7;
  float    b8;
  double   b9;
  stream >> b0 >> b1 >> b2 >> b3 >> b4 >> b5 >> b6 >> b7 >> b8 >> b9;

  MTL_EQUAL(a0, b0);
  MTL_EQUAL(a1, b1);
  MTL_EQUAL(a2, b2);
  MTL_EQUAL(a3, b3);
  MTL_EQUAL(a4, b4);
  MTL_EQUAL(a5, b5);
  MTL_EQUAL(a6, b6);
  MTL_EQUAL(a7, b7);
  MTL_EQUAL(a8, b8);
  MTL_EQUAL(a9, b9);

  stream.Save(tempFile);
  int8_t   c0;
  int16_t  c1;
  int32_t  c2;
  int64_t  c3;
  uint8_t  c4;
  uint16_t c5;
  uint32_t c6;
  uint64_t c7;
  float    c8;
  double   c9;
  BinaryStream stream2;
  stream2.Load(tempFile);
  stream2 >> c0 >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9;

  MTL_EQUAL(a0, c0);
  MTL_EQUAL(a1, c1);
  MTL_EQUAL(a2, c2);
  MTL_EQUAL(a3, c3);
  MTL_EQUAL(a4, c4);
  MTL_EQUAL(a5, c5);
  MTL_EQUAL(a6, c6);
  MTL_EQUAL(a7, c7);
  MTL_EQUAL(a8, c8);
  MTL_EQUAL(a9, c9);

  File::RemoveIfExists(tempFile);
}

TEST(Test_vector)
{
  std::vector<int> u(7, 1234);
  std::shared_ptr<std::vector<float>> nullPtr;
  std::shared_ptr<std::vector<int>> uu(new std::vector<int>(u));

  BinaryStream stream;
  stream << u << nullPtr << uu;

  std::vector<int> v;
  std::shared_ptr<std::vector<float>> empty;
  std::shared_ptr<std::vector<int>> vv;
  stream >> v >> empty >> vv;

  MTL_VERIFY(empty.get() == nullptr);

  for (uint32_t i = 0; i < u.size(); i++)
  {
    MTL_EQUAL(u[i], v[i]);
    MTL_EQUAL(u[i], (*vv)[i]);
  }
}

TEST(Test_DynamicVector)
{
  DynamicVector<double> empty;
  DynamicVector<int> u(7, 1234);
  BinaryStream stream;
  stream << empty << u;

  DynamicVector<double> bogus;
  DynamicVector<int> v;
  stream >> bogus >> v;

  MTL_EQUAL((int)bogus.Size(), 0);

  for (uint32_t i = 0; i < u.Size(); i++)
    MTL_EQUAL(u[i], v[i]);

  stream.ReadPosition(0);
  std::vector<int> vv;
  stream >> bogus >> vv;

  for (uint32_t i = 0; i < u.Size(); i++)
    MTL_EQUAL(u[i], vv[i]);
}

TEST(Test_Strings)
{
  std::vector<std::string> strings8;
  std::vector<std::wstring> strings16;
  strings16.push_back(L"Boggis");
  strings16.push_back(L"Bunce");
  strings16.push_back(L"Bean");
  strings16.push_back(L"");
  strings16.push_back(L"!@#$%^&*()_-+=");

  for (const auto& s : strings16)
    strings8.push_back(ToUTF8(s));

  BinaryStream stream;
  stream << strings8 << strings16;

  std::vector<std::string> strings_8;
  std::vector<std::wstring> strings_16;
  stream >> strings_8 >> strings_16;

  for (uint32_t i = 0; i < strings16.size(); i++)
  {
    MTL_VERIFY(strings_8[i] == strings8[i]);
    MTL_VERIFY(strings_16[i] == strings16[i]);
  }
}
