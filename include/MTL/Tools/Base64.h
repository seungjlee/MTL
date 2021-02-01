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


#ifndef MTL_BASE_64_H
#define MTL_BASE_64_H

#include <vector>
#include <string>
#include <cstdint>

namespace MTL
{
static const char Alphabet64[] =
  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  "abcdefghijklmnopqrstuvwxyz"
  "0123456789+/";

static uint8_t InverseBase64[128];
static int InitializeInverseBase64()
{
  for (uint8_t i = 0; i < sizeof(Alphabet64); i++)
    InverseBase64[Alphabet64[i]] = i;

  return 0;
}
static int Initialize_Inverse_Base_64 = InitializeInverseBase64();

class Base64
{
public:
  static std::string Encode(const std::vector<char>& v)
  {
    return Encode((uint8_t*)v.data(), v.size());
  }
  static std::string Encode(const std::vector<uint8_t>& v)
  {
    return Encode(v.data(), v.size());
  }
  static std::string Encode(const char* pSrc, size_t size)
  {
    return Encode((uint8_t*)pSrc, size);
  }
  static std::string Encode(const uint8_t* pSrc, size_t size)
  {
    std::string base64;
    Encode(base64, pSrc, size);
    return base64;
  }
  static void Encode(std::string& base64, const uint8_t* pSrc, size_t size)
  {
    uint8_t buffer4[4];
    size_t size64 = ((size + 2) / 3) * 4;
    base64.resize(size64);

    int64_t srcBuffer;
    uint8_t* pSrcBuffer = (uint8_t*)&srcBuffer;

    const uint8_t* pSrcEnd = pSrc + size;
    char* pDst = &base64[0];
    for (uint32_t i = 0; i < size64; i += 4, pSrc += 3)
    {
      srcBuffer = 0;
      pSrcBuffer[0] = pSrc[0];
      srcBuffer <<= 6;
      buffer4[0] = pSrcBuffer[1];

      if (pSrc + 1 < pSrcEnd)
      {
        srcBuffer <<= 2;
        pSrcBuffer[0] = pSrc[1];
        srcBuffer <<= 4;
        buffer4[1] = pSrcBuffer[1] & 0x3f;

        if (pSrc + 2 < pSrcEnd)
        {
          srcBuffer <<= 4;
          pSrcBuffer[0] = pSrc[2];
          buffer4[3] = pSrcBuffer[0] & 0x3f;
          srcBuffer <<= 2;
          buffer4[2] = pSrcBuffer[1] & 0x3f;

          pDst[i + 0] = Alphabet64[buffer4[0]];
          pDst[i + 1] = Alphabet64[buffer4[1]];
          pDst[i + 2] = Alphabet64[buffer4[2]];
          pDst[i + 3] = Alphabet64[buffer4[3]];
        }
        else
        {
          srcBuffer <<= 6;
          buffer4[2] = pSrcBuffer[1] & 0x3f;

          pDst[i + 0] = Alphabet64[buffer4[0]];
          pDst[i + 1] = Alphabet64[buffer4[1]];
          pDst[i + 2] = Alphabet64[buffer4[2]];
          pDst[i + 3] = '=';
        }
      }
      else
      {
        srcBuffer <<= 6;
        buffer4[1] = pSrcBuffer[1] & 0x3f;

        pDst[i + 0] = Alphabet64[buffer4[0]];
        pDst[i + 1] = Alphabet64[buffer4[1]];
        pDst[i + 2] = '=';
        pDst[i + 3] = '=';
      }
    }
  }

  static std::vector<uint8_t> Decode(const std::string& src)
  {
    return Decode(src.data(), src.size());
  }
  static std::vector<uint8_t> Decode(const char* pSrc, size_t size)
  {
    std::vector<uint8_t> decoded;
    Decode(decoded, pSrc, size);
    return decoded;
  }
  static void Decode(std::vector<uint8_t>& decoded, const char* pSrc, size_t size)
  {
    int padding = 0;
    if (pSrc[size - 1] == '=')
    {
      padding = 1;
      if (pSrc[size - 2] == '=')
        padding = 2;
    }

    int decodeSize = int(3 * size / 4) - padding;

    int64_t srcBuffer;
    uint8_t* pSrcBuffer = (uint8_t*)&srcBuffer;

    decoded.resize(decodeSize);
    uint8_t* pDst = decoded.data();

    uint32_t i = 0;
    for (; i < size - 4; i += 4, pDst += 3)
    {
      pSrcBuffer[0] = InverseBase64[pSrc[i + 0]];
      srcBuffer <<= 6;
      pSrcBuffer[0] |= InverseBase64[pSrc[i + 1]];
      srcBuffer <<= 6;
      pSrcBuffer[0] |= InverseBase64[pSrc[i + 2]];
      srcBuffer <<= 6;
      pSrcBuffer[0] |= InverseBase64[pSrc[i + 3]];
      pDst[0] = pSrcBuffer[2];
      pDst[1] = pSrcBuffer[1];
      pDst[2] = pSrcBuffer[0];
    }

    pSrcBuffer[0] = InverseBase64[pSrc[i + 0]];
    srcBuffer <<= 6;
    pSrcBuffer[0] |= InverseBase64[pSrc[i + 1]];
    srcBuffer <<= 6;
    pSrcBuffer[0] |= InverseBase64[pSrc[i + 2]];
    srcBuffer <<= 6;
    pSrcBuffer[0] |= InverseBase64[pSrc[i + 3]];
    pDst[0] = pSrcBuffer[2];
    if (padding < 2)
      pDst[1] = pSrcBuffer[1];
    if (padding == 0)
      pDst[2] = pSrcBuffer[0];
  }
};

static std::string ToString(const std::vector<uint8_t>& v)
{
  return std::string((char*)v.data(), v.size());
}

}  // namespace MTL

#endif  // MTL_BASE_64_H
