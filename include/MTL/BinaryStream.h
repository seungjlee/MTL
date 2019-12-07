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


#ifndef MTL_BINARY_STREAM_H
#define MTL_BINARY_STREAM_H

#include <MTL/DynamicVector.h>
#include <MTL/StringHelpers.h>
#include <vector>
#include <cstdint>
#include <assert.h>

namespace MTL
{
class BinaryStream
{
public:
  BinaryStream()
    : ReadPosition_(0)
  {
  }

  void Read(uint8_t* p, size_t size) const
  {
    assert(Data_.begin() + size <= Data_.end());
    memcpy(p, &Data_[ReadPosition_], size);
    const_cast<BinaryStream*>(this)->ReadPosition_ += size;
  }
  void Write(const uint8_t* p, size_t size)
  {
    Data_.insert(Data_.end(), p, p + size);
  }

  template <class T> BinaryStream& operator<<(const T& val)
  {
    Write((uint8_t*)&val, sizeof(T));
    return *this;
  }
  template <class T> const BinaryStream& operator>>(T& val) const 
  {
    Read((uint8_t*)&val, sizeof(T));
    return *this;
  }

  template <class T> BinaryStream& operator<<(const std::basic_string<T>& str)
  {
    *this << uint64_t(str.size());
    Write((uint8_t*)str.data(), str.size() * sizeof(str[0]));

    return *this;
  }
  template <class T> const BinaryStream& operator>>(std::basic_string<T>& str) const
  {
    uint64_t size;
    *this >> size;
    str.resize(size);
    Read((uint8_t*)str.data(), str.size() * sizeof(str[0]));

    return *this;
  }
  template <class T> BinaryStream& operator<<(const std::vector<T>& v)
  {
    *this << uint64_t(v.size());
    for (const auto& x : v)
      *this << x;

    return *this;
  }
  template <class T> const BinaryStream& operator>>(std::vector<T>& v) const
  {
    uint64_t size;
    *this >> size;
    v.resize(size);
    for (auto& x : v)
      *this >> x;

    return *this;
  }
  template <class T> BinaryStream& operator<<(const DynamicVector<T>& v)
  {
    *this << uint64_t(v.Size());
    for (const auto& x : v)
      *this << x;

    return *this;
  }
  template <class T> const BinaryStream& operator>>(DynamicVector<T>& v) const
  {
    uint64_t size;
    *this >> size;
    v.Resize(size);
    for (auto& x : v)
      *this >> x;

    return *this;
  }

  const std::vector<uint8_t>& Data() const  { return Data_; }
  size_t ReadPosition() const         { return ReadPosition_;     }
  void ReadPosition(size_t position)  { ReadPosition_ = position; }

  void Save(const std::wstring& fileName) const
  {
    Save(ToUTF8(fileName));
  }
  void Load(const std::wstring& fileName)
  {
    Load(ToUTF8(fileName));
  }
  void Save(const std::string& fileName) const
  {
#ifdef WIN32
    std::FILE* f;
    fopen_s(&f, fileName.c_str(), "wb");
#else
    std::FILE* f = fopen(fileName.c_str(), "wb");
#endif
    uint64_t size = Data_.size();
    fwrite(&size, sizeof(size), 1, f);
    fwrite(Data_.data(), 1, Data_.size(), f);
    fclose(f);
  }
  void Load(const std::string& fileName)
  {
#ifdef WIN32
    std::FILE* f;
    fopen_s(&f, fileName.c_str(), "rb");
#else
    std::FILE* f = fopen(fileName.c_str(), "rb");
#endif
    uint64_t size;
    size_t count = fread(&size, sizeof(size), 1, f);
    if (count != 1)
      MTL_THROW("Invalid binary file!");

    Data_.resize(size);
    count = fread(Data_.data(), 1, Data_.size(), f);
    if (count != Data_.size())
      MTL_THROW("Invalid binary file!");

    ReadPosition_ = 0;
    fclose(f);
  }

private:
  std::vector<uint8_t> Data_;
  size_t ReadPosition_;
};

}  // namespace MTL

#endif  // MTL_BINARY_STREAM_H
