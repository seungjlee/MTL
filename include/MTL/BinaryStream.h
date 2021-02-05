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

#include <MTL/Exception.h>
#include <assert.h>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>

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
    assert(ReadPosition_ + size <= Data_.size());
    memcpy(p, Data_.data() + ReadPosition_, size);
    const_cast<BinaryStream*>(this)->ReadPosition_ += size;
  }
  void Write(const uint8_t* p, size_t size)
  {
    Data_.insert(Data_.end(), p, p + size);
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
    std::ofstream fileStream(fileName.c_str(),
                             std::ofstream::trunc | std::ofstream::binary);
    uint64_t size = Data_.size();
    fileStream.write((char*)&size, sizeof(size));
    fileStream.write((char*)Data_.data(), Data_.size());
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

template <class T>
static BinaryStream& operator<<(BinaryStream& stream, const T& val)
{
  static_assert(std::is_trivial<T>::value && std::is_standard_layout<T>::value, "Unspecialized type must be POD.");
  stream.Write((uint8_t*)&val, sizeof(T));
  return stream;
}
template <class T>
static const BinaryStream& operator>>(const BinaryStream& stream, T& val) 
{
  static_assert(std::is_trivial<T>::value && std::is_standard_layout<T>::value, "Unspecialized type must be POD.");
  stream.Read((uint8_t*)&val, sizeof(T));
  return stream;
}

template <class T>
static BinaryStream& operator<<(BinaryStream& stream, const std::basic_string<T>& str)
{
  stream << uint64_t(str.size());
  stream.Write((uint8_t*)str.data(), str.size() * sizeof(str[0]));

  return stream;
}
template <class T>
static const BinaryStream& operator>>(const BinaryStream& stream, std::basic_string<T>& str)
{
  uint64_t size;
  stream >> size;
  str.resize(size);
  stream.Read((uint8_t*)str.data(), str.size() * sizeof(str[0]));

  return stream;
}
template <class T>
static BinaryStream& operator<<(BinaryStream& stream, const std::vector<T>& v)
{
  stream << uint64_t(v.size());
  for (const auto& x : v)
    stream << x;

  return stream;
}
template <class T>
static const BinaryStream& operator>>(const BinaryStream& stream, std::vector<T>& v)
{
  uint64_t size;
  stream >> size;
  v.resize(size);
  for (auto& x : v)
    stream >> x;

  return stream;
}
template <class T>
static BinaryStream& operator<<(BinaryStream& stream, const std::set<T>& s)
{
  stream << uint64_t(s.size());
  for (const auto& x : s)
    stream << x;

  return stream;
}
template <class T>
static const BinaryStream& operator>>(const BinaryStream& stream, std::set<T>& s)
{
  uint64_t size;
  stream >> size;
  for (int i = 0; i < size; i++)
  {
    T x;
    stream >> x;
    s.insert(x);
  }

  return stream;
}
template <class T>
static BinaryStream& operator<<(BinaryStream& stream, const std::shared_ptr<T>& p)
{
  if (p == nullptr)
    stream << uint8_t(0);
  else
    stream << uint8_t(1) << *p;

  return stream;
}
template <class T>
static const BinaryStream& operator>>(const BinaryStream& stream, std::shared_ptr<T>& p)
{
  uint8_t byte;
  stream >> byte;
  if (byte != 0)
  {
    p.reset(new T());
    stream >> *p;
  }
  else
  {
    p.reset();
  }
  
  return stream;
}

}  // namespace MTL

#endif  // MTL_BINARY_STREAM_H
