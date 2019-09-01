//
// Math Template Library
//
// Copyright (c) 2014: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_CPU_H
#define MTL_CPU_H

#include "Definitions.h"
#include <bitset>
#include <array>
#include <vector>
#include <omp.h>

#ifndef __ARM_ARCH

#ifdef WIN32
#include <intrin.h>
#else
#include <cpuid.h>
#include <string.h>
#endif

#endif

namespace MTL
{

class CPU;
static CPU* Instance_ = 0;

class CPU
{
public:
  static CPU& Instance()
  {
    if (!Instance_)
    {
      static CPU instance;
      Instance_ = &instance;
    }

    return *Instance_;
  }

  U64 NumberOfThreads() const  { return NumberOfThreads_; }
  void NumberOfThreads(U64 n)  { NumberOfThreads_ = n;    }
  U64 NumberOfCores() const    { return NumberOfCores_;   }

  bool SSE()  const  { return CPU_Rep_.f_1_EDX_[25] || false; }
  bool SSE2() const  { return CPU_Rep_.f_1_EDX_[26] || false; }
  bool AVX()  const  { return CPU_Rep_.f_1_ECX_[28] || false; }
  bool AVX2() const  { return CPU_Rep_.f_7_EBX_[5]  || false; }
  bool FMA()  const  { return CPU_Rep_.f_1_ECX_[12] || false; }

  bool IsIntel() const         { return CPU_Rep_.isIntel_;        }

  // Intel's hyperthreading for example.
  bool Multithreading() const  { return CPU_Rep_.multithreading_; }

private:
  CPU() : NumberOfThreads_(1)
  {
#ifndef __ARM_ARCH
    NumberOfCores_ = omp_get_num_procs();
    omp_set_num_threads((int)NumberOfCores_);

    // We usually want the actual number of cores.
    if (IsIntel() && Multithreading())
      NumberOfCores_ /= 2;

    NumberOfThreads(NumberOfCores_);
#endif
  }

  ~CPU() {}
  CPU(const CPU&);
  void operator=(const CPU&);

  // Adapted from Microsoft's example code.
  class InstructionSet_Internal
  {
  public:
    InstructionSet_Internal()
      : nIds_(0),
        nExIds_(0),
        isIntel_(false),
        isAMD_(false),
        f_1_ECX_(0),
        f_1_EDX_(0),
        f_7_EBX_(0),
        f_7_ECX_(0),
        f_81_ECX_(0),
        f_81_EDX_(0),
        data_(),
        extdata_()
    {
#ifndef __ARM_ARCH
      std::array<int, 4> cpui;

#ifdef WIN32
      // Calling __cpuid with 0x0 as the function_id argument
      // gets the number of the highest valid function ID.
      __cpuid(cpui.data(), 0);
      nIds_ = cpui[0];
#else
      nIds_ = __get_cpuid_max(0, 0);
#endif

      for (int i = 0; i <= nIds_; ++i)
      {
#ifdef WIN32
        __cpuidex(cpui.data(), i, 0);
#else
        __get_cpuid(i, (unsigned int*)&cpui[0], (unsigned int*)&cpui[1],
                    (unsigned int*)&cpui[2], (unsigned int*)&cpui[3]);
#endif
        data_.push_back(cpui);
      }

      // Capture vendor string
      char vendor[0x20];
      memset(vendor, 0, sizeof(vendor));
      *reinterpret_cast<int*>(vendor) = data_[0][1];
      *reinterpret_cast<int*>(vendor + 4) = data_[0][3];
      *reinterpret_cast<int*>(vendor + 8) = data_[0][2];
      vendor_ = vendor;

      if (vendor_ == "GenuineIntel")
      {
        isIntel_ = true;
      }
      else if (vendor_ == "AuthenticAMD")
      {
        isAMD_ = true;
      }

      // load bitset with flags for function 0x00000001
      if (nIds_ >= 1)
      {
        f_1_ECX_ = data_[1][2];
        f_1_EDX_ = data_[1][3];
      }

      // load bitset with flags for function 0x00000007
      if (nIds_ >= 7)
      {
        f_7_EBX_ = data_[7][1];
        f_7_ECX_ = data_[7][2];
      }

      // Calling __cpuid with 0x80000000 as the function_id argument
      // gets the number of the highest valid extended ID.
#ifdef WIN32
      __cpuid(cpui.data(), 0x80000000);
      nExIds_ = cpui[0];
#else
      nExIds_ = __get_cpuid_max(0x80000000, 0);
#endif

      char brand[0x40];
      memset(brand, 0, sizeof(brand));

      for (int i = 0x80000000; i <= nExIds_; ++i)
      {
#ifdef WIN32
        __cpuidex(cpui.data(), i, 0);
#else
        __get_cpuid(i, (unsigned int*)&cpui[0], (unsigned int*)&cpui[1],
                    (unsigned int*)&cpui[2], (unsigned int*)&cpui[3]);
#endif
        extdata_.push_back(cpui);
      }

      // load bitset with flags for function 0x80000001
      if (nExIds_ >= 0x80000001)
      {
        f_81_ECX_ = extdata_[1][2];
        f_81_EDX_ = extdata_[1][3];
      }

      // Interpret CPU brand string if reported
      if (nExIds_ >= 0x80000004)
      {
        memcpy(brand, extdata_[2].data(), sizeof(cpui));
        memcpy(brand + 16, extdata_[3].data(), sizeof(cpui));
        memcpy(brand + 32, extdata_[4].data(), sizeof(cpui));
        brand_ = brand;
      }

      // Check for multithreading feature like Hyperthreading.
      int info1[4];
#ifdef WIN32
      __cpuid(info1, 1);
#else
      __get_cpuid(1, (unsigned int*)&info1[0], (unsigned int*)&info1[1],
                  (unsigned int*)&info1[2], (unsigned int*)&info1[3]);
#endif
      extdata_.push_back(cpui);
      multithreading_ = (info1[3]& (1 << 28)) || false;
#endif
    };

    int nIds_;
    int nExIds_;
    std::string vendor_;
    std::string brand_;
    bool isIntel_;
    bool isAMD_;
    std::bitset<32> f_1_ECX_;
    std::bitset<32> f_1_EDX_;
    std::bitset<32> f_7_EBX_;
    std::bitset<32> f_7_ECX_;
    std::bitset<32> f_81_ECX_;
    std::bitset<32> f_81_EDX_;
    std::vector<std::array<int, 4>> data_;
    std::vector<std::array<int, 4>> extdata_;
    int cores_;
    bool multithreading_;
  };

  U64 NumberOfThreads_;
  U64 NumberOfCores_;
  U64 NumberOfLogicalCores_;
  const InstructionSet_Internal CPU_Rep_;
};

}  // namespace MTL

#endif  // MTL_CPU_H
