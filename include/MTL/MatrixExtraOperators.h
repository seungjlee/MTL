//
// Math Template Library
//
// Copyright (c) 2018: Seung Jae Lee, https://github.com/seungjlee/MTL
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


#ifndef MTL_MATRIX_EXTRA_OPERATORS_H
#define MTL_MATRIX_EXTRA_OPERATORS_H

#include <MTL/ArrayExtra.h>
#include <MTL/Matrix.h>

template<MTL::I32 M, MTL::I32 N, class T>
MTL_INLINE MTL::Matrix<M,N,T>& operator^=(MTL::Matrix<M,N,T>& dst, const MTL::Matrix<M,N,T>& src)
{
  MTL::ArrayExtra<M*N,T>::XOR(dst.Data()[0], src.Data()[0]);
  return dst;
}
template<MTL::I32 M, MTL::I32 N, class T>
MTL_INLINE MTL::Matrix<M,N,T> operator^(const MTL::Matrix<M,N,T>& a, const MTL::Matrix<M,N,T>& b)
{
  MTL::Matrix c = a;
  c += b;
  return c;
}

#endif // MTL_MATRIX_EXTRA_OPERATORS_H
