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


#ifndef MTL_SPARSE_H
#define MTL_SPARSE_H

#include <MTL/Davis_LDL_COLAMD.h>

namespace MTL
{

template <class T> class CompressedSparseMatrix;

// Helper class for LDLt solver.
template <class T = F64>
class SymbolicLDLt
{
public:
  SymbolicLDLt(const CompressedSparseMatrix<T>& A)
  {
    assert(A.Cols() == A.Rows());

    I32 n = A.Cols();

    const I32* Ap = A.Ap();
    const I32* Ai = A.Ai();

    I32 stats[COLAMD_STATS];
    P_.Resize(n + 1);
    symamd(n, Ai, Ap, P_.begin(), NULL, stats, &calloc, &free);

    buffer_.Resize(5, n);
    Parent_  = buffer_[0];
    Lnz_     = buffer_[1];
    Flag_    = buffer_[2];
    Pattern_ = buffer_[3];
    Pinv_    = buffer_[4];

    L_ = CompressedSparseMatrix<T>(n, n, n + 1);

    for (int i = 0; i < n; i++)
      Pinv_[P_[i]] = i;

    LDLt_Symbolic(n, Ap, Ai, L_.Ap(), Parent_, Lnz_, Flag_, P_.Begin(), Pinv_);

    L_.Allocate(L_.Ap()[n]);
  }

  CompressedSparseMatrix<T>& L()  { return L_; }
  const DynamicVector<I32>& P() const  { return P_; }
  const I32* Parent() const  { return Parent_; }
  const I32* Pinv() const    { return Pinv_;   }
  I32* Lnz()      { return Lnz_;  }
  I32* Flag()     { return Flag_; }
  I32* Pattern()  { return Pattern_; }

private:
  CompressedSparseMatrix<T> L_;
  DynamicVector<I32> P_;
  I32* Parent_;
  I32* Lnz_;
  I32* Flag_;
  I32* Pattern_;
  I32* Pinv_;

  DynamicMatrix<I32> buffer_;
};

}  // namespace MTL

#endif // MTL_SPARSE_H
