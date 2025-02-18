//
// Math Template Library
//
// Copyright (c) 2025: Seung Jae Lee, https://github.com/seungjlee/MTL
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
#include <MTL/MemoryPool.h>

using namespace MTL;

TEST(TestMemoryPoolBasic)
{
    MemoryPool pool(1024);  // 1KB pool
    
    // Test single allocation
    void* p1 = pool.allocate(128);
    MTL_VERIFY(p1 != nullptr);
    
    // Test deallocation and reuse
    pool.deallocate(p1);
    void* p2 = pool.allocate(128);
    MTL_VERIFY(p2 == p1);  // Should reuse same block
}

TEST(TestMemoryPoolMultiple)
{
    MemoryPool pool(1024);
    
    void* blocks[4];
    // Allocate multiple blocks
    for(int i = 0; i < 4; i++) {
        blocks[i] = pool.allocate(128);
        MTL_VERIFY(blocks[i] != nullptr);
    }
    
    // Verify blocks don't overlap
    for(int i = 0; i < 3; i++) {
        MTL_VERIFY((char*)blocks[i+1] >= ((char*)blocks[i] + 128));
    }
}

TEST(TestMemoryPoolMergeRegions1)
{
    PoolBlock block(1024);
    
    // Allocate 3 regions of 256 bytes each
    void* p1 = block.allocate(256);
    void* p2 = block.allocate(256);
    void* p3 = block.allocate(256);
    MTL_VERIFY(p1 != nullptr);
    MTL_VERIFY(p2 != nullptr);
    MTL_VERIFY(p3 != nullptr);
    
    // Free p1 and p2 - should merge
    block.deallocate(p1);
    block.deallocate(p2);
    
    // Try to allocate 512 bytes - should succeed using merged space
    void* large = block.allocate(512);
    MTL_VERIFY(large == p1);
    
    // Try to allocate remaining 256 - should fail
    void* tooLarge = block.allocate(512);
    MTL_VERIFY(tooLarge == nullptr);
    
    // Free everything
    block.deallocate(large);
    block.deallocate(p3);
    
    // Should be able to allocate full size
    void* full = block.allocate(1024);
    MTL_VERIFY(full == p1);
}

TEST(TestMemoryPoolMergeRegions2)
{
    PoolBlock block(1024);
    
    // Test sequential allocations
    void* regions[4];
    for(int i = 0; i < 4; i++) {
        regions[i] = block.allocate(128);
        MTL_VERIFY(regions[i] != nullptr);
    }
    
    // Free middle regions - should merge
    block.deallocate(regions[1]);
    block.deallocate(regions[2]);
    
    // Try to allocate merged space
    void* merged = block.allocate(256);
    MTL_VERIFY(merged == regions[1]);
    
    // Free all regions
    block.deallocate(regions[0]);
    block.deallocate(merged);
    block.deallocate(regions[3]);
    
    // Should be able to allocate full block
    void* full = block.allocate(1024);
    MTL_VERIFY(full == regions[0]);
    MTL_VERIFY(full != nullptr);
}