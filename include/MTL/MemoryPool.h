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

#ifndef MTL_MEMORY_POOL_H
#define MTL_MEMORY_POOL_H

#include <vector>
#include <list>
#include <mutex>

namespace MTL
{

struct MemoryRegion
{
    size_t offset;
    size_t size;
    bool used;
};

class PoolBlock
{
    std::vector<uint8_t> data_;
    std::list<MemoryRegion> regions_;
    
public:
    explicit PoolBlock(size_t size) : data_(size)
    {
        regions_.push_back({0, size, false});
    }
    
    void* allocate(size_t size)
    {
        for (auto it = regions_.begin(); it != regions_.end(); ++it) {
            if (!it->used && it->size >= size) {
                if (it->size > size) {
                    regions_.insert(std::next(it), 
                        {it->offset + size, it->size - size, false});
                }
                it->size = size;
                it->used = true;
                return data_.data() + it->offset;
            }
        }
        return nullptr;
    }
    
    bool deallocate(void* ptr)
    {
        size_t offset = static_cast<uint8_t*>(ptr) - data_.data();
        for (auto it = regions_.begin(); it != regions_.end(); ++it) {
            if (it->offset == offset && it->used) {
                it->used = false;
                merge_regions();
                return true;
            }
        }
        return false;
    }

private:
    void merge_regions()
    {
        auto it = regions_.begin();
        while (it != regions_.end()) {
            auto next = std::next(it);
            if (next != regions_.end() && !it->used && !next->used) {
                it->size += next->size;
                regions_.erase(next);
            } else {
                ++it;
            }
        }
    }
};

class MemoryPool
{
    std::list<PoolBlock> blocks_;
    std::mutex mutex_;
    size_t blockSize_;

public:
    explicit MemoryPool(int initialBlocks = 1, size_t blockSize = 1024*1024)
      : blockSize_(blockSize)
    {
        for (int i = 0; i < initialBlocks; ++i) {
            blocks_.emplace_back(blockSize_);
        }
    }

    void* allocate(size_t size)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& block : blocks_) {
            if (void* ptr = block.allocate(size))
                return ptr;
        }
        blocks_.emplace_back(std::max(size, blockSize_));
        return blocks_.back().allocate(size);
    }

    void deallocate(void* ptr)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        for (auto& block : blocks_) {
            if (block.deallocate(ptr))
                return;
        }
    }
};

}  // namespace MTL

#endif  // MTL_MEMORY_POOL_H