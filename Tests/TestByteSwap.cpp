#include <MTL/Tools/Test.h>
#include <MTL/ByteSwap.h>

#include <numeric>
#include <vector>

using namespace MTL;

template<typename T> static void TestSwap(T bytes)
{
  T swapped = ByteSwap(bytes);
  const uint8_t* input = (uint8_t*)&bytes;
  const uint8_t* output = (uint8_t*)&swapped;
  for (uint32_t i = 0; i < sizeof(T); i++)
    MTL_EQUAL(output[i], input[sizeof(T) - 1 - i]);
}

TEST(TestByteSwap)
{
  std::vector<uint8_t> bytes(8);
  std::iota(bytes.begin(), bytes.end(), (uint8_t)1);

  TestSwap(*(int16_t*)bytes.data());
  TestSwap(*(int32_t*)bytes.data());
  TestSwap(*(int64_t*)bytes.data());
}
