//
// Math Template Library
//
// Copyright (c) 2026: Seung Jae Lee, https://github.com/seungjlee/MTL
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

#include <MTL/Macros.h>
#include <MTL/Tools/Event.h>
#include <MTL/Tools/Signal.h>
#include <MTL/Tools/Test.h>

#include <atomic>
#include <string>
#include <thread>
#include <vector>

using namespace MTL;

TEST(TestSignalBasicEmission)
{
  Signal<void(int)> signal;

  std::vector<int> order;
  signal.Connect([&](int x) { order.push_back(x + 1); });
  signal.Connect([&](int x) { order.push_back(x + 2); });
  signal.Connect([&](int x) { order.push_back(x + 3); });

  MTL_EQUAL(signal.ConnectionCount(), size_t(3));

  signal(10);

  // Slots run in connection order.
  MTL_EQUAL(order.size(), size_t(3));
  MTL_EQUAL(order[0], 11);
  MTL_EQUAL(order[1], 12);
  MTL_EQUAL(order[2], 13);

  // Emit is an alias of operator().
  signal.Emit(20);
  MTL_EQUAL(order.size(), size_t(6));
  MTL_EQUAL(order[3], 21);
}

TEST(TestSignalArgumentsByReference)
{
  Signal<void(const std::string&, int&)> signal;

  signal.Connect([](const std::string& name, int& counter) { counter += (int)name.length(); });
  signal.Connect([](const std::string& name, int& counter) { counter += (int)name.length(); });

  int counter = 0;
  signal("four", counter);
  MTL_EQUAL(counter, 8);  // Both slots saw the same referenced counter.
}

TEST(TestSignalLastValueResult)
{
  Signal<int(int)> signal;

  // No slots: empty result.
  MTL_VERIFY(!signal(5).has_value());

  signal.Connect([](int x) { return x + 1; });
  signal.Connect([](int x) { return x * 10; });

  // The last connected slot's value wins (boost::signals2 last_value semantics).
  std::optional<int> result = signal(5);
  MTL_VERIFY(result.has_value());
  MTL_EQUAL(*result, 50);
}

TEST(TestSignalDisconnect)
{
  Signal<void()> signal;
  int calls1 = 0, calls2 = 0;

  SignalConnection connection1 = signal.Connect([&] { calls1++; });
  SignalConnection connection2 = signal.Connect([&] { calls2++; });
  MTL_VERIFY(connection1.IsConnected());
  MTL_VERIFY(connection2.IsConnected());

  signal();
  MTL_EQUAL(calls1, 1);
  MTL_EQUAL(calls2, 1);

  connection1.Disconnect();
  MTL_VERIFY(!connection1.IsConnected());
  MTL_VERIFY(connection2.IsConnected());
  MTL_EQUAL(signal.ConnectionCount(), size_t(1));

  signal();
  MTL_EQUAL(calls1, 1);  // No longer called.
  MTL_EQUAL(calls2, 2);

  // Double-disconnect is a safe no-op.
  connection1.Disconnect();
  signal();
  MTL_EQUAL(calls1, 1);
  MTL_EQUAL(calls2, 3);
}

TEST(TestSignalScopedConnection)
{
  Signal<void()> signal;
  int calls = 0;

  {
    ScopedSignalConnection scoped = signal.Connect([&] { calls++; });
    MTL_VERIFY(scoped.IsConnected());
    signal();
    MTL_EQUAL(calls, 1);
  }
  // Scope exit disconnected the slot.
  signal();
  MTL_EQUAL(calls, 1);
  MTL_EQUAL(signal.ConnectionCount(), size_t(0));

  // Move semantics: the moved-from connection must not disconnect.
  ScopedSignalConnection outer;
  {
    ScopedSignalConnection inner = signal.Connect([&] { calls++; });
    outer = std::move(inner);
  }
  signal();
  MTL_EQUAL(calls, 2);  // Still connected through 'outer'.
  outer.Disconnect();
  signal();
  MTL_EQUAL(calls, 2);

  // Release gives up ownership without disconnecting.
  SignalConnection kept;
  {
    ScopedSignalConnection scoped = signal.Connect([&] { calls++; });
    kept = scoped.Release();
  }
  signal();
  MTL_EQUAL(calls, 3);
  kept.Disconnect();
}

TEST(TestSignalConnectionOutlivesSignal)
{
  SignalConnection connection;
  {
    Signal<void()> signal;
    connection = signal.Connect([] {});
    MTL_VERIFY(connection.IsConnected());
  }
  // The signal is gone; the handle must be inert, not crash.
  MTL_VERIFY(!connection.IsConnected());
  connection.Disconnect();
}

TEST(TestSignalConnectDuringEmission)
{
  Signal<void()> signal;
  int outerCalls = 0, innerCalls = 0;

  signal.Connect([&] {
    outerCalls++;
    if (outerCalls == 1)
      signal.Connect([&] { innerCalls++; });
  });

  signal();
  // The slot connected during the emission must not run in that emission.
  MTL_EQUAL(outerCalls, 1);
  MTL_EQUAL(innerCalls, 0);

  signal();
  MTL_EQUAL(outerCalls, 2);
  MTL_EQUAL(innerCalls, 1);
}

TEST(TestSignalSelfDisconnectDuringEmission)
{
  Signal<void()> signal;
  int calls = 0;

  SignalConnection connection;
  connection = signal.Connect([&] {
    calls++;
    connection.Disconnect();  // Self-disconnect from inside the slot: must not deadlock.
  });

  signal();
  MTL_EQUAL(calls, 1);
  signal();
  MTL_EQUAL(calls, 1);
  MTL_EQUAL(signal.ConnectionCount(), size_t(0));
}

TEST(TestSignalDisconnectLaterSlotDuringEmission)
{
  Signal<void()> signal;
  int calls1 = 0, calls2 = 0;

  SignalConnection connection2;
  signal.Connect([&] {
    calls1++;
    connection2.Disconnect();  // Disconnect a slot that has not yet run in this emission.
  });
  connection2 = signal.Connect([&] { calls2++; });

  signal();
  MTL_EQUAL(calls1, 1);
  MTL_EQUAL(calls2, 0);  // Must not run: it was disconnected before its turn.
}

TEST(TestSignalRecursiveEmission)
{
  Signal<void(int)> signal;
  std::vector<int> values;

  signal.Connect([&](int depth) {
    values.push_back(depth);
    if (depth > 0)
      signal(depth - 1);  // Recursive emission of the same signal from inside a slot.
  });

  signal(3);
  MTL_EQUAL(values.size(), size_t(4));
  MTL_EQUAL(values[0], 3);
  MTL_EQUAL(values[3], 0);
}

// The race this class exists to prevent: after Disconnect returns on one thread, the
// slot must not still be executing on another. (With boost::signals2, disconnect does
// not wait for in-flight invocations - the cause of use-after-scope bugs in handlers.)
TEST(TestSignalDisconnectBlocksUntilSlotCompletes)
{
  Signal<void()> signal;

  Event slotStarted;
  std::atomic<bool> insideSlot(false);
  std::atomic<bool> tornStateObserved(false);

  SignalConnection connection = signal.Connect([&] {
    insideSlot = true;
    slotStarted.Signal();
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    insideSlot = false;
  });

  std::thread emitter([&] { signal(); });

  MTL_VERIFY(slotStarted.Wait(5000));  // The slot is now running on the emitter thread.
  connection.Disconnect();             // Must block until the slot finishes.
  if (insideSlot)
    tornStateObserved = true;

  emitter.join();
  MTL_VERIFY(!tornStateObserved);
}

TEST(TestSignalSameSlotSerialized)
{
  // Concurrent emissions must not run the same slot in parallel.
  Signal<void()> signal;

  std::atomic<int> concurrent(0);
  std::atomic<int> maxConcurrent(0);
  std::atomic<int> calls(0);

  signal.Connect([&] {
    int now = ++concurrent;
    int expectedMax = maxConcurrent.load();
    while (now > expectedMax && !maxConcurrent.compare_exchange_weak(expectedMax, now)) {}
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    --concurrent;
    ++calls;
  });

  static const int kThreads = 4;
  static const int kEmitsPerThread = VALUE_DEBUG_RELEASE(5, 25);
  std::vector<std::thread> threads;
  for (int t = 0; t < kThreads; t++)
    threads.emplace_back([&] {
      for (int i = 0; i < kEmitsPerThread; i++)
        signal();
    });
  for (std::thread& thread : threads)
    thread.join();

  MTL_EQUAL(calls.load(), kThreads * kEmitsPerThread);
  MTL_EQUAL(maxConcurrent.load(), 1);
}

TEST(TestSignalConcurrentStress)
{
  // Emitters, connectors, and disconnectors all running at once: must not crash,
  // deadlock, or lose the accounting invariants checked below.
  Signal<void(int)> signal;

  std::atomic<long long> receivedSum(0);
  std::atomic<bool> stop(false);

  // A permanent slot so emissions always have at least one receiver.
  signal.Connect([&](int x) { receivedSum += x; });

  static const int kEmitThreads = 2;
  static const int kEmitsPerThread = VALUE_DEBUG_RELEASE(500, 5000);

  std::vector<std::thread> threads;
  std::atomic<long long> emittedSum(0);
  for (int t = 0; t < kEmitThreads; t++)
    threads.emplace_back([&] {
      for (int i = 1; i <= kEmitsPerThread; i++)
      {
        signal(i);
        emittedSum += i;
      }
    });

  // Churn: connect + disconnect transient slots while emissions are running.
  threads.emplace_back([&] {
    while (!stop)
    {
      std::atomic<int> transientCalls(0);
      SignalConnection connection = signal.Connect([&](int) { transientCalls++; });
      std::this_thread::sleep_for(std::chrono::microseconds(50));
      connection.Disconnect();
      // The blocking guarantee makes it safe for transientCalls to leave scope now.
    }
  });

  for (int t = 0; t < kEmitThreads; t++)
    threads[t].join();
  stop = true;
  threads.back().join();

  // The permanent slot saw every emission exactly once.
  MTL_EQUAL(receivedSum.load(), emittedSum.load());
  MTL_EQUAL(signal.ConnectionCount(), size_t(1));
}

TEST(TestSignalDisconnectAllAndDestruction)
{
  int calls = 0;
  Signal<void()> signal;
  signal.Connect([&] { calls++; });
  signal.Connect([&] { calls++; });

  signal();
  MTL_EQUAL(calls, 2);

  signal.DisconnectAll();
  MTL_EQUAL(signal.ConnectionCount(), size_t(0));
  signal();
  MTL_EQUAL(calls, 2);

  // Destruction with live connections and a running emitter thread.
  Event slotRunning;
  std::atomic<bool> slotFinished(false);
  std::thread emitter;
  {
    Signal<void()> dying;
    dying.Connect([&] {
      slotRunning.Signal();
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
      slotFinished = true;
    });
    emitter = std::thread([&] { dying(); });
    MTL_VERIFY(slotRunning.Wait(5000));
    // ~Signal now runs: it must block until the in-flight slot completes.
  }
  MTL_VERIFY(slotFinished.load());
  emitter.join();
}
