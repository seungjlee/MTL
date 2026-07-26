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

#ifndef MTL_SIGNAL_H
#define MTL_SIGNAL_H

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <vector>

namespace MTL
{

//
// MTL::Signal<R(Args...)> - a thread-safe multicast callback (observer) with explicit,
// documented concurrency semantics:
//
//   * Connect/Disconnect/Emit are safe to call concurrently from any thread.
//   * Slots are invoked in connection order, outside of the signal's internal list lock,
//     so slots may freely Connect or Disconnect slots of the same signal.
//   * A slot connected during an emission is not invoked by that emission.
//   * A slot disconnected during an emission is not invoked once Disconnect returns.
//   * Disconnect BLOCKS until an invocation of that slot running on another thread has
//     completed. After Disconnect returns, the slot is guaranteed to not be executing
//     and to never execute again. (This closes the classic use-after-scope race where a
//     handler's owner is destroyed while the handler is still running.)
//   * A slot disconnecting ITSELF from within its own invocation is detected and returns
//     without blocking.
//   * Emissions of the same slot are serialized: the same slot is never run concurrently
//     from multiple emitting threads. Distinct slots may run concurrently.
//   * Recursive emission of the same signal from within a slot is supported.
//   * For a non-void result type, Emit returns the value of the last slot invoked wrapped
//     in std::optional (empty when no slot ran).
//
// Constraint: do not block-Disconnect slot B from inside slot A of the same signal while
// another thread inside slot B may block-Disconnect slot A; as with any blocking
// synchronization, opposing waits can deadlock. Self-disconnection is always safe.
//
// Exceptions thrown by a slot propagate to the emitter; remaining slots are not invoked
// for that emission.
//

namespace SignalImplementation
{

struct SlotControl
{
  std::atomic<bool> Connected{true};
  std::recursive_mutex RunMutex;                 // Held while the slot executes.
  std::atomic<std::thread::id> RunningThread{};  // Thread currently executing the slot.

  // Marks the slot disconnected. Blocks until an in-flight invocation on another
  // thread has completed. Non-blocking when called from within the slot itself.
  void Disconnect()
  {
    Connected = false;
    if (RunningThread.load() == std::this_thread::get_id())
      return;  // Self-disconnection from inside the slot: safe, do not block.
    std::lock_guard<std::recursive_mutex> waitForCompletion(RunMutex);
  }
};

// Restores the running-thread marker even if the slot throws.
class RunningThreadScope
{
public:
  explicit RunningThreadScope(std::atomic<std::thread::id>& runningThread)
    : RunningThread_(runningThread), Previous_(runningThread.exchange(std::this_thread::get_id()))
  {
  }
  ~RunningThreadScope() { RunningThread_.store(Previous_); }

private:
  std::atomic<std::thread::id>& RunningThread_;
  std::thread::id Previous_;
};

}  // namespace SignalImplementation

// Handle to a single connection. Copyable; all copies refer to the same connection.
// Outliving the signal is safe: operations on an expired connection are no-ops.
class SignalConnection
{
public:
  SignalConnection() = default;
  explicit SignalConnection(std::shared_ptr<SignalImplementation::SlotControl> control)
    : Control_(std::move(control))
  {
  }

  // See SlotControl::Disconnect for the blocking guarantee.
  void Disconnect()
  {
    if (auto control = Control_.lock())
      control->Disconnect();
  }

  bool IsConnected() const
  {
    auto control = Control_.lock();
    return control && control->Connected;
  }

private:
  std::weak_ptr<SignalImplementation::SlotControl> Control_;
};

// RAII connection: disconnects (with the blocking guarantee) when destroyed.
class ScopedSignalConnection
{
public:
  ScopedSignalConnection() = default;
  ScopedSignalConnection(SignalConnection connection)  // Implicit from Connect result.
    : Connection_(std::move(connection))
  {
  }
  ScopedSignalConnection(ScopedSignalConnection&& other) noexcept
    : Connection_(other.Release())
  {
  }
  ScopedSignalConnection& operator=(ScopedSignalConnection&& other) noexcept
  {
    if (this != &other)
    {
      Disconnect();
      Connection_ = other.Release();
    }
    return *this;
  }
  ScopedSignalConnection(const ScopedSignalConnection&) = delete;
  ScopedSignalConnection& operator=(const ScopedSignalConnection&) = delete;

  ~ScopedSignalConnection() { Disconnect(); }

  void Disconnect()
  {
    Connection_.Disconnect();
    Connection_ = SignalConnection();
  }

  // Gives up ownership without disconnecting.
  SignalConnection Release()
  {
    SignalConnection connection = Connection_;
    Connection_ = SignalConnection();
    return connection;
  }

  bool IsConnected() const { return Connection_.IsConnected(); }

private:
  SignalConnection Connection_;
};

template <class Signature>
class Signal;

template <class R, class... Args>
class Signal<R(Args...)>
{
public:
  using SlotFunction = std::function<R(Args...)>;

  Signal() = default;
  Signal(const Signal&) = delete;
  Signal& operator=(const Signal&) = delete;

  // Destroying the signal disconnects all slots with the blocking guarantee.
  ~Signal() { DisconnectAll(); }

  SignalConnection Connect(SlotFunction slot)
  {
    auto entry = std::make_shared<Entry>(std::move(slot));
    std::lock_guard<std::mutex> lock(Mutex_);
    Slots_.push_back(entry);
    // Aliasing constructor: the handle keeps the control block reachable while
    // sharing the entry's lifetime.
    return SignalConnection(
      std::shared_ptr<SignalImplementation::SlotControl>(entry, &entry->Control));
  }

  // Invokes all connected slots in connection order.
  // For non-void R, returns the last invoked slot's result (empty if none ran).
  auto operator()(Args... args)
  {
    std::vector<std::shared_ptr<Entry>> snapshot = Snapshot();

    if constexpr (std::is_void_v<R>)
    {
      for (const auto& entry : snapshot)
        Invoke(*entry, args...);
    }
    else
    {
      std::optional<R> lastResult;
      for (const auto& entry : snapshot)
      {
        std::optional<R> result = Invoke(*entry, args...);
        if (result)
          lastResult = std::move(result);
      }
      return lastResult;
    }
  }

  auto Emit(Args... args) { return (*this)(args...); }

  // Number of currently connected slots.
  size_t ConnectionCount() const
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    size_t count = 0;
    for (const auto& entry : Slots_)
      if (entry->Control.Connected)
        count++;
    return count;
  }

  // Disconnects every slot with the same blocking guarantee as SignalConnection::Disconnect.
  void DisconnectAll()
  {
    // Note: identifiers named "slots", "signals" or "emit" must not be used in this
    // header - Qt defines them as macros, which would break any Qt translation unit
    // that includes it.
    std::vector<std::shared_ptr<Entry>> allSlots;
    {
      std::lock_guard<std::mutex> lock(Mutex_);
      allSlots.swap(Slots_);
    }
    for (const auto& entry : allSlots)
      entry->Control.Disconnect();
  }

private:
  struct Entry
  {
    explicit Entry(SlotFunction&& function) : Function(std::move(function)) {}
    SignalImplementation::SlotControl Control;
    SlotFunction Function;
  };

  // Copies the current slot list and prunes disconnected entries.
  std::vector<std::shared_ptr<Entry>> Snapshot()
  {
    std::lock_guard<std::mutex> lock(Mutex_);
    size_t connected = 0;
    for (auto& entry : Slots_)
      if (entry->Control.Connected)
        Slots_[connected++] = std::move(entry);
    Slots_.resize(connected);
    return Slots_;
  }

  // Runs one slot with the running-thread bookkeeping. Returns the result for non-void R,
  // wrapped in std::optional (empty when the slot was skipped as disconnected).
  auto Invoke(Entry& entry, Args&... args)
  {
    if constexpr (std::is_void_v<R>)
    {
      if (!entry.Control.Connected)
        return;
      std::lock_guard<std::recursive_mutex> run(entry.Control.RunMutex);
      if (!entry.Control.Connected)  // Re-check: may have disconnected while we waited.
        return;
      SignalImplementation::RunningThreadScope scope(entry.Control.RunningThread);
      entry.Function(args...);
    }
    else
    {
      std::optional<R> result;
      if (!entry.Control.Connected)
        return result;
      std::lock_guard<std::recursive_mutex> run(entry.Control.RunMutex);
      if (!entry.Control.Connected)
        return result;
      SignalImplementation::RunningThreadScope scope(entry.Control.RunningThread);
      result = entry.Function(args...);
      return result;
    }
  }

  mutable std::mutex Mutex_;
  std::vector<std::shared_ptr<Entry>> Slots_;
};

}  // namespace MTL

#endif  // MTL_SIGNAL_H
