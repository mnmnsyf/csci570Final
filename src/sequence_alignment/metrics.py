"""Runtime and memory measurement helpers."""

import os
import threading
import time
from dataclasses import dataclass

import psutil


@dataclass(frozen=True)
class RunMetrics:
    elapsed_seconds: float
    memory_before_bytes: int
    memory_after_bytes: int
    peak_memory_bytes: int

    @property
    def elapsed_ms(self):
        return self.elapsed_seconds * 1000

    @property
    def peak_memory_kb(self):
        return self.peak_memory_bytes / 1024


def run_with_metrics(operation, sample_interval_seconds=0.001):
    process = psutil.Process(os.getpid())
    memory_before = process.memory_info().rss
    peak_memory = memory_before
    stop_event = threading.Event()

    def sample_memory():
        nonlocal peak_memory
        while not stop_event.is_set():
            peak_memory = max(peak_memory, process.memory_info().rss)
            stop_event.wait(sample_interval_seconds)

    sampler = threading.Thread(target=sample_memory, daemon=True)
    start_time = time.perf_counter()
    sampler.start()
    try:
        value = operation()
    finally:
        elapsed_seconds = time.perf_counter() - start_time
        stop_event.set()
        sampler.join()

    memory_after = process.memory_info().rss
    peak_memory = max(peak_memory, memory_after)
    return value, RunMetrics(elapsed_seconds, memory_before, memory_after, peak_memory)
