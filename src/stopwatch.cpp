#include "stopwatch.hpp"

Stopwatch::Stopwatch() : zero(std::chrono::high_resolution_clock::now()) {}

void Stopwatch::reset() { zero = std::chrono::high_resolution_clock::now(); }

std::string Stopwatch::getSeconds() {
  auto cur = std::chrono::high_resolution_clock::now();
  return std::to_string(
             std::chrono::duration_cast<std::chrono::seconds>(cur - zero)
                 .count()) +
         "s";
}
