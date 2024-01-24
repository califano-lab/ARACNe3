#include "stopwatch.hpp"

Watch::Watch() : zero(std::chrono::high_resolution_clock::now()) {}

void Watch::reset() { zero = std::chrono::high_resolution_clock::now(); }

std::string Watch::getSeconds() {
  auto cur = std::chrono::high_resolution_clock::now();
  return std::to_string(
             std::chrono::duration_cast<std::chrono::seconds>(cur - zero)
                 .count()) +
         "s";
}
