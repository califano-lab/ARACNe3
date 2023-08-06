#include "stopwatch.hpp"

Watch::Watch() : zero(std::chrono::high_resolution_clock::now()) {}

void Watch::setWatch() { zero = std::chrono::high_resolution_clock::now(); }

void Watch::printWatch(std::ostream &ofs, const std::string message) {
  auto cur = std::chrono::high_resolution_clock::now();
  ofs << message;
  ofs << std::chrono::duration_cast<std::chrono::seconds>(cur - zero).count()
      << "s" << std::endl;
}
