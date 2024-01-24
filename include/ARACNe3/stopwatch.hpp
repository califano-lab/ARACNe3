#pragma once

#include <chrono>
#include <string>

class Stopwatch {
public:
  Stopwatch();
  void reset();
  std::string getSeconds();

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> zero;
};
