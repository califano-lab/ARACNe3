#ifndef WATCH_H
#define WATCH_H

#include <chrono>
#include <ostream>
#include <string>

// TODO: Add Documentation
class Watch {
public:
  Watch();
  void reset();
  std::string getSeconds();

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> zero;
};

#endif // WATCH_H
