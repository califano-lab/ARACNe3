#ifndef WATCH_H
#define WATCH_H

#include <chrono>
#include <string>
#include <ostream>

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
