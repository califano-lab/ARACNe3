#pragma once

#include <string>

class CmdLineParser {
private:
  char **begin, **end;
  int argc;

public:
  CmdLineParser(int argc, char **argv);
  char *getOpt(const std::string& option);
  bool optExists(const std::string& option);
};
