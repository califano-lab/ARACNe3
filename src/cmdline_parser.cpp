#include "cmdline_parser.hpp"

CmdLineParser::CmdLineParser(int argc, char **argv) {
  this->begin = argv;
  this->end = argv+argc;
  this->argc = argc;
}


char *CmdLineParser::getOpt(const std::string& option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    return *itr;
  return nullptr;
}

bool CmdLineParser::optExists(const std::string& option) {
  return std::find(begin, end, option) != end;
}
