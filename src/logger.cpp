#include "logger.hpp"

#include <ctime>

Logger::Logger(const std::string &log_file_name) : log_output(log_file_name) {
  std::time_t t = std::time(nullptr);
  log_output << "---------" << std::put_time(std::localtime(&t), "%c %Z")
             << "---------" << std::endl;
}

Logger::Logger(const std::string &log_file_name, int argc, char *argv[]) : log_output(log_file_name) {
  for (uint16_t i = 0; i < argc; ++i)
    log_output << std::string(argv[i]) << " ";
  log_output << std::endl;

  std::time_t t = std::time(nullptr);
  log_output << "\n---------" << std::put_time(std::localtime(&t), "%c %Z")
             << "---------" << std::endl;
}

void Logger::write(const std::string &message) { log_output << message; }

void Logger::writeWithTime(const std::string &message) {
  std::time_t t = std::time(nullptr);
  log_output << std::put_time(std::localtime(&t), "[%H:%M:%S]") << " "
             << message;
}

void Logger::writeLineWithTime(const std::string &message) {
  std::time_t t = std::time(nullptr);
  log_output << std::put_time(std::localtime(&t), "[%H:%M:%S]") << " "
             << message << std::endl;
}

Logger::~Logger() { log_output.close(); }
