#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <fstream>
#include <string>

/**
 * Logger class for writing messages to a log file.
 */
class Logger {
public:
  /**
   * Constructor that opens a log file and writes the time.
   *
   * @param log_file_name Path to the log file.
   */
  Logger(const std::string &log_file_name);

  /**
   * Constructor that opens a log file and prints the initial command (it is
   * assumed that argc and argv are derived from main). Then writes the time.
   *
   * @param log_file_name Path to the log file.
   * @param argc argc from main().
   * @param argv argv from main().
   */
  Logger(const std::string &log_file_name, int argc, char *argv[]);

  /**
   * Function to write a message to the log.
   *
   * @param message The message to be logged.
   */
  void write(const std::string &message);

  /**
   * Writes a message line to the log.
   *
   * @param message The message to be logged.
   */
  void writeLine(const std::string &message);

  /**
   * Writes a message to the log, prepended with a timestamp.
   *
   * @param message The message to be logged.
   */
  void writeWithTime(const std::string &message);

  /**
   * Writes a message line to the log, prepended with a timestamp.
   *
   * @param message The message to be logged.
   */
  void writeLineWithTime(const std::string &message);

  /**
   * Destructor that closes the log file.
   */
  ~Logger();

private:
  std::ofstream log_output;
};

#endif // LOGGER_HPP
