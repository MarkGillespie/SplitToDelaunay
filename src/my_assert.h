#pragma once

#include <iostream>
#include <sstream>
#include <string>

using std::string;

string getFilename(string filePath, bool withExtension = true,
                   char seperator = '/');

class my_invalid_argument : public std::runtime_error {
    std::string msg;

  public:
    my_invalid_argument(const std::string& arg, const char* file, int line);
    ~my_invalid_argument() throw() {}
    const char* what() const throw() { return msg.c_str(); }
};

class my_runtime_error : public std::runtime_error {
    std::string msg;

  public:
    my_runtime_error(const std::string& arg, const char* file, int line);
    ~my_runtime_error() throw() {}
    const char* what() const throw() { return msg.c_str(); }
};

#define non_nan(num)                                                           \
    if (num != num) {                                                          \
        throw_my_runtime_error("Encountered NaN.");                            \
    }

#define happy(num)                                                             \
    if (num != num) {                                                          \
        throw_my_runtime_error("Encountered NaN.");                            \
    } else if (!std::isfinite(num)) {                                          \
        throw_my_runtime_error("Encountered Inf.");                            \
    }

#define my_assert(arg, msg)                                                    \
    if (!(arg)) throw_my_runtime_error(msg);

#define throw_my_invalid_argument(arg)                                         \
    throw my_invalid_argument(arg, __FILE__, __LINE__);

#define throw_my_runtime_error(arg)                                            \
    throw my_runtime_error(arg, __FILE__, __LINE__);
