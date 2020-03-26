#include "my_assert.h"

string getFilename(string filePath, bool withExtension, char seperator) {
    // Get last dot position
    size_t dotPos = filePath.rfind('.');
    size_t sepPos = filePath.rfind(seperator);

    if (sepPos != string::npos) {
        return filePath.substr(
            sepPos + 1,
            filePath.size() -
                (withExtension || dotPos != std::string::npos ? 1 : dotPos));
    }
    return "";
}


my_invalid_argument::my_invalid_argument(const std::string& arg,
                                         const char* file, int line)
    : std::runtime_error(arg) {

    std::ostringstream o;
    o << arg << " At " << getFilename(file) << ":" << line;
    msg = o.str();
}

my_runtime_error::my_runtime_error(const std::string& arg, const char* file,
                                   int line)
    : std::runtime_error(arg) {

    std::ostringstream o;
    o << arg << " At " << getFilename(file) << ":" << line;
    msg = o.str();
}
