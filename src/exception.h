#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <string>
#include <sstream>

class Exception
{
public:
    Exception(const char *message, size_t error_code = -1);
    const char *message();
    size_t error_code();

private:
    const char *m_message;
    size_t m_error_code;
};

#endif