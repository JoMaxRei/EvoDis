#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <string>
#include <sstream>

class Exception
{
public:
    Exception(char *message, int error_code = -1);
    char *message();
    int error_code();

private:
    char *m_message;
    int m_error_code;
};

#endif