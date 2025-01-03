#include "exception.h"

Exception::Exception(char *message, int error_code) : m_message(message), m_error_code(error_code)
{
}

char *Exception::message()
{
    return m_message;
}

int Exception::error_code()
{
    return m_error_code;
}