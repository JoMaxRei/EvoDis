#include "exception.h"

Exception::Exception(const char *message, size_t error_code) : m_message(message), m_error_code(error_code)
{
}

const char *Exception::message()
{
    return m_message;
}

size_t Exception::error_code()
{
    return m_error_code;
}