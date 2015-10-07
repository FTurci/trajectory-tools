#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <string>
#include <sstream>
#include <exception>
#include <type_traits>

template <typename T >
T StringToNum ( std::string &Text )
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}


/*** Variadic means of creating strings. ***/
inline std::string stringify(std::string arg)
{
    return arg;
}
inline std::string stringify(const char* arg)
{
    return std::string(arg);
}
// Conversions for numeric types.
template <class T,
          typename std::enable_if<std::is_arithmetic<T>::value, T>::type = 0 >
std::string stringify(T const& arg)
{
    return std::to_string(arg);
}
// The variadic list.
template<typename T, typename... Args>
std::string stringify(T const& first, Args const&... rest)
{
    std::string result = stringify( first );
    return result + stringify( rest... );
}


// Generic exception class which takes a variadic list of arguments and converts them to a string message using the above framework.
class Exception : public std::exception
{
public:
    template<typename... Args>
    Exception (Args const&... args) throw() : message( stringify(args...) )
    {
    }

    virtual const char* what() const throw()
    {
        return this->message.c_str();
    }

private:
    std::string message;
};


#endif
