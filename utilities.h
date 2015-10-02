#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <string>
#include <sstream>  
#include <exception>

template <typename T >
T StringToNum ( std::string &Text )
{                               
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

inline std::string const& to_string(std::string const& s) { return s; }

template<typename... Args>
std::string stringify(Args const&... args)
{
    std::string result;
    using ::to_string;
    using std::to_string;
    int unpack[]{0, (result += to_string(args), 0)...};
    static_cast<void>(unpack);
    return result;
}

class Exception : public std::exception
{
public:
    template<typename... Args>
    Exception (Args const&... args) throw()
    {
        this->message =  stringify(args...);
    }
    
    virtual const char* what() const throw()
    {
        return this->message.c_str();
    }
    
private:
    std::string message;
};

#endif 

