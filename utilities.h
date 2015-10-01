#ifndef __utilities_H
#define __utilities_H

//#include <iostream>
//#include <fstream>
//#include <vector>
#include <string>
#include <sstream>  


template <typename T >
T StringToNum ( std::string &Text )
{                               
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

inline std::string const& to_string(std::string const& s) { return s; }

template<typename... Args>
std::string stringer(Args const&... args)
{
    std::string result;
    using ::to_string;
    using std::to_string;
    int unpack[]{0, (result += to_string(args), 0)...};
    static_cast<void>(unpack);
    return result;
}

#endif 

