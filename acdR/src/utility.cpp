
#include "utility.hpp"

// find numbers in a string and return the numbers in a vector of int
std::vector<int> extract_integers(std::string const & str)
{
    std::vector<int> results;
    const std::regex regex(R"(\d+)");   // matches a sequence of digits
    std::smatch match;
    std::string s(str);

    for( std::smatch match; std::regex_search(s, match, regex); )
    {
        results.emplace_back( std::stoi(match.str()) );
        s = match.suffix();
    }

    return results;
}
