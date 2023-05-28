#pragma once

#include <vector>

template<typename T>
struct Triag_Element{
    public:

    std::vector<T> x;
    std::vector<T> y;
    std::vector<std::size_t> poly;

};