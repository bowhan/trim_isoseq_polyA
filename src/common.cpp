#include <iostream>
#include "common.hpp"
#include "kernel_color.h"

void Warning(const std::string& s) {
    std::cerr << KERNAL_BOLDMAGENTA << "[Warning] " << s << KERNAL_RESET << std::endl;
}

void Error(const std::string& s) {
    std::cerr << KERNAL_BOLDRED << "[Error] " << s << KERNAL_RESET << std::endl;
    exit(EXIT_FAILURE);
}
