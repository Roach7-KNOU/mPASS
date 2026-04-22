#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <string>

void readPdbFile(std::ifstream& ifs);
void readAtomPropertyFile(std::ifstream& ifs);
void readGridPropertyFile(std::ifstream& ifs);

inline std::string trimLeft(const std::string& str) {
    std::string::size_type n = str.find_first_not_of(" \t\v\n");
    return n == std::string::npos ? str : str.substr(n, str.length());
}

inline std::string trimRight(const std::string& str) {
    std::string::size_type n = str.find_last_not_of(" \t\v\n");
    return n == std::string::npos ? str : str.substr(0, n + 1);
}

inline std::string trim(const std::string& str) {
    return trimLeft(trimRight(str));
}

#endif
