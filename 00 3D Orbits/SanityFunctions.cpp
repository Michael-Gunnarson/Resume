#include "SanityFunctions.hpp"

void _print(std::string text) {
    std::cout << text << "\n";
}

void _print(int num) {
    std::cout << num << "\n";
}

void _print(bool num) {
    std::cout << num << "\n";
}

/*void _print(float num) {
    std::cout << num << "\n";
}*/

void _print(long double num) {
    std::cout << num << "\n";
}

void _print(std::vector<std::string> list) {
    int n = list.size();
    for(int i=0; i<n; i++) {
        std::cout << list.at(i) << "\n";
    }
}

void _print(std::vector<int> list) {
    int n = list.size();
    for(int i=0; i<n; i++) {
        std::cout << list.at(i) << "\n";
    }
}

void _print() {
    std::cout << "\n";
}

void _print_(std::string text) {
    std::cout << text;
}

void _print_(long double num) {
    std::cout << num;
}