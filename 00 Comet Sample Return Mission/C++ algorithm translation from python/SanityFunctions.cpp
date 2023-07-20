#include "SanityFunctions.hpp"

void print(std::string text) {
    std::cout << text << "\n";
}

void print(int num) {
    std::cout << num << "\n";
}

void print(bool num) {
    std::cout << num << "\n";
}

/*void print(float num) {
    std::cout << num << "\n";
}*/

void print(long double num) {
    std::cout << num << "\n";
}

void print(std::vector<std::string> list) {
    int n = list.size();
    for(int i=0; i<n; i++) {
        std::cout << list.at(i) << "\n";
    }
}

void print(std::vector<int> list) {
    int n = list.size();
    for(int i=0; i<n; i++) {
        std::cout << list.at(i) << "\n";
    }
}

void printN(std::string text) {
    std::cout << text;
}