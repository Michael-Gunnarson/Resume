#ifndef SanityFunctions_hpp
#define SanityFunctions_hpp

#include <iostream>
#include <string>
#include <vector>

#define PI 3.1415926535897932384626433832795029

void _print(std::string text);
void _print(int num);
void _print(bool num);
//void print(float num);
void _print(long double num);
void _print(std::vector<std::string> list);
void _print();

void _print_(std::string text); // print with no new line
void _print_(long double num);


#endif