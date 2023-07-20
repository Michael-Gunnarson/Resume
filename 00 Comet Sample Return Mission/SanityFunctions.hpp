#ifndef SanityFunctions_hpp
#define SanityFunctions_hpp

#include <iostream>
#include <string>
#include <vector>

void print(std::string text);
void print(int num);
void print(bool num);
//void print(float num);
void print(long double num);
void print(std::vector<std::string> list);

void printN(std::string text); // print with no line


#endif