#include <tgmath.h>
#include "VectorMath.hpp"
#include "SanityFunctions.hpp"



// https://www.tutorialspoint.com/cplusplus/cpp_overloading.htm

// overload operators
Vector Vector::operator+(const Vector& v) {
    Vector s;
    s.x = x + v.x;
    s.y = y + v.y;
    s.z = z + v.z;
    return s;
}

Vector Vector::operator+(const long double& c) {
    Vector s;
    s.x = x + c;
    s.y = y + c;
    s.z = z + c;
    return s;
}

Vector Vector::operator-(const Vector& v) {
    Vector s;
    s.x = x - v.x;
    s.y = y - v.y;
    s.z = z - v.z;
    return s;
}

Vector Vector::operator-(const long double& c) {
    Vector s;
    s.x = x - c;
    s.y = y - c;
    s.z = z - c;
    return s;
}

Vector Vector::operator%(const Vector& v) {
    Vector s;
    s.x = y*v.z - z*v.y;
    s.y = z*v.x - x*v.z;
    s.z = x*v.y - y*v.x;
    return s;
}

Vector Vector::operator*(const long double& c) {
    Vector s;
    s.x = c * x;
    s.y = c * y;
    s.z = c * z;
    return s;
}

long double Vector::operator*(const Vector& v) {
    return x*v.x + y*v.y + z*v.z;
}



long double Vector::length() {
    return pow(*this * *this,0.5); // i think this is the single worst line of code i've ever written, i'm so proud
}

void Vector::print() {
    _print(x);
    _print(y);
    _print(z);
    _print();
}

void Vector::print(std::string s) {
    // s is intended to be "varName ="
    std::string tab;
    tab = "    "; // four whitespace
    _print(s);
    _print_(tab);
    _print(x);
    _print_(tab);
    _print(y);
    _print_(tab);
    _print(z);
    _print();
}

long double angle(Vector v1, Vector v2) {
    return acos(cosTheta(v1,v2));
}

long double sinTheta(Vector v1, Vector v2) {
    return (v1 % v2).length()/(v1.length()*v2.length());
}

long double cosTheta(Vector v1, Vector v2) {
    return (v1 * v2)/(v1.length()*v2.length());
}

long double distance(Vector v1, Vector v2) {
    long double r1,r2;
    r1 = v1.length();
    r2 = v2.length();
    return pow(pow(r1,2.0) + pow(r2,2.0) - 2.0*r1*r2*cosTheta(v1,v2),0.5);
}

Matrix3 Matrix3::t() {
    Matrix3 output;
    output.row1 = {row1.at(0),row2.at(0),row3.at(0)};
    output.row2 = {row1.at(1),row2.at(1),row3.at(1)};
    output.row3 = {row1.at(2),row2.at(2),row3.at(2)};
    return output;
}

void Matrix3::print() {
    std::string t;
    t = "  ";
    _print_(row1.at(0));
    _print_(t);
    _print_(row1.at(1));
    _print_(t);
    _print(row1.at(2));

    _print_(row2.at(0));
    _print_(t);
    _print_(row2.at(1));
    _print_(t);
    _print(row2.at(2));

    _print_(row3.at(0));
    _print_(t);
    _print_(row3.at(1));
    _print_(t);
    _print(row3.at(2));
    _print();
}

void Matrix3::print(std::string s) {
    // s is intended to be "varName ="
    std::string t,tab;
    t = "  "; // spaces between numbers
    tab = "    "; // indentation

    _print(s);

    _print_(tab);
    _print_(row1.at(0));
    _print_(t);
    _print_(row1.at(1));
    _print_(t);
    _print(row1.at(2));

    _print_(tab);
    _print_(row2.at(0));
    _print_(t);
    _print_(row2.at(1));
    _print_(t);
    _print(row2.at(2));

    _print_(tab);
    _print_(row3.at(0));
    _print_(t);
    _print_(row3.at(1));
    _print_(t);
    _print(row3.at(2));
    _print();
}

Matrix3 C(long double LAN, long double ARG, long double INC) {
    /*
    LAN : Longitude of ascending node
    ARG : argument of the periapsis
    INC : orbital inclination
    (units in radians)

    returns Matrix3 object, a cosine matrix to convert from
    universal ref (IJK). frame to perifocal ref frame (ijk) 
    {ijk} = C * {IJK}
    {IJK} = C.t() * {ijk}
    */

    long double C11,C12,C13,C21,C22,C23,C31,C32,C33;
    C11 = cos(LAN)*cos(ARG) - sin(LAN)*sin(ARG)*cos(INC);
    C12 = sin(LAN)*cos(ARG) + cos(LAN)*sin(ARG)*cos(INC);
    C13 = sin(ARG)*sin(INC);
    C21 = -cos(LAN)*sin(ARG) - sin(LAN)*cos(ARG)*cos(INC);
    C22 = -sin(LAN)*sin(ARG) + cos(LAN)*cos(ARG)*cos(INC);
    C23 = cos(ARG)*sin(INC);
    C31 = sin(LAN)*sin(INC);
    C32 = -cos(LAN)*sin(INC);
    C33 = cos(INC);

    Matrix3 C;
    C.row1 = {C11,C12,C13};
    C.row2 = {C21,C22,C23};
    C.row3 = {C31,C32,C33};

    return C;
}

Vector Matrix3::operator*(const Vector& v) {
    Vector r; // return
    r.x = row1.at(0)*v.x + row1.at(1)*v.y + row1.at(2)*v.z;
    r.y = row2.at(0)*v.x + row2.at(1)*v.y + row2.at(2)*v.z;
    r.z = row3.at(0)*v.x + row3.at(1)*v.y + row3.at(2)*v.z;
    return r;
}