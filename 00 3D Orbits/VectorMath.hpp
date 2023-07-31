/*
The purpose of this vector script is limited and not general in scope
Therefore, only a few math functions will be implemented.
Although vectors are nice, I'm only working on 3d vectors.  I'm also
only doing 3x3 matrix * 1x3 vector for coordinate conversions.
*/

#include <vector>
#include <string>

class Vector {
    public:
        // data
        long double x;
        long double y;
        long double z;

        // overload operators
        // ADD, SUBTRACT, MULT, CROSS
        Vector operator+(const Vector& v);
        Vector operator+(const long double& c);

        Vector operator-(const Vector& v);
        Vector operator-(const long double& c);

        long double operator*(const Vector& v); // dot product
        Vector operator*(const long double& c); // multiply by constant
        
        Vector operator%(const Vector& v); // cross product

        // common vector needs
        long double length();
        void print();
        void print(std::string s);
};

// find angle between two vectors
long double angle(Vector v1, Vector v2);
long double sinTheta(Vector v1, Vector v2);
long double cosTheta(Vector v1, Vector v2);
long double distance(Vector v1, Vector v2);

// matrix for coordinate frame transformations mostly
class Matrix3 {
    public:
        std::vector<long double> row1;
        std::vector<long double> row2;
        std::vector<long double> row3;

        Matrix3 t();
        void print();
        void print(std::string s);

        Vector operator*(const Vector& v);

};

// function to generate direction cosine matrix to convert from
// universal ref (IJK). frame to perifocal ref frame (ijk) 
// {ijk} = C * {IJK}
// {IJK} = C.t() * {ijk}
// units in radians
Matrix3 C(long double LAN, long double ARG,long double INC);