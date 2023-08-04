#include "VectorMath.hpp"
#include "SanityFunctions.hpp"

// g++ VectorMath.cpp SanityFunctions.cpp VectorMathTest.cpp
// .\a.exe

int main() {
    /*CartesianVector x = {1,2,3};
    CartesianVector y = {4,5,6};
    CartesianVector z,cross;
    long double dotProduct,l,l2;
    z = x + y;
    _print(z.x);
    _print(z.y);
    _print(z.z);
    z = x - y;
    _print(z.x);
    _print(z.y);
    _print(z.z);
    dotProduct = x * y;
    _print(dotProduct);
    l = x.length();
    l2 = y.length();
    _print(l);
    _print(l2);
    cross = x%y;
    _print(cross.x);
    _print(cross.y);
    _print(cross.z);
*/
    Matrix3 A;
    A = {{1,2,3},{4,5,6},{7,8,9}};
    A.print("A=");

    Vector b = {10,11,12};
    b.print("b =");

    Vector x = A*b; 
    x.print("x =");

    Matrix3 Cosine = C(0.1,0.5,1);
    Cosine.print("C =");

    Matrix3 InvCosine = Cosine.t();
    InvCosine.print("C' =");

    return 0;
}