#ifndef Celestialbodies_hpp
#define Celestialbodies_hpp

#include "Orbits.hpp"


class Body {
    public:
        long double mu; // gravitational parameter
        long double SOI; // sphere of influence
        long double r; // mean radius
        OG3D orbit; 

};

class Star {
    public:
        long double mu;
        long double r;
};

long double SOI(Body planet,Body Star); // SOI wrt sun
Star sun();

Body mercury();
Body venus();
Body earth();
Body mars();
Body jupiter();
Body saturn();
Body uranus();
Body neptune();

Body ceres();
Body pluto();
Body eris();
Body makemake();
Body haumea();

#endif