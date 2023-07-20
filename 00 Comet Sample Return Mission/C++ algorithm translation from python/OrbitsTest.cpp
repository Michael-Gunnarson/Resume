/*
Code written by:
Michael Gunnarson

July 17, 2023
Test of algorithm to find coplanar orbital elements using
two body assumtion and smallest number of parameters that
define an orbit.

Preliminary implementation does not include
the position orientation of the object
like python does.

*/

#include "CoplanarOrbits.hpp"
#include "SanityFunctions.hpp"

int main(){

    // basic orbit class text
    Orbits earth, GEO;

    long double oneAU = 1.495978707e11; // meters
    long double aEarth = 1.000*oneAU; // m
    long double eEarth = 0.0167;
    long double rEarth = 6.378e6;

    long double muSun = 1.327e20; // m^3/s^2
    long double muEarth = 3.986e14; // m
    
    // earth around sun 
    earth.e = eEarth;
    earth.a = aEarth;
    earth.mu = muSun;
    earth.fill();
    
    earth.printAll();
    earth.printAll_();

    GEO.e = 0.0;
    GEO.rp = 35786.0*1000.0 + rEarth; // https://en.wikipedia.org/wiki/Geosynchronous_orbit 
    GEO.mu = muEarth;
    GEO.fill();

    GEO.printAll();
    GEO.printAll_();
    

    //return
    return 0;
}