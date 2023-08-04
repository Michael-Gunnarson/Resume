#ifndef Orbits_hpp 
#define Orbits_hpp

// <string> and <vector> included in SanityFunctions.hpp and VectorMath.hpp
#include <tgmath.h>
#include "SanityFunctions.hpp" // print functions
#include "VectorMath.hpp" // state vectors and matrix multiplication

/*
Classes:
State: contains cartesian state vectors, position and velocity
OG2D: Orbit geometry, 2-dimensional perifocal frame
OP2D: Orbit position, contains orbit geometry 2D and positioning
information on it
OG3D: Orbit geometry with 3D specifiers, inclination, argument
 of the periapsis, and longitude of the asccending node
 OP3D: orbit geometry and positioning for a 3D orbit

*/

class State {
    public:
        Vector pos;
        Vector vel;
        long double mu;
        
        void print();
        void print(std::string s);
};

class OG2D {
    public:
        // variables
        // set initial value to dummy number, to allow initializer
        // to know which vars were inputted
        long double E = 42069; // specific energy, m^2/s^2
        long double e = 42069; // eccentricity, unitless
        long double p = 42069; // semilatus rectum, m
        long double a = 42069; // semimajor axis, m
        long double H = 42069; // angular momentum, m^2/s
        long double rp = 42069; // radius periapsis, m
        long double ra = 42069; // radius apoapsis, m
        long double vp = 42069; // velocity periapsis, m/s
        long double va = 42069; // velocity apoapsis, m/s
        long double mu = 42069; // gravitational parameter of the body you're orbiting, m^3/s^2


        void print(); // print all strings to cmnd line for troubleshooting
        void print(std::string s);
        void fill(); // loop until params full

    protected:
        // storage tracking
        // initializes as empty
        std::string E_;
        std::string e_;
        std::string p_;
        std::string a_;
        std::string H_;
        std::string rp_;
        std::string ra_;
        std::string vp_;
        std::string va_;
        std::string mu_;

        // variable names for printing:
        std::string E__ = "Specific Energy ";
        std::string e__ = "Eccentricity ";
        std::string p__ = "Semilatus Rectum ";
        std::string a__ = "Semimajor Axis ";
        std::string H__ = "Angular Momentum ";
        std::string rp__ = "Radius Periapsis ";
        std::string ra__ = "Radius Apoapsis ";
        std::string vp__ = "Velocity Periapsis ";
        std::string va__ = "Velocity Apoapsis ";
        std::string mu__ = "Gravitational Parameter ";

        // units for printing:
        std::string __E = "m^2/s^2";
        std::string __e = "";
        std::string __p = "m";
        std::string __a = "m";
        std::string __H = "m^2/s";
        std::string __rp = "m";
        std::string __ra = "m";
        std::string __vp = "m/s";
        std::string __va = "m/s";
        std::string __mu = "m^3/s^2";

        
        // managerial functions
        void initialize(); // check what variables have been entered and add the strings as needed
        bool checkState(); // check if solution has converged
        std::vector<std::string> getString(); // get all tracking strings
        std::vector<std::string> cutString(std::vector<std::string> pass); // get all nonempty strings from vec of strings
        bool in(std::vector<std::string> pass,std::string check); // check if string in vector 
        virtual std::vector<long double> getVars();
        virtual std::vector<std::string> getNames();
        virtual std::vector<std::string> getUnits();

        // function definitions
        // calculate orbital parameters
        void _E();
        void _e();
        void _p();
        void _a();
        void _H();
        void _rp();
        void _ra();
        void _vp();
        void _va();
        void _mu();
        
};


class OP2D: public OG2D {
    public:
        // operating functions
        OP2D operator=(const OG2D& o2);
        OP2D operator=(const OP2D& o2);

        // add position variables to OG2D class
        long double r = 42069; // radius, m
        long double v = 42069; // velocity m/s
        long double nu = 42069; // true anomaly, radians
        long double fi = 42069; // elevation angle, radians


        void fill(); // need to update with added function calls
        State getState();
        
    protected:
        // add storage tracking vars
        std::string r_;
        std::string v_;
        std::string nu_;
        std::string fi_;

        // add names for printing
        std::string r__ = "Radius ";
        std::string v__ = "Velocity ";
        std::string nu__ = "True Anomaly ";
        std::string fi__ = "Elevation Angle ";

        // add units for printing
        std::string __r = "m";
        std::string __v = "m/s";
        std::string __nu = "rad";
        std::string __fi = "rad";

        // state vectors:
        virtual Vector getPos();
        virtual Vector getVel();

        void initialize(); // update vars check
        std::vector<std::string> getString(); // update which tracking strings to grab
        
        // update all printing retrievals
        std::vector<long double> getVars() override;
        std::vector<std::string> getNames() override;
        std::vector<std::string> getUnits() override;

        // update funciton definitions:
        void _E();
        void _H();

        // add function definitions
        void _r();
        void _v();
        void _nu();
        void _fi();
};

// UPDATE ORBIT ORIENTATION WITH TRACKING STRINGS, INITIALIZERS, ETC
class OG3D: public OG2D {
    public:
        // overload operators as a way of type casting
        OG3D operator=(const OG2D& o2);

        // add orientation variables to OG2D class
        long double i = 42069; // inclination, radians
        long double arg = 42069; // argument of the periapsis, radians
        long double LAN = 42069; // longitude of the ascending node, radians
        // so these need to be populated by the user for the most part
        // if i can get the lambert theorem solver working then it
        // should automatically fill these in for transfer orbits.

    protected:
        // add names for printing
        std::string i__ = "Inclination ";
        std::string ARG__ = "Argument of the Periapsis ";
        std::string LAN__ = "Longitude of the Ascending Node ";

        // add units for printing
        std::string __i = "rad";
        std::string __ARG = "rad";
        std::string __LAN = "rad";


        // update all printing retrievals
        std::vector<long double> getVars() override;
        std::vector<std::string> getNames() override;
        std::vector<std::string> getUnits() override;
};

// UPDATE ORIENTATION PLUS POSITION WITH INITIALIZERS, TRACKING STRINGS
// PRINTING FUNCTIONS, ETC.

class OP3D: public OP2D {
    /* Take the perifocal orientation (OP2D) and the 3D geometry
    (OG3D) to combine into one class */
    public:

        OG3D operator=(const OP3D& o2);

        long double i = 42069; // inclination, radians
        long double arg = 42069; // argument of the periapsis, radians
        long double LAN = 42069; // longitude of the ascending node, radians

    protected:
        // add names for printing
        std::string i__ = "Inclination ";
        std::string ARG__ = "Argument of the Periapsis ";
        std::string LAN__ = "Longitude of the Ascending Node ";

        // add units for printing
        std::string __i = "rad";
        std::string __ARG = "rad";
        std::string __LAN = "rad";

        
        Vector getPos() override;
        Vector getVel() override;

        // update all printing retrievals
        std::vector<long double> getVars() override;
        std::vector<std::string> getNames() override;
        std::vector<std::string> getUnits() override;
};

class OP3DIter: public OP3D {
    public:

        //overload operator as cast
        OP3DIter operator=(const OP3D& o2);

        std::vector<long double> r;
        std::vector<long double> v;
        std::vector<long double> nu;
        std::vector<long double> fi;

        void fill(); // chnage the looping algorithm to handle vectors
        std::vector<State> getState(); // change getState to acquire all states
        void print(); // getVars, getNames, and getUnits should all still be the same

        long int iter = 0; // global variable for vector tracking
        int vecLength;

    protected:

        Vector getPos() override;
        Vector getVel() override;

        void initVectors(); // fill empty vectors with dummy number
        void initialize(); // check for dummy number at current iteration

        // update funciton definitions:
        void _E();
        void _H();

        // add function definitions
        void _r();
        void _v();
        void _nu();
        void _fi();
};

// operator shorthand create iteration object
OP3DIter operator+(const OP3D&, const OP3D&);

// Celestial Bodies

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