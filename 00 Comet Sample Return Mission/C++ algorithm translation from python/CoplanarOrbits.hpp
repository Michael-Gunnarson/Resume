#ifndef CoplanarOrbits_hpp
#define CoplanarOrbits_hpp

#include <string>
#include <vector>

class Orbits {
    public:
        // variables
        // set initial value to dummy number, to allow initializer
        // to know which vars were inputted
        long double E = 42069; // specific energy
        long double e = 42069; // eccentricity
        long double p = 42069; // semilatus rectum
        long double a = 42069; // semimajor axis
        long double H = 42069; // angular momentum
        long double rp = 42069; // radius periapsis
        long double ra = 42069; // radius apoapsis
        long double vp = 42069; // velocity periapsis
        long double va = 42069; // velocity apoapsis
        long double mu = 42069; // gravitational parameter of the body you're orbiting
        /*
        // variables for position configuration, to implement later
        float v; // velocity
        float r; // radius
        float nu; // true anomaly
        float phi; // elevation angle
        float delta; // turning angle
        */
        void printAll_(); // print all strings to cmnd line for troubleshooting
        void printAll(); // print all variables to cmnd line for troubleshooting
        void fill(); // loop until params full

    private:
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
        // managerial functions
        bool checkState();
        std::vector<std::string> getString(); // get all strings in class
        std::vector<std::string> cutString(std::vector<std::string> pass); // get all nonempty strings  
        bool in(std::vector<std::string> pass,std::string check); // check if string in vector 
        void initialize(); // check what variables have been entered and add the strings as needed
        
};

#endif