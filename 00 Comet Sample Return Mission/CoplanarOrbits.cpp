/* 
Code written by:
Michael Gunnarson

July 17, 2023
This is mostly an exercise to learn C++, I am
porting my main coplanar orbits solving algorithm
from python to C++.  This is a proof of concept
requiring object oriented programming and figuring out
how to flag down variables.

mathematics are derived from Introduction to Spaceflight by Francis J. Hale
 */

#include <vector>
#include <string>
#include <tgmath.h>
#include "CoplanarOrbits.hpp"
#include "SanityFunctions.hpp" // need print functions because i don't like cout

void Orbits::initialize() {
    if(E != 42069) {
        E_ = "E";
    }
    if(e != 42069) {
        e_ = "e";
    }
    if(p != 42069) {
        p_ = "p";
    }
    if(a != 42069) {
        a_ = "a";
    }
    if(H != 42069) {
        H_ = "H";
    }
    if(rp != 42069) {
        rp_ = "rp";
    }
    if(ra != 42069) {
        ra_ = "ra";
    }
    if(vp != 42069) {
        vp_ = "vp";
    }
    if(va != 42069) {
        va_ = "va";
    }
    if(mu != 42069) {
        mu_ = "mu";
    }
}

void Orbits::fill() {
    /* Loop through solving funcitons
    until each parameter is filled. 
    This is accomplished through storing strings
    and updating them as parameters are solved*/
    int i = 0;
    int n = 14*5;
    //int n = 5;
    bool loopFlag = false;
    initialize();
    while(loopFlag==false and i < n) {
        // try to solve all parameters
        _E();
        _e();
        _p();
        _a();
        _H();
        _rp();
        _ra();
        _vp();
        _va();
        // check to see if all parameters have been updated
        loopFlag = checkState();
        i++;
    }
    if (i==n) {
        print("solution did not converge");
    }
}

bool Orbits::checkState() {
    /* check if getString and cutString yield the same results. 
    If so, then you have completed the algorithm and can safely exit*/
    std::vector<std::string> all = getString();
    std::vector<std::string> cut = cutString(all);
    bool result = (all == cut);
    return result;

}

std::vector<std::string> Orbits::getString() {
    /* get all strings within the object*/
    std::vector<std::string> s;//(9);
    s.push_back(this->E_);
    s.push_back(this->e_);
    s.push_back(this->p_);
    s.push_back(this->a_);
    s.push_back(this->H_);
    s.push_back(this->rp_);
    s.push_back(this->ra_);
    s.push_back(this->vp_);
    s.push_back(this->va_);
    s.push_back(this->mu_);
    return s;
}

std::vector<std::string> Orbits::cutString(std::vector<std::string> pass) {
    /* take vector of strings and find all non-empty cases */
    std::vector<std::string> s;
    int n = pass.size();
    for(int i=0; i<n; i++) {
        if(pass.at(i) != "") {
            s.push_back(pass.at(i));
        }
    }
    return s;
}

bool Orbits::in(std::vector<std::string> pass,std::string check) {
    /* iterate through vector of strings, finding if it matches any cases */
    bool flag = false;
    int n = pass.size();
    for(int i=0; i<n; i++) {
        if(pass.at(i) == check) {
            flag = true;
        }
    }
    return flag;
}

void Orbits::printAll_() {
    /* print all strings in object */
    printN("Specific Energy ");
    print(E_);
    printN("Eccentricity ");
    print(e_);
    printN("Semilatus Rectum ");
    print(p_);
    printN("Semimajor Axis ");
    print(a_);
    printN("Angular Momentum ");
    print(H_);
    printN("Radius Periapsis ");
    print(rp_);
    printN("Radius Apoapsis ");
    print(ra_);
    printN("Velocity Periapsis ");
    print(vp_);
    printN("Velocity Apoapsis ");
    print(va_);
    printN("Gravitational Parameter ");
    print(mu_);
    std::string blank = "";
    print(blank);

    
}

void Orbits::printAll() {
    /* print all strings in object */
    printN("Specific Energy ");
    print(E);
    printN("Eccentricity ");
    print(e);
    printN("Semilatus Rectum ");
    print(p);
    printN("Semimajor Axis ");
    print(a);
    printN("Angular Momentum ");
    print(H);
    printN("Radius Periapsis ");
    print(rp);
    printN("Radius Apoapsis ");
    print(ra);
    printN("Velocity Periapsis ");
    print(vp);
    printN("Velocity Apoapsis ");
    print(va);
    printN("Gravitational Parameter ");
    print(mu);
    std::string blank = "";
    print(blank);
    
}

/*
write strings:
E_ will contain "E"
e_ will contain "e"
etc etc
*/

void Orbits::_E() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"E")) {
        {}
    } else if(in(s,"mu") and in(s,"a")) {
        E = -mu/(2.0*a); // calc, etc
        E_ = "E"; // update string flag that specific energy has been calculated
    } else if(in(s,"vp") and in(s,"mu") and in(s,"rp")) {
        E = pow(vp,2.0)/2.0 - mu/rp;
        E_ = "E";
    } else if(in(s,"va") and in(s,"mu") and in(s,"ra")) {
        E = pow(va,2.0)/2.0 - mu/ra;
        E_ = "E";
    }
}

void Orbits::_e() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"e")) {
        {}
    } else if(in(s,"rp") and in(s,"ra")) {
        e = (ra-rp)/(ra+rp);
        e_ = "e";
    } else if(in(s,"E") and in(s,"mu") and in(s,"H")) {
        e = pow(1.0+2.0*E*pow(H,2.0)/pow(mu,2.0),0.5);
        e_ = "e";
    }
}

void Orbits::_p() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"p")) {
        {}
    } else if(in(s,"a") and in(s,"e")) {
        p = a*(1.0-pow(e,2.0));
        p_ = "p";
    } else if(in(s,"H") and in(s,"mu")) {
        p = pow(H,2.0)/mu;
        p_ = "p";
    }
}

void Orbits::_a() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"a")) {
        {}
    } else if(in(s,"ra") and in(s,"rp")) {
        a = (ra+rp)/2.0;
        a_ = "a";
    } else if(in(s,"ra") and in(s,"e")) {
        a = ra/(1.0+e);
        a_ = "a";
    } else if(in(s,"rp") and in(s,"e")) {
        a = rp/(1.0-e);
        a_ = "a";
    } else if(in(s,"mu") and in(s,"E")) {
        a = -mu/(2.0*E);
        a_ = "a";
    }
}

void Orbits::_H() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"H")) {
        {}
    } else if(in(s,"mu") and in(s,"p")) {
        H = pow(mu*p,0.5);
        H_ = "H";
    } else if(in(s,"ra") and in(s,"va")) {
        H = ra*va;
        H_ = "H";
    } else if(in(s,"rp") and in(s,"vp")) {
        H = rp*vp;
        H_ = "H";
    } 
}

void Orbits::_rp() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"rp")) {
        {}
    } else if(in(s,"a") and in(s,"e")) {
        rp = a*(1.0-e);
        rp_ = "rp";
    } else if(in(s,"a") and in(s,"ra")) {
        rp = 2.0*a-ra;
        rp_ = "rp";
    } else if(in(s,"p") and in(s,"e")) {
        rp = p/(1.0+e);
        rp_ = "rp";
    }
}

void Orbits::_ra() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"ra")) {
        {}
    } else if(in(s,"p") and in(s,"e")) {
        ra = p/(1.0-e);
        ra_ = "ra";
    } else if(in(s,"a") and in(s,"e")) {
        ra = a*(1.0+e);
        ra_ = "ra";
    } else if(in(s,"a") and in(s,"rp")) {
        ra = 2.0*a-rp;
        ra_ = "ra";
    }
}

void Orbits::_vp() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"vp")) {
        {}
    } else if(in(s,"H") and in(s,"rp")) {
        vp = H/rp;
        vp_ = "vp";
    } else if(in(s,"va") and in(s,"rp") and in(s,"ra")) {
        vp = ra*va/rp;
        vp_ = "vp";
    }
}

void Orbits::_va() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"va")) {
        {}
    } else if(in(s,"H") and in(s,"ra")) {
        va = H/ra;
        va_ = "va";
    } else if(in(s,"vp") and in(s,"rp") and in(s,"ra")) {
        va = rp*vp/ra;
        va_ = "va";
    }
}