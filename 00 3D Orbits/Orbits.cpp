#include "Orbits.hpp"
/*
Code written by:

Michael Gunnarson
started July 17, 2023

This started as an exercise to learn the basics of OOP C++.
Although still helpful in that regard, this project has
turned into more.  I have decided to make a minimum viable 
product for fully 3D orbits and trajectory design.  Part of
this is a redesign of my senior design code in Professor Joshi's
Space Mission Design course, TA'd by Adam Zufall.  The other
part of this project is the realization that having an easy
trajectory creator in C++ for speed, if integratable into python,
might be interfaceable with NASA's GMAT software.  

The ability to take an intermediary step in GMAT and calculate
what kind of state vector I need to get the trajectory I want is
helpful.  If a lambert solver is used then the trajectories won't
actually be that off.  As far as I understand, GMAT requires you to 
know many parts of an orbit in order to model it.  I am hoping
this software will give me those orbits, or close enough to it
to apply the VF13ad optimizer to get better fuel efficiency.  
We will have to see.

The other purpose of this program is to give to future UC Davis
students in Joshi's class, which is why I created the CelestialBodies
portion of the script.  Adam Zufall had mentioned that 
if I could get it to work it would be a wonderful tool to give
to students in the future.  So if you're reading this, just 
know a UC Davis alum made it through the class and you can
too.  I know it's wishful thinking, though, I doubt anyone 
will read my notes lmao.  I wish you the best.

-Michael Gunnarson
*/

// STATE VECTOR PUBLIC //
void State::print() {
    pos.print("Position:");
    vel.print("Velocity:");
}

void State::print(std::string s) {
    _print(s);
    pos.print("Position:");
    vel.print("Velocity:");
}

// ORBIT 2D GEOMETRY OBJECT PUBLIC //

void OG2D::print() {
    /* print all variables in object */
    std::vector<long double> vars;
    std::vector<std::string> names,units;
    std::string space = " ";

    // fill vectors
    vars = getVars();
    names = getNames();
    units = getUnits();

    // check for size
    if((vars.size()==names.size()) and (names.size()==units.size())) {
        for(int i=0; i<vars.size(); i++){
            _print_(names.at(i)); // no newline
            _print_(vars.at(i));
            _print_(space);
            _print(units.at(i)); // with newline
        }
    }
    _print(); // user print from SanityFunctions
}

void OG2D::print(std::string s) {
    // s is intended to be "objName:"
    /* print all variables in object */
    std::vector<long double> vars;
    std::vector<std::string> names,units;
    std::string space = " ";
    std::string tab = "    ";

    // fill vectors
    vars = getVars();
    names = getNames();
    units = getUnits();

    // check for size
    if((vars.size()==names.size()) and (names.size()==units.size())) {
        _print(s);
        for(int i=0; i<vars.size(); i++){
            _print_(tab);
            _print_(names.at(i)); // no newline
            if(vars.at(i) != 42069) {
                _print_(vars.at(i));
            } else {
                _print_("variable not filled in");
            }
            _print_(space);
            _print(units.at(i)); // with newline
        }
    }
    _print(); // user print from SanityFunctions
}

void OG2D::fill() {
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
        _print("solution did not converge");
    }
}

// ORBIT 2D GEOMETRY OBJECT PROTECTED //

void OG2D::initialize() {
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

bool OG2D::checkState() {
    /* check if getString and cutString yield the same results. 
    If so, then you have completed the algorithm and can safely exit*/
    std::vector<std::string> all = getString();
    std::vector<std::string> cut = cutString(all);
    bool result = (all == cut);
    return result;

}

std::vector<std::string> OG2D::getString() {
    /* get all strings within the object*/
    std::vector<std::string> s(10);
    s.at(0) = E_;
    s.at(1) = e_;
    s.at(2) = p_;
    s.at(3) = a_;
    s.at(4) = H_;
    s.at(5) = rp_;
    s.at(6) = ra_;
    s.at(7) = vp_;
    s.at(8) = va_;
    s.at(9) = mu_;
    return s;
}

std::vector<std::string> OG2D::cutString(std::vector<std::string> pass) {
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

bool OG2D::in(std::vector<std::string> pass,std::string check) {
    /* iterate through vector of strings, finding if it matches any cases */
    bool flag = false;
    int n = pass.size();
    for(int i=0; i<n; i++) {
        if(pass.at(i) == check) {
            flag = true;
            break;
        }
    }
    return flag;
}

std::vector<long double> OG2D::getVars() {
    std::vector<long double> s(10);
    s.at(0) = E;
    s.at(1) = e;
    s.at(2) = p;
    s.at(3) = a;
    s.at(4) = H;
    s.at(5) = rp;
    s.at(6) = ra;
    s.at(7) = vp;
    s.at(8) = va;
    s.at(9) = mu;
    
    return s;
}

std::vector<std::string> OG2D::getNames() {
    std::vector<std::string> s(10);
    s.at(0) = E__;
    s.at(1) = e__;
    s.at(2) = p__;
    s.at(3) = a__;
    s.at(4) = H__;
    s.at(5) = rp__;
    s.at(6) = ra__;
    s.at(7) = vp__;
    s.at(8) = va__;
    s.at(9) = mu__;
    
    return s;
}

std::vector<std::string> OG2D::getUnits() {
    std::vector<std::string> s(10);
    s.at(0) = __E;
    s.at(1) = __e;
    s.at(2) = __p;
    s.at(3) = __a;
    s.at(4) = __H;
    s.at(5) = __rp;
    s.at(6) = __ra;
    s.at(7) = __vp;
    s.at(8) = __va;
    s.at(9) = __mu;
    
    return s;
}


void OG2D::_E() {
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

void OG2D::_e() {
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

void OG2D::_p() {
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

void OG2D::_a() {
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

void OG2D::_H() {
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

void OG2D::_rp() {
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

void OG2D::_ra() {
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

void OG2D::_vp() {
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

void OG2D::_va() {
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



// ORBIT 2D POSITION PUBLIC //

// operating functions
OP2D OP2D::operator=(const OG2D& o2) { // return 
    OP2D o;
    o.E = o2.E;
    o.e = o2.e;
    o.p = o2.p;
    o.a = o2.a;
    o.H = o2.H;
    o.rp = o2.rp;
    o.ra = o2.ra;
    o.vp = o2.vp;
    o.va = o2.va;
    o.mu = o2.mu;
    return o;
}

void OP2D::fill() {
    // update with new vars
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
        _r();
        _v();
        _nu();
        _fi();
        // check to see if all parameters have been updated
        loopFlag = checkState();
        i++;
    }
    if (i==n) {
        _print("solution did not converge");
    }

}

State OP2D::getState() {
    State s;
    s.pos = getPos();
    s.vel = getVel();
    return s;
}

// ORBIT 2D POSITION PROTECTED //

Vector OP2D::getPos() {
    Vector peri,pos;
    long double INC,ARG,LAN;
    Matrix3 cosine;
    peri.x = r*cos(nu);
    peri.y = r*sin(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    INC = 0;
    ARG = 0;
    LAN = 0;
    cosine = C(LAN,ARG,INC);
    pos = cosine.t() * peri;
    return pos;
}

Vector OP2D::getVel() {
    Vector peri,vel;
    long double vr,vt,INC,ARG,LAN;
    Matrix3 cosine;
    vr = v*sin(fi);
    vt = v*cos(fi);
    peri.x = vr*cos(nu) - vt*sin(nu);
    peri.y = vr*sin(nu) + vt*cos(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    INC = 0;
    ARG = 0;
    LAN = 0;
    cosine = C(LAN,ARG,INC);
    vel = cosine.t() * peri;
    return vel;
}

void OP2D::initialize() {
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
    if(r != 42069) {
        r_ = "r";
    }
    if(v != 42069) {
        v_ = "v";
    }
    if(nu != 42069) {
        nu_ = "nu";
    }
    if(fi != 42069) {
        fi_ = "fi";
    }
}

std::vector<std::string> OP2D::getString() {
    /* get all strings within the object*/
    std::vector<std::string> s(14);
    s.at(0) = E_;
    s.at(1) = e_;
    s.at(2) = p_;
    s.at(3) = a_;
    s.at(4) = H_;
    s.at(5) = rp_;
    s.at(6) = ra_;
    s.at(7) = vp_;
    s.at(8) = va_;
    s.at(9) = mu_;
    s.at(10) = r_;
    s.at(11) = v_;
    s.at(12) = nu_;
    s.at(13) = fi_;

    return s;
}

std::vector<long double> OP2D::getVars() {
    std::vector<long double> s(14);
    s.at(0) = E;
    s.at(1) = e;
    s.at(2) = p;
    s.at(3) = a;
    s.at(4) = H;
    s.at(5) = rp;
    s.at(6) = ra;
    s.at(7) = vp;
    s.at(8) = va;
    s.at(9) = mu;
    s.at(10) = r;
    s.at(11) = v;
    s.at(12) = nu;
    s.at(13) = fi;

    return s;
}

std::vector<std::string> OP2D::getNames() {
    std::vector<std::string> s(14);
    s.at(0) = E__;
    s.at(1) = e__;
    s.at(2) = p__;
    s.at(3) = a__;
    s.at(4) = H__;
    s.at(5) = rp__;
    s.at(6) = ra__;
    s.at(7) = vp__;
    s.at(8) = va__;
    s.at(9) = mu__;
    s.at(10) = r__;
    s.at(11) = v__;
    s.at(12) = nu__;
    s.at(13) = fi__;

    return s;
}

std::vector<std::string> OP2D::getUnits() {
    std::vector<std::string> s(14);
    s.at(0) = __E;
    s.at(1) = __e;
    s.at(2) = __p;
    s.at(3) = __a;
    s.at(4) = __H;
    s.at(5) = __rp;
    s.at(6) = __ra;
    s.at(7) = __vp;
    s.at(8) = __va;
    s.at(9) = __mu;
    s.at(10) = __r;
    s.at(11) = __v;
    s.at(12) = __nu;
    s.at(13) = __fi;

    return s;
}

void OP2D::_E() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"E")) {
        {}
    } else if(in(s,"v") and in(s,"mu") and in(s,"r")) {
        E = pow(this->v,2.0)/2.0 - mu/this->r;
        E_ = "E";
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

void OP2D::_H() {
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
    } else if(in(s,"r") and in(s,"v") and in(s,"fi")) {
        H = r*v*cos(fi);
        H_ = "H";
    }
}

void OP2D::_r() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"r")) {
        {}
    } else if(in(s,"p") and in(s,"e") and in(s,"nu")) {
        r = p/(1.0 + e*cos(nu));
        r_ = "r";
    }
}

void OP2D::_v() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"v")) {
        {}
    } else if(in(s,"E") and in(s,"mu") and in(s,"r")) {
        v = pow(2.0*(E + mu/r),0.5);
        v_ = "v";
    }
}

void OP2D::_nu() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"nu")) {
        {}
    } else if(in(s,"r") and in(s,"p") and in(s,"e")) {
        nu = acos(1.0/e * (p/r -1));
        nu_ = "nu";
    }
}

void OP2D::_fi() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"fi")) {
        {}
    } else if(in(s,"H") and in(s,"r") and in(s,"v")) {
        fi = acos(H/(r*v));
        fi_ = "fi";
    }
}

// ORBIT 3D GEOMETRY PUBLIC //
OG3D OG3D::operator=(const OG2D& o2) {
    OG3D o;
    o.E = o2.E;
    o.e = o2.e;
    o.p = o2.p;
    o.a = o2.a;
    o.H = o2.H;
    o.rp = o2.rp;
    o.ra = o2.ra;
    o.vp = o2.vp;
    o.va = o2.va;
    o.mu = o2.mu;
    return o;
}

// ORBIT 3D GEOMETRY PROTECTED //
std::vector<long double> OG3D::getVars() {
    std::vector<long double> s(13);
    s.at(0) = E;
    s.at(1) = e;
    s.at(2) = p;
    s.at(3) = a;
    s.at(4) = H;
    s.at(5) = rp;
    s.at(6) = ra;
    s.at(7) = vp;
    s.at(8) = va;
    s.at(9) = mu;
    s.at(10) = INC;
    s.at(11) = ARG;
    s.at(12) = LAN;

    return s;
}

std::vector<std::string> OG3D::getNames() {
    std::vector<std::string> s(13);
    s.at(0) = E__;
    s.at(1) = e__;
    s.at(2) = p__;
    s.at(3) = a__;
    s.at(4) = H__;
    s.at(5) = rp__;
    s.at(6) = ra__;
    s.at(7) = vp__;
    s.at(8) = va__;
    s.at(9) = mu__;
    s.at(10) = INC__;
    s.at(11) = ARG__;
    s.at(12) = LAN__;

    return s;
}

std::vector<std::string> OG3D::getUnits() {
    std::vector<std::string> s(13);
    s.at(0) = __E;
    s.at(1) = __e;
    s.at(2) = __p;
    s.at(3) = __a;
    s.at(4) = __H;
    s.at(5) = __rp;
    s.at(6) = __ra;
    s.at(7) = __vp;
    s.at(8) = __va;
    s.at(9) = __mu;
    s.at(10) = __INC;
    s.at(11) = __ARG;
    s.at(12) = __LAN;

    return s;
}

// ORBIT 3D POSITION PROTECTED //
Vector OP3D::getPos() {
    Vector peri,pos;
    Matrix3 cosine;
    peri.x = r*cos(nu);
    peri.y = r*sin(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    cosine = C(LAN,ARG,INC);
    pos = cosine.t() * peri;
    return pos;
}

Vector OP3D::getVel() {
    Vector peri,vel;
    long double vr,vt;
    Matrix3 cosine;
    vr = v*sin(fi);
    vt = v*cos(fi);
    peri.x = vr*cos(nu) - vt*sin(nu);
    peri.y = vr*sin(nu) + vt*cos(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    cosine = C(LAN,ARG,INC);
    vel = cosine.t() * peri;
    return vel;
}

std::vector<long double> OP3D::getVars() {
    std::vector<long double> s(17);
    s.at(0) = E;
    s.at(1) = e;
    s.at(2) = p;
    s.at(3) = a;
    s.at(4) = H;
    s.at(5) = rp;
    s.at(6) = ra;
    s.at(7) = vp;
    s.at(8) = va;
    s.at(9) = mu;
    s.at(10) = r;
    s.at(11) = v;
    s.at(12) = nu;
    s.at(13) = fi;
    s.at(14) = INC;
    s.at(15) = ARG;
    s.at(16) = LAN;

    return s;
}

std::vector<std::string> OP3D::getNames() {
    std::vector<std::string> s(17);
    s.at(0) = E__;
    s.at(1) = e__;
    s.at(2) = p__;
    s.at(3) = a__;
    s.at(4) = H__;
    s.at(5) = rp__;
    s.at(6) = ra__;
    s.at(7) = vp__;
    s.at(8) = va__;
    s.at(9) = mu__;
    s.at(10) = r__;
    s.at(11) = v__;
    s.at(12) = nu__;
    s.at(13) = fi__;
    s.at(14) = INC__;
    s.at(15) = ARG__;
    s.at(16) = LAN__;

    return s;
}

std::vector<std::string> OP3D::getUnits() {
    std::vector<std::string> s(17);
    s.at(0) = __E;
    s.at(1) = __e;
    s.at(2) = __p;
    s.at(3) = __a;
    s.at(4) = __H;
    s.at(5) = __rp;
    s.at(6) = __ra;
    s.at(7) = __vp;
    s.at(8) = __va;
    s.at(9) = __mu;
    s.at(10) = __r;
    s.at(11) = __v;
    s.at(12) = __nu;
    s.at(13) = __fi;
    s.at(14) = __INC;
    s.at(15) = __ARG;
    s.at(16) = __LAN;

    return s;
}