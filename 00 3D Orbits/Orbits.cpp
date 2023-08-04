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

OP2D OP2D::operator=(const OP2D& o2) { // return 
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
    o.r = o2.r;
    o.v = o2.v;
    o.nu = o2.nu;
    o.fi = o2.fi;
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
    s.mu = mu;
    return s;
}

// ORBIT 2D POSITION PROTECTED //

Vector OP2D::getPos() {
    Vector peri,pos;
    long double i,arg,LAN;
    Matrix3 cosine;
    peri.x = r*cos(nu);
    peri.y = r*sin(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    i = 0;
    arg = 0;
    LAN = 0;
    cosine = C(LAN,arg,i);
    pos = cosine.t() * peri;
    return pos;
}

Vector OP2D::getVel() {
    Vector peri,vel;
    long double vr,vt,i,arg,LAN;
    Matrix3 cosine;
    vr = v*sin(fi);
    vt = v*cos(fi);
    peri.x = vr*cos(nu) - vt*sin(nu);
    peri.y = vr*sin(nu) + vt*cos(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    i = 0;
    arg = 0;
    LAN = 0;
    cosine = C(LAN,arg,i);
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
    s.at(10) = i;
    s.at(11) = arg;
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
    s.at(10) = i__;
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
    s.at(10) = __i;
    s.at(11) = __ARG;
    s.at(12) = __LAN;

    return s;
}

// ORBIT 3D POSITION PUBLIC //
OG3D OP3D::operator=(const OP3D& o2) {
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

// ORBIT 3D POSITION PROTECTED //
Vector OP3D::getPos() {
    Vector peri,pos;
    Matrix3 cosine;
    peri.x = r*cos(nu);
    peri.y = r*sin(nu);
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    cosine = C(LAN,arg,i);
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
    cosine = C(LAN,arg,i);
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
    s.at(14) = i;
    s.at(15) = arg;
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
    s.at(14) = i__;
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
    s.at(14) = __i;
    s.at(15) = __ARG;
    s.at(16) = __LAN;

    return s;
}

OP3DIter OP3DIter::operator=(const OP3D& o2) {
    OP3DIter o;
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
    o.r.at(iter) = o2.r;
    o.v.at(iter) = o2.v;
    o.nu.at(iter) = o2.nu;
    o.fi.at(iter) = o2.fi;
    o.iter = 1;
    return o;
}


void OP3DIter::fill() {
    initVectors();
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
    iter++;
    for(int i=1; i<r.size();i++) {
        bool loopFlag = false;
        // RESET TRACKING STRINGS
        r_ = "";
        v_ = "";
        nu_ = "";
        fi_ = "";
        initialize();
        while(loopFlag == false and i < n) {
            _r();
            _v();
            _nu();
            _fi();
            loopFlag = checkState();
            if(loopFlag) {
                iter++;
            }
        }
    }

}

std::vector<State> OP3DIter::getState() {
    std::vector<State> s(vecLength);
    for(int i=0; i<vecLength; i++) {
        State ss;
        ss.mu = mu;
        ss.pos = getPos();
        ss.vel = getVel();
        s.at(i) = ss;
    }
    return s;
}

/*
void OP3DIter::print() {

}*/

Vector OP3DIter::getPos() {
    Vector peri,pos;
    Matrix3 cosine;
    peri.x = r.at(iter)*cos(nu.at(iter));
    peri.y = r.at(iter)*sin(nu.at(iter));
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    cosine = C(LAN,arg,i);
    pos = cosine.t() * peri;
    return pos;
}

Vector OP3DIter::getVel() {
    Vector peri,vel;
    long double vr,vt;
    Matrix3 cosine;
    vr = v.at(iter)*sin(fi.at(iter));
    vt = v.at(iter)*cos(fi.at(iter));
    peri.x = vr*cos(nu.at(iter)) - vt*sin(nu.at(iter));
    peri.y = vr*sin(nu.at(iter)) + vt*cos(nu.at(iter));
    peri.z = 0;
    // CONVERT TO HELIOCENTRIC FRAME
    cosine = C(LAN,arg,i);
    vel = cosine.t() * peri;
    return vel; 
}

void OP3DIter::initVectors() {
    std::vector<std::string> s(3);
    if(r.size()>0) {
        vecLength = r.size();
        s = {"v","fi","nu"};
    } else if(v.size()>0) {
        vecLength = v.size();
        s = {"r","fi","nu"};
    } else if(nu.size()>0) {
        vecLength = nu.size();
        s = {"r","v","fi",};
    } else if(fi.size()>0) {
        vecLength = fi.size();
        s = {"r","v","nu"};
    } else {
        _print("Error: please create vector of evaluation points");
    }

    if(in(s,"r")) {
        for(int i=0; i<vecLength; i++) {
            r.at(i) = 42069;
        }
    }
    if(in(s,"v")) {
        for(int i=0; i<vecLength; i++) {
            v.at(i) = 42069;
        }
    }
    if(in(s,"nu")) {
        for(int i=0; i<vecLength; i++) {
            nu.at(i) = 42069;
        }
    }
    if(in(s,"fi")) {
        for(int i=0; i<vecLength; i++) {
            fi.at(i) = 42069;
        }
    }
}

void OP3DIter::initialize() {
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
    if(r.at(iter) != 42069) {
        r_ = "r";
    }
    if(v.at(iter) != 42069) {
        v_ = "v";
    }
    if(nu.at(iter) != 42069) {
        nu_ = "nu";
    }
    if(fi.at(iter) != 42069) {
        fi_ = "fi";
    }
}

void OP3DIter::_E() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"E")) {
        {}
    } else if(in(s,"v") and in(s,"mu") and in(s,"r")) {
        E = pow(v.at(iter),2.0)/2.0 - mu/r.at(iter);
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

void OP3DIter::_H() {
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
        H = r.at(iter)*v.at(iter)*cos(fi.at(iter));
        H_ = "H";
    }
}

void OP3DIter::_r() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"r")) {
        {}
    } else if(in(s,"p") and in(s,"e") and in(s,"nu")) {
        r.at(iter) = p/(1.0 + e*cos(nu.at(iter)));
        r_ = "r";
    }
}

void OP3DIter::_v() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"v")) {
        {}
    } else if(in(s,"E") and in(s,"mu") and in(s,"r")) {
        v.at(iter) = pow(2.0*(E + mu/r.at(iter)),0.5);
        v_ = "v";
    }
}

void OP3DIter::_nu() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"nu")) {
        {}
    } else if(in(s,"r") and in(s,"p") and in(s,"e")) {
        nu.at(iter) = acos(1.0/e * (p/r.at(iter) -1));
        nu_ = "nu";
    }
}

void OP3DIter::_fi() {
    std::vector<std::string> s;
    s = getString();
    if(in(s,"fi")) {
        {}
    } else if(in(s,"H") and in(s,"r") and in(s,"v")) {
        fi.at(iter) = acos(H/(r.at(iter)*v.at(iter)));
        fi_ = "fi";
    }
}

OP3DIter operator+(const OP3D& o1, const OP3D& o2) {
    if(o1.E != 42069 and (o1.E == o2.E)) {
        OP3DIter o;
        o = o1; // cast, fills first pos. iter = 1
        o.r.at(o.iter) = o2.r;
        o.v.at(o.iter) = o2.v;
        o.nu.at(o.iter) = o2.nu;
        o.fi.at(o.iter) = o2.fi;
        return o;
    } else {
        OP3DIter o;
        return o;
    }
}



// Celestial Bodies

// https://en.wikipedia.org/wiki/Gravitational_constant
#define G 6.67430e-11

// https://en.wikipedia.org/wiki/Astronomical_unit
#define AU 1.495978707e11 // meters

#define deg PI/180

long double SOI(Body planet, Star star) {
    //https://en.wikipedia.org/wiki/Sphere_of_influence_(astrodynamics)
    long double soi = pow(planet.orbit.a * (planet.mu/star.mu),2.0/5.0);
    return soi;
}

// https://en.wikipedia.org/wiki/Sun
// accessed 7/28/2023
Star sun() {
    // assume sun does not revolve around barycenter
    Star Sun;
    Sun.mu = 1.9885e30*G;
    Sun.r = 696342.0*1000.0;
    return Sun;
}

// https://ssd.jpl.nasa.gov/planets/phys_par.html
// https://ssd.jpl.nasa.gov/planets/approx_pos.html
// accessed 7/28/2023
Body mercury() {
    Star Sun = sun();
    Body Mercury;
    Mercury.mu = 0.330103e24 * G;
    Mercury.r = 2439.4*1000;
    Mercury.orbit.a = 0.38709927*AU;
    Mercury.orbit.e = 0.20563593;
    Mercury.orbit.mu = Sun.mu;
    Mercury.orbit.fill();
    Mercury.orbit.i = 7.00497902*deg;
    // https://en.wikipedia.org/wiki/Longitude_of_the_periapsis
    Mercury.orbit.arg = (48.33076593-77.45779628)*deg; // website gives longitude of periapsis, not argument of periapsis
    Mercury.orbit.LAN = 48.33076593*deg;
    Mercury.SOI = SOI(Mercury,Sun);
    return Mercury;
}

Body venus() {
    Star Sun = sun();
    Body Venus;
    Venus.mu = 4.86731e24*G;
    Venus.r = 6051.8*1000;
    Venus.orbit.a = 0.72333566*AU;
    Venus.orbit.e = 0.00677672;
    Venus.orbit.mu = Sun.mu;
    Venus.orbit.fill();
    Venus.orbit.i = 3.39467605*deg;
    Venus.orbit.arg = (76.67984255-131.60246718)*deg;
    Venus.orbit.LAN = 76.67984255*deg;
    Venus.SOI = SOI(Venus,Sun);
    return Venus;
}

Body earth() {
    Star Sun = sun();
    Body Earth;
    Earth.mu = 5.97217e24*G;
    Earth.r = 6371.0084e3;
    Earth.orbit.a = 1.00000261*AU; // values for the Earth/Moon Barycenter
    Earth.orbit.e = 0.01671123;
    Earth.orbit.mu = Sun.mu;
    Earth.orbit.fill();
    Earth.orbit.i = -0.00001531*deg;
    Earth.orbit.arg = (0.0-102.93768193)*deg;
    Earth.orbit.LAN = 0.0*deg;
    Earth.SOI = SOI(Earth,Sun);
    return Earth;
}

Body mars() {
    Star Sun = sun();
    Body Mars;
    Mars.mu = 0.641691e24*G;
    Mars.r = 3389.50e3;
    Mars.orbit.a = 1.52371034*AU;
    Mars.orbit.e = 0.09339410;
    Mars.orbit.mu = Sun.mu;
    Mars.orbit.fill();
    Mars.orbit.i = 1.84969142*deg;
    Mars.orbit.arg = (49.55953891+23.94362959)*deg;
    Mars.orbit.LAN = 49.55953891*deg;
    Mars.SOI = SOI(Mars,Sun);
    return Mars;
}

Body jupiter() {
    Star Sun = sun();
    Body Jupiter;
    Jupiter.mu = 1898.125e24*G;
    Jupiter.r = 69911e3;
    Jupiter.orbit.a = 5.20288700*AU;
    Jupiter.orbit.e = 0.04838624;
    Jupiter.orbit.mu = Sun.mu;
    Jupiter.orbit.fill();
    Jupiter.orbit.i = 1.30439695*deg;
    Jupiter.orbit.arg = (100.47390909-14.72847983)*deg;
    Jupiter.orbit.LAN = 100.47390909*deg;
    Jupiter.SOI = SOI(Jupiter,Sun);
    return Jupiter;
}

Body saturn() {
    Star Sun = sun();
    Body Saturn;
    Saturn.mu = 568.317e24*G;
    Saturn.r = 58232e3;
    Saturn.orbit.a = 9.53667594*AU;
    Saturn.orbit.e = 0.05386179;
    Saturn.orbit.mu = Sun.mu;
    Saturn.orbit.fill();
    Saturn.orbit.i = 2.48599187*deg;
    Saturn.orbit.arg = (113.66242448-92.59887831)*deg;
    Saturn.orbit.LAN = 113.66242448*deg;
    Saturn.SOI = SOI(Saturn,Sun);
    return Saturn;
}

Body uranus() {
    Star Sun = sun();
    Body Uranus;
    Uranus.mu = 86.8099e24*G;
    Uranus.r = 25362e3;
    Uranus.orbit.a = 19.18916464*AU;
    Uranus.orbit.e = 0.04725744;
    Uranus.orbit.mu = Sun.mu;
    Uranus.orbit.fill();
    Uranus.orbit.i = 0.77263783*deg;
    Uranus.orbit.arg = (74.01692503-170.95427630)*deg;
    Uranus.orbit.LAN = 74.01692503*deg;
    Uranus.SOI = SOI(Uranus,Sun);
    return Uranus;
}

Body neptune () {
    Star Sun = sun();
    Body Neptune;
    Neptune.mu = 102.4092e24*G;
    Neptune.r = 24622e3;
    Neptune.orbit.a = 30.06992276*AU;
    Neptune.orbit.e = 0.00859048;
    Neptune.orbit.mu = Sun.mu;
    Neptune.orbit.fill();
    Neptune.orbit.i = 1.77004347*deg;
    Neptune.orbit.arg = (131.78422574-44.96476227)*deg;
    Neptune.orbit.LAN = 131.78422574*deg;
    Neptune.SOI = SOI(Neptune,Sun);
    return Neptune;
}

// https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/
Body ceres() {
    Star Sun = sun();
    Body Ceres;
    Ceres.mu = 938.416e18*G;
    Ceres.r = 469.7e3;
    Ceres.orbit.a = 2.767181746561148*AU;
    Ceres.orbit.e = 0.07881745099351356;
    Ceres.orbit.mu = Sun.mu;
    Ceres.orbit.fill();
    Ceres.orbit.i = 10.58634325178498*deg;
    Ceres.orbit.arg = 73.47046097923743*deg;
    Ceres.orbit.LAN = 80.26014850244651*deg;
    Ceres.SOI = SOI(Ceres,Sun);
    return Ceres;
}

Body pluto() {
    Star Sun = sun();
    Body Pluto;
    Pluto.mu = 13029e18*G;
    Pluto.r = 1188.3e3;
    Pluto.orbit.a = 39.4450697257358*AU;
    Pluto.orbit.e = 0.250248713478499;
    Pluto.orbit.mu = Sun.mu;
    Pluto.orbit.fill();
    Pluto.orbit.i = 17.089000919562*deg;
    Pluto.orbit.arg = 112.5971416774872*deg;
    Pluto.orbit.LAN = 110.3769579554089*deg;
    Pluto.SOI = SOI(Pluto,Sun);
    return Pluto;
}

Body eris() {
    Star Sun = sun();
    Body Eris;
    Eris.mu = 16600e18*G;
    Eris.r = 1200e3;
    Eris.orbit.a = 68.14536545526268*AU;
    Eris.orbit.e = 0.4319581224809352;
    Eris.orbit.mu = Sun.mu;
    Eris.orbit.fill();
    Eris.orbit.i = 43.76049565740999*deg;
    Eris.orbit.arg = 150.9970340507818*deg;
    Eris.orbit.LAN = 36.07068767579804*deg;
    Eris.SOI = SOI(Eris,Sun);
    return Eris;
}

Body makemake() {
    Star Sun = sun();
    Body Makemake;
    Makemake.mu = 3100e18*G;
    Makemake.r = 714e3;
    Makemake.orbit.a = 45.28559041889222*AU;
    Makemake.orbit.e = 0.1657693944892303;
    Makemake.orbit.mu = Sun.mu;
    Makemake.orbit.fill();
    Makemake.orbit.i = 29.02544951221016*deg;
    Makemake.orbit.arg = 296.1558178544629*deg;
    Makemake.orbit.LAN = 79.31408143490017*deg;
    Makemake.SOI = SOI(Makemake,Sun);
    return Makemake;
}

Body haumea() {
    Star Sun = sun();
    Body Haumea;
    Haumea.mu = 4006e18*G;
    Haumea.r = 715e3;
    Haumea.orbit.a = 42.89195837479598*AU;
    Haumea.orbit.e = 0.199923092706661;
    Haumea.orbit.mu = Sun.mu;
    Haumea.orbit.fill();
    Haumea.orbit.i = 28.20947938353613*deg;
    Haumea.orbit.arg = 240.7021053657648*deg;
    Haumea.orbit.LAN = 121.9749971146863*deg;
    Haumea.SOI = SOI(Haumea,Sun);
    return Haumea;
}