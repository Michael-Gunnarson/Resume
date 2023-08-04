#include "Transfers.hpp"

OP3D stateToOrbit(State s) {
    Vector iHat,jHat,kHat; // unit vectors
    iHat = {1,0,0};
    jHat = {0,1,0};
    kHat = {0,0,1};

    Vector r,v,H,e,nHat; // orbital vectors
    r = s.pos;
    v = s.vel;
    H = r%v;
    e = v%H/s.mu - r.hat();
    nHat = kHat%H.hat();

    long double i,arg,LAN,nu;
    i = acos(kHat*H.hat());
    arg = acos(nHat*e.hat());
    LAN = acos(iHat*nHat);
    nu = acos(e.hat()*r.hat());

    OP3D orbit;
    orbit.r = r.length();
    orbit.v = v.length();
    orbit.H = H.length();
    orbit.e = e.length();
    orbit.i = i;
    orbit.arg = arg;
    orbit.LAN = LAN;
    orbit.nu = nu;
    orbit.fill();
    return orbit;
}

TransferRange Transfer(OP3D start, OP3D end, Star star, long double maxdV) {
    TransferRange Transfer;
    State s1,s2;
    s1 = start.getState();
    s2 = end.getState();

    Vector r1,r2;
    r1 = s1.pos;
    r2 = s2.pos; // velocity not continuous across transfers

    // assume direct path from r1 to r2
    // non-direct paths require HHat = -r1%r2/(r1%r2).length();
    Vector HHat;
    HHat = r1%r2/(r1%r2).length();

    // geometric properties:
    long double c,s,ct;
    c = distance(r1,r2);
    s = (r1.length() + r2.length() + c)/2.0;
    ct = cosTheta(r1,r2);

    Transfer.e = _minEnergyEllipse( s, c, ct, r1, r2, HHat, start, end, star);
    Transfer.p = _minEnergyParabola(s,c,r1,r2,HHat,start,end,star);
    std::vector<Maneuver> vec = _rangeHyperbola(s,c,r1,r2,HHat,start,end,star,maxdV);
    std::vector<Maneuver> hRange = _extractHyperbola(vec);
    // unpack
    Transfer.hdVMin = hRange.at(0);
    Transfer.hTOFMin = hRange.at(1);
    Transfer.hTOFdVMin = hRange.at(2);
    return Transfer;
}

TransferRange Transfer(OP3D start, OP3D end, Body body, long double maxdV) {
    TransferRange Transfer;
    State s1,s2;
    s1 = start.getState();
    s2 = end.getState();

    Vector r1,r2;
    r1 = s1.pos;
    r2 = s2.pos; // velocity not continuous across transfers

    // assume direct path from r1 to r2
    // non-direct paths require HHat = -r1%r2/(r1%r2).length();
    Vector HHat;
    HHat = r1%r2/(r1%r2).length();

    // geometric properties:
    long double c,s,ct;
    c = distance(r1,r2);
    s = (r1.length() + r2.length() + c)/2.0;
    ct = cosTheta(r1,r2);

    Transfer.e = _minEnergyEllipse( s, c, ct, r1, r2, HHat, start, end, body);
    Transfer.p = _minEnergyParabola(s,c,r1,r2,HHat,start,end,body);
    std::vector<Maneuver> vec = _rangeHyperbola(s,c,r1,r2,HHat,start,end,body,maxdV);
    std::vector<Maneuver> hRange = _extractHyperbola(vec);
    // unpack
    Transfer.hdVMin = hRange.at(0);
    Transfer.hTOFMin = hRange.at(1);
    Transfer.hTOFdVMin = hRange.at(2);
    return Transfer;
}



std::vector<Maneuver> _rangeHyperbola(long double s, long double c,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Star star,long double maxdV) {
    long double dv = 0.0;
    long double E = 1.0;
    long int i = 0;
    std::vector<Maneuver> vec;
    while(dv < maxdV and i<10000) {
        Maneuver h = _constrainedHyperbola(s,c,E,r1,r2,HHat,start,end,star);
        dv = h.dv.dvTot;
        vec.at(i) = h;
        i++;
    }
    return vec;
}


std::vector<Maneuver> _rangeHyperbola(long double s, long double c,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Body body,long double maxdV) {
    long double dv = 0.0;
    long double E = 1.0;
    long int i = 0;
    std::vector<Maneuver> vec;
    while(dv < maxdV and i<10000) {
        Maneuver h = _constrainedHyperbola(s,c,E,r1,r2,HHat,start,end,body);
        dv = h.dv.dvTot;
        vec.at(i) = h;
        i++;
    }
    return vec;
}



std::vector<Maneuver> _extractHyperbola(std::vector<Maneuver> vec) {
    // find min delta V:
    // find min TOF:
    // find min TOF/dv:
    long double dv,TOF,TOFdV,dvMin,TOFMin,TOFdVMin;
    long int dvInt,TOFInt,tdvInt;

    //initialize
    dvMin = vec.at(0).dv.dvTot;
    dvInt = 0;
    TOFMin = vec.at(0).TOF;
    TOFInt = 0;
    TOFdVMin = TOFMin/dvMin;
    tdvInt = 0;
    
    // find min
    for(long int i=1; i<vec.size(); i++) {
        dv = vec.at(i).dv.dvTot;
        TOF = vec.at(i).TOF;
        TOFdV = TOF/dv;
        if(dv < dvMin) {
            dvMin = dv;
            dvInt = i;
        }
        if(TOF < TOFMin) {
            TOFMin = TOF;
            TOFInt = i;
        }
        if(TOFdV < TOFdVMin) {
            TOFdVMin = TOFdV;
            tdvInt = i;
        }
    }
    std::vector<Maneuver> hRange(3);
    hRange.at(0) = vec.at(dvInt);
    hRange.at(1) = vec.at(TOFInt);
    hRange.at(2) = vec.at(tdvInt);
    return hRange;
}



Maneuver _minEnergyEllipse(long double s, long double c,
long double ct, Vector r1, Vector r2,Vector HHat,OP3D start, 
OP3D end,Star star) {
    OP2D e,ePos2;
    e.a = s/2.0;
    e.p = r1.length()*r2.length()/c*(1.0-ct);
    e.mu = star.mu;
    ePos2 = e;
    e.r = r1.length();
    ePos2.r = r2.length();
    e.fill();
    ePos2.fill();
    
    Vector H = HHat*e.H;
    Vector v1 = r1*H.length()*tan(e.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(ePos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = star.mu;
    OP3D eMin1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = star.mu;
    OP3D eMin2 = stateToOrbit(S2);

    OP3DIter eMinPos = eMin1 + eMin2;
    
    Maneuver _e;
    _e.transferOrbit = eMin1; // implicit cast from OP3D to OG3D
    _e.transferOrbitPositions = eMinPos;
    _e.dv = findDV(S1,S2,start,end);
    _e.TOF = TOF(eMinPos);
    return _e;
}


Maneuver _minEnergyEllipse(long double s, long double c,
long double ct, Vector r1, Vector r2,Vector HHat,OP3D start, 
OP3D end,Body body) {
    OP2D e,ePos2;
    e.a = s/2.0;
    e.p = r1.length()*r2.length()/c*(1.0-ct);
    e.mu = body.mu;
    ePos2 = e;
    e.r = r1.length();
    ePos2.r = r2.length();
    e.fill();
    ePos2.fill();
    
    Vector H = HHat*e.H;
    Vector v1 = r1*H.length()*tan(e.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(ePos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = body.mu;
    OP3D eMin1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = body.mu;
    OP3D eMin2 = stateToOrbit(S2);

    OP3DIter eMinPos = eMin1 + eMin2;
    
    Maneuver _e;
    _e.transferOrbit = eMin1; // implicit cast from OP3D to OG3D
    _e.transferOrbitPositions = eMinPos;
    _e.dv = findDV(S1,S2,start,end);
    _e.TOF = TOF(eMinPos);
    return _e;
}

Maneuver _minEnergyParabola(long double s,long double c, Vector r1, 
Vector r2, Vector HHat, OP3D start, OP3D end, Star star) {
    OP2D p, pPos2;
    p.mu = star.mu;
    p.E = 0;
    p.p = 4.0*(s-r1.length()) * (s-r2.length())/pow(c,2.0) * pow((pow(s/2,0.5) + pow((s-c)/2,0.5)),2.0);
    // determination of round trip planetary reconnaissance trajectories
    pPos2 = p;
    p.r = r1.length();
    pPos2.r = r2.length();
    p.fill();
    pPos2.fill();

    Vector H = HHat*p.H;
    Vector v1 = r1*H.length()*tan(p.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(pPos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = star.mu;
    OP3D pMin1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = star.mu;
    OP3D pMin2 = stateToOrbit(S2);

    OP3DIter pMinPos = pMin1 + pMin2;
    
    Maneuver _p;
    _p.transferOrbit = pMin1; // implicit cast from OP3D to OG3D
    _p.transferOrbitPositions = pMinPos;
    _p.dv = findDV(S1,S2,start,end);
    _p.TOF = TOF(pMinPos);
    return _p;
}

Maneuver _minEnergyParabola(long double s,long double c, Vector r1, 
Vector r2, Vector HHat, OP3D start, OP3D end, Body body) {
    OP2D p, pPos2;
    p.mu = body.mu;
    p.E = 0;
    p.p = 4.0*(s-r1.length()) * (s-r2.length())/pow(c,2.0) * pow((pow(s/2,0.5) + pow((s-c)/2,0.5)),2.0);
    // determination of round trip planetary reconnaissance trajectories
    pPos2 = p;
    p.r = r1.length();
    pPos2.r = r2.length();
    p.fill();
    pPos2.fill();

    Vector H = HHat*p.H;
    Vector v1 = r1*H.length()*tan(p.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(pPos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = body.mu;
    OP3D pMin1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = body.mu;
    OP3D pMin2 = stateToOrbit(S2);

    OP3DIter pMinPos = pMin1 + pMin2;
    
    Maneuver _p;
    _p.transferOrbit = pMin1; // implicit cast from OP3D to OG3D
    _p.transferOrbitPositions = pMinPos;
    _p.dv = findDV(S1,S2,start,end);
    _p.TOF = TOF(pMinPos);
    return _p;
}

Maneuver _constrainedHyperbola(long double s, long double c, long double E,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Star star) {
    OP2D h,hPos2;
    h.mu = star.mu;
    h.E = E;
    h.a = -h.mu/(2.0*h.E);

    long double alpha,beta;
    alpha = 2*asinh(pow(s/(2.0*h.a),0.5));
    beta = 2*asinh(pow((s-c)/(2.0*h.a),0.5));

    h.p = (4.0*h.a*(s-r1.length())*(s-r2.length())/pow(c,2.0))*pow(sinh((alpha + beta)/2.0),2.0);
    
    hPos2 = h;
    h.r = r1.length();
    hPos2.r = r2.length();
    
    h.fill();
    hPos2.fill();

    Vector H = HHat*h.H;
    Vector v1 = r1*H.length()*tan(h.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(hPos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = star.mu;
    OP3D h1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = star.mu;
    OP3D h2 = stateToOrbit(S2);

    OP3DIter hPos = h1 + h2;
    
    Maneuver _h;
    _h.transferOrbit = h1; // implicit cast from OP3D to OG3D
    _h.transferOrbitPositions = hPos;
    _h.dv = findDV(S1,S2,start,end);
    _h.TOF = TOF(hPos);
    return _h;
}

Maneuver _constrainedHyperbola(long double s, long double c, long double E,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Body body) {
    OP2D h,hPos2;
    h.mu = body.mu;
    h.E = E;
    h.a = -h.mu/(2.0*h.E);

    long double alpha,beta;
    alpha = 2*asinh(pow(s/(2.0*h.a),0.5));
    beta = 2*asinh(pow((s-c)/(2.0*h.a),0.5));

    h.p = (4.0*h.a*(s-r1.length())*(s-r2.length())/pow(c,2.0))*pow(sinh((alpha + beta)/2.0),2.0);
    
    hPos2 = h;
    h.r = r1.length();
    hPos2.r = r2.length();
    
    h.fill();
    hPos2.fill();

    Vector H = HHat*h.H;
    Vector v1 = r1*H.length()*tan(h.fi)/(r1*r1) - r1%H/(r1*r1);
    Vector v2 = r2*H.length()*tan(hPos2.fi)/(r2*r2) - r2%H/(r2*r2);

    // transfer orbit
    State S1;
    S1.pos = r1;
    S1.vel = v1;
    S1.mu = body.mu;
    OP3D h1 = stateToOrbit(S1); // makes sure that nu is in correct plane

    State S2;
    S2.pos = r2;
    S2.pos = v2;
    S2.mu = body.mu;
    OP3D h2 = stateToOrbit(S2);

    OP3DIter hPos = h1 + h2;
    
    Maneuver _h;
    _h.transferOrbit = h1; // implicit cast from OP3D to OG3D
    _h.transferOrbitPositions = hPos;
    _h.dv = findDV(S1,S2,start,end);
    _h.TOF = TOF(hPos);
    return _h;
}

dVStor findDV(State start, State end, OP3D pre, OP3D post) {
    dVStor v;
    v.dv1 = dv(start, pre);
    v.dv2 = dv(end,post);
    v.dvTot = v.dv1 + v.dv2;
    return v;
}

long double dv(State s, OP3D patch) {
    State s2 = patch.getState();
    long double dv;
    dv = abs(s.vel.length() - s2.vel.length());
    return dv;
}

long double TOF(OP3DIter o) {
    long double TOF;
    if(o.e < 1.0 and o.e <= 0.0) {
        // elliptic/circular case
        std::vector<long double> M = _M(o);
        TOF = pow(pow(o.a,3.0)/o.mu,0.5) * (M.at(o.vecLength)-M.at(0));
    } else if(o.E == 0.0) {
        // parabolic case
        std::vector<long double> M = _M(o);
        TOF = pow(o.H,3.0)/(2.0*pow(o.mu,2.0)) * (tan(o.nu.at(o.vecLength)/2.0) + 1.0/3.0*pow(tan(o.nu.at(o.vecLength)/2.0),3.0) - tan(o.nu.at(o.vecLength)/2.0) - 1.0/3.0*pow(tan(o.nu.at(o.vecLength)/2.0),3.0));
    } else if(o.E > 0.0) {
        // hyperbolic case
        std::vector<long double> M = _M(o);
        TOF = pow(pow(-o.a,3.0)/o.mu,0.5) * (M.at(o.vecLength) - M.at(0));
    }
    return TOF;
}

std::vector<long double> _M(OP3DIter o) {
    std::vector<long double> s(2);
    if(o.e < 1.0 and o.e >= 0.0) {
        // elliptic + circular case
        for(int i=0; i<o.vecLength; i++) {
            long double u = acos((o.e + cos(o.nu.at(i)))/(1.0 + o.e*cos(o.nu.at(i))));
            long double M = u - o.e*sin(u);
            s.at(i) = M;
        }
    } else if(o.e>0) {
        // hyperbolic case
        for(int i=0; i<o.vecLength; i++) {
            long double F = acosh((o.e + cos(o.nu.at(i)))/(1.0 + o.e*cos(o.nu.at(i))));
            long double M = o.e*sinh(F) - F;
            s.at(i) = M;
        }
    }
    return s;
}