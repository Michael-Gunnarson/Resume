#ifndef Transfers_hpp
#define Transfers_hpp

#include "Orbits.hpp"

// we want TOF measurements, flyby's, and general transfers

class dVStor {
    public:
        long double dv1;
        long double dv2;
        long double dvTot;
};

class Maneuver {
    public:
        dVStor dv;
        long double TOF;
        OG3D transferOrbit;
        OP3DIter transferOrbitPositions;
};

class TransferRange {
    public:
        Maneuver e;
        Maneuver p;
        Maneuver hdVMin;
        Maneuver hTOFMin;
        Maneuver hTOFdVMin;
};


// convert state vector to OP3D
OP3D stateToOrbit(State s);

TransferRange Transfer(OP3D start, OP3D end, Star star,long double maxdV);
TransferRange Transfer(OP3D start, OP3D end, Body body,long double maxdV);


std::vector<Maneuver> _rangeHyperbola(long double s, long double c, 
Vector r1, Vector r2, Vector HHat,OP3D start, OP3D end, Star star, long double maxdV);
std::vector<Maneuver> _rangeHyperbola(long double s, long double c, 
Vector r1, Vector r2, Vector HHat,OP3D start, OP3D end, Body body, long double maxdV);


std::vector<Maneuver> _extractHyperbola(std::vector<Maneuver> vec);


Maneuver _minEnergyEllipse(long double s, long double c,
long double ct, Vector r1, Vector r2,Vector HHat, OP3D start, 
OP3D end, Star star);
Maneuver _minEnergyEllipse(long double s, long double c,
long double ct, Vector r1, Vector r2,Vector HHat, OP3D start, 
OP3D end, Body body);


Maneuver _minEnergyParabola(long double s, long double c,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Star star);
Maneuver _minEnergyParabola(long double s, long double c,
Vector r1, Vector r2, Vector HHat, OP3D start, OP3D end, Body body);


Maneuver _constrainedHyperbola(long double s, long double c, long double E,
Vector r1,Vector r2, Vector HHat, OP3D start, OP3D end, Star star);
Maneuver _constrainedHyperbola(long double s, long double c, long double E,
Vector r1,Vector r2, Vector HHat, OP3D start, OP3D end, Body body);

dVStor findDV(State start, State end, OP3D pre, OP3D post);
long double dv(State s, OP3D patch);

long double TOF(OP3DIter);
std::vector<long double> _M(OP3DIter);


#endif