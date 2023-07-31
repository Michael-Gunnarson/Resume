#include "CelestialBodies.hpp"

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
    Mercury.orbit.INC = 7.00497902*deg;
    // https://en.wikipedia.org/wiki/Longitude_of_the_periapsis
    Mercury.orbit.ARG = (48.33076593-77.45779628)*deg; // website gives longitude of periapsis, not argument of periapsis
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
    Venus.orbit.INC = 3.39467605*deg;
    Venus.orbit.ARG = (76.67984255-131.60246718)*deg;
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
    Earth.orbit.INC = -0.00001531*deg;
    Earth.orbit.ARG = (0.0-102.93768193)*deg;
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
    Mars.orbit.INC = 1.84969142*deg;
    Mars.orbit.ARG = (49.55953891+23.94362959)*deg;
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
    Jupiter.orbit.INC = 1.30439695*deg;
    Jupiter.orbit.ARG = (100.47390909-14.72847983)*deg;
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
    Saturn.orbit.INC = 2.48599187*deg;
    Saturn.orbit.ARG = (113.66242448-92.59887831)*deg;
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
    Uranus.orbit.INC = 0.77263783*deg;
    Uranus.orbit.ARG = (74.01692503-170.95427630)*deg;
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
    Neptune.orbit.INC = 1.77004347*deg;
    Neptune.orbit.ARG = (131.78422574-44.96476227)*deg;
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
    Ceres.orbit.INC = 10.58634325178498*deg;
    Ceres.orbit.ARG = 73.47046097923743*deg;
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
    Pluto.orbit.INC = 17.089000919562*deg;
    Pluto.orbit.ARG = 112.5971416774872*deg;
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
    Eris.orbit.INC = 43.76049565740999*deg;
    Eris.orbit.ARG = 150.9970340507818*deg;
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
    Makemake.orbit.INC = 29.02544951221016*deg;
    Makemake.orbit.ARG = 296.1558178544629*deg;
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
    Haumea.orbit.INC = 28.20947938353613*deg;
    Haumea.orbit.ARG = 240.7021053657648*deg;
    Haumea.orbit.LAN = 121.9749971146863*deg;
    Haumea.SOI = SOI(Haumea,Sun);
    return Haumea;
}