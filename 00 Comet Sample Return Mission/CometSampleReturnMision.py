import CoplanarOrbits as ob
import numpy
import matplotlib.pyplot as plt


'''
Script is built using my work in progress orbital solvers.
Solves the problem 3 times, one for each mass estimation

A few parameters are subject to change that are both
easily changable and may affect science instruments.


1. parking orbit around earth at launch (near bottom of script)
    change parkAltEarth
    currently 300 km
2. closest passing of mars (hopefully above atmosphere)
    change closestMars
    currently 20 km
3. interception points with comet 46P
    change frac1, frac2
    frac1 currently 2/3, spend 2/3 of time of science before apoapsis
    frac2 currently 1/3, spend 1/3 of time of science after apoapsis
4. return to earth parking orbit:
    change earthAltitudeArrival
    currently 300 km
    
'''


# constants
G = 6.674E-11
oneAU = 1.496E11
oneDay = 60*60*24 # seconds


# orbital parameters

muEarth = 3.986e14
aEarth = 1.000*149597870.70*1000
eEarth = 0 # force circular orbit
rEarth = 6.378E6




muMars = 4.297e13 # intro to spaceflight by francis j hale
aMars = 1.524 # Au
aMars = aMars*oneAU # meters
rMars = 3.393e6 # meters



muSun = 1.327E20


# comet 46 orbital parameters

a_approx_46P = 3.09270540*oneAU


radius_46 = [1000,700,500]
density_46 = [1000,500,200]
def volume(r):
    import numpy
    vol = 4/3*numpy.pi*r**3
    return vol

a_approx_46P = 3.09270540*oneAU # from wikipedia?  check source!
radiusApoapsis_46P = 5.129946*oneAU # same deal as above


mu_46 = []
c46b1 = []
c46b2 = []
for i in range(3):
    vol = volume(radius_46[i])
    mu_46.append(G* density_46[i] * vol)
    #c46b.append(ob.celestialBody(mu_46[i],mu_star = muSun,semimajorAxis = 5.129946*oneAU))
    c46b1.append(ob.celestialBody(mu_46[i],radius_46[1],muSun,a_approx_46P))
    ### DISTANCE SHOULD BE AT POINT OF INTERSECTION AND DEPARTURE
    print('For gravitaitonal parameter of ',mu_46[i],'m^3/s^2')
    print('the average SOI through its orbit is',c46b1[i].SOI)
    print('\n')




totDeltaV = []
totTOF = []
for c46b in c46b1: # run calculations for each case

    print('the following calculations are for gravitational parameter of ',c46b.mu,' m^3/s^2')
    print('\n')

    dv = 0
    tof = 0

    EarthOrbit = ob.Orbit(mu=muSun,semimajorAxis=aEarth,eccentricity=0) # model as spherical
    Comet46P = ob.Orbit(mu=muSun,semimajorAxis=a_approx_46P,radiusApoapsis=radiusApoapsis_46P)
    MarsOrbit = ob.Orbit(mu=muSun,semimajorAxis=aMars,eccentricity=0)


    Mars = ob.celestialBody(muMars,rMars,muSun,MarsOrbit.ra)
    Earth = ob.celestialBody(muEarth,rEarth,muSun,EarthOrbit.ra)
    
    
    # find rough distance traversed by comet 46P during time it takes to do science
    # 3-4 weeks to find place
    # 2-3 weeks for drilling
    # call that 7 weeks
    timeOfScience = 7*7*24*60*60 # seconds
    orbitalAngularVelocity = numpy.pi*2/Comet46P.period
    averageAngleTraversed = orbitalAngularVelocity*timeOfScience
    frac1 = 2/3
    frac2 = 1/3
    # intersectAt = numpy.pi - frac1*averageAngleTraversed # radians, how far before apoapsis do you intersect?
    departAt = numpy.pi + frac2*averageAngleTraversed
    #intersectAti = numpy.linspace(90,180,500)
    #intersectAti = numpy.deg2rad(intersectAti)
    # intersectAt = 3.0912264387025825
    #intersectAt = 3.000103076197092 # for direct transfer eccentricity = 0.7
    #departAt = 4.23844896023525
    #intersectAt = 3.71111111111
    #departAt = 4.0
    #intersectAt = 2.841944329785664
    #departAt = numpy.pi

    #intersectAt = numpy.pi
    #departAt = 3.5957798856495486
    #departAt = 3.631926274360186
    #intersectAt = 2.3482611754105527
    #departAt = 5.267721015110158
    intersectAt = 2.0943951
    departAt = 3.8079911

    

    

    '''
    ##########
    for 12000 points between 180 and 300 deg:
    min dif  0.04952669850411206 deg for depart angle  4.23844896023525 radians
                                                   242.84523710309193 deg
    at intersect angle  3.092544814288811 radians
    ##########
    '''


    ###################################
    #       ORDER OF OPERATIONS       #
    #################################33

    '''
    This script will solve the problem bakcwards, beginning with solar centric
    transfer to earth.  This will force arrival when the earth needs to be there
    further working backwards will need to take time of flight into account
    so that we can make sure the earth is where we need it to be.  The solution will
    be out of order and difficult to follow, but time is the most important variable

    oberth escape from earth parking orbit
    coast to mars
    mars flyby
    coast after flyby (or fire at end, unclear what tof will need)
    transfer after flyby to comet 46 p
    interception at 46 p
    science maneuvers: plane change, synchrosity is a to do
    send lander down
    oberth escape 46 p
    transfer from 46P back to earth
    Hyperbolic arrival + capture

    each part of the flight needs to have 4 calculations:

    tof
    deltaV
    angle traversed in comet46 reference plane
    angle traversed of earth
    '''

    ### comet 46 p arrival and departure orbits
    scComet1 = ob.addOrbitPosition(Comet46P,trueAnomaly=ob.angle(intersectAt,'r')) # comet 46 P orbit at point of intersection
    scComet2 = ob.addOrbitPosition(Comet46P,trueAnomaly=ob.angle(departAt,'r')) # comet 46 P point of leaving
    scComet2 = ob.secondQuadrant(scComet2)

    sevenWeeks = 7*7*24*60*60 # expected scienc tof
    tof_science = ob.TOF(scComet1,scComet2)
    print('tof science',tof_science)
    tof_diff = tof_science-sevenWeeks
    tof_percDiff = ob.percentDifference(tof_science,sevenWeeks)
    print('difference between expected time of science and actual is ',tof_diff,' seconds or ',tof_diff/oneDay,' days')
    print('percent difference is ',tof_percDiff,'%')
    print('we have ',tof_diff/oneDay,' extra days to do science')
    print('\n')
    tof = tof + tof_science






    ### solar centric transfer from comet46P to earth
    scReturnCruise = ob.Orbit(mu=muSun,radiusApoapsis=scComet2.r,radiusPeriapsis=EarthOrbit.ra)
    scReturnCruise1 = ob.addOrbitPosition(scReturnCruise,radius=scComet2.r)
    scReturnCruise2 = ob.addOrbitPosition(scReturnCruise,radius=EarthOrbit.ra)

    # since we're in the second half of the plane, true Anomaly should be greater
    # than 180 degrees
    scReturnCruise1 = ob.secondQuadrant(scReturnCruise1)
    #scReturnCruise2 = ob.changeOrbitPosition(scReturnCruise2,trueAnomaly=ob.angle(360,'d'))
    scReturnCruise2 = ob.secondQuadrant(scReturnCruise2)

    #tof_returnCruise = ob.TOF(scReturnCruise1,scReturnCruise2)
    tof_returnCruise = scReturnCruise.period/2
    tof = tof + tof_returnCruise

    thetaReturn = abs(scReturnCruise2.fi.angle-scComet2.fi.angle)
    dv_returnCruise = ob.lawOfCosines(scReturnCruise1.v,scComet2.v,thetaReturn)
    # instead of adding to dv, save it for deltaV leaving the comet
    vInfComet = dv_returnCruise
    
    returnCruiseOffset = ob.offsetAngle(scComet2,scReturnCruise1)
    nuReturnCruise46Ref = ob.nu46Ref(scReturnCruise1,scReturnCruise2,returnCruiseOffset)

    # 46 reference frame
    earthReturnCruiseEnd = nuReturnCruise46Ref[1] # force connection between s.c. and earth
    earthReturnCruiseStart = earthReturnCruiseEnd - (EarthOrbit.angularVelocity*tof_returnCruise) # work backwards

    title = 'return cruise'
    ob.dispParams(tof_returnCruise,0,0,title) # oberth in, flyby out, dv = 0
    







    ### earth position for science maneuvers:
    earthScienceEnd = earthReturnCruiseStart # science ends when return starts
    earthScienceStart = earthScienceEnd - (EarthOrbit.angularVelocity*tof_science) # work backwards

    NuScience46Ref = [intersectAt,departAt]


    ### solar centric transfer to mars
    earthTransfersEnd = earthScienceStart
    '''
    #earthTransfersStart = earthTransfersEnd - (EarthOrbit.angularVelocity*tof_trasnferTot)

    ### transfer from earth to mars after oberth
    # start going through math chronologically, check difference at end, see what parameters can be adjusted
    scCruise = ob.Orbit(mu=muSun,radiusPeriapsis=EarthOrbit.rp,radiusApoapsis=1.2*MarsOrbit.ra)
    # fast transfer
    scCruise1 = ob.addOrbitPosition(scCruise,radius=EarthOrbit.rp)
    scCruise2 = ob.addOrbitPosition(scCruise,radius=MarsOrbit.ra)

    tof_cruise = ob.TOF(scCruise1,scCruise2)
    tof = tof + tof_cruise
    tof_transferTot = tof_cruise

    thetaCruise = scCruise1.fi.angle - 0 # elevation angle zero for circular orbits
    dv_cruise = ob.lawOfCosines(scCruise1.v,EarthOrbit.va,thetaCruise)
    vInfEarth = dv_cruise

    title = 'transfer from earth to mars'
    ob.dispParams(tof_cruise,0,0,title) # oberth in, flyby out, dv = 0

    dNuCruise = scCruise2.nu.angle-scCruise1.nu.angle





    ### mars flyby
    # velocity is not collinear at intersection, flyby
    # because target is circular, angular positioning does not matter
    scTarget = ob.addOrbitPosition(MarsOrbit,trueAnomaly = ob.angle(0,'r')) # sun centric target, must both have position at intersection
    scSatellite = scCruise2 # sun centric transfer
    cbTarget = Mars
    direction = "counterClockwise"
    closestMars = 20 * 1000 + Mars.radius # 20 km, 20,000 meters + planet radius.
                                          # mars atmosphere starts at 11.1 km altitude according to https://en.wikipedia.org/wiki/Atmosphere_of_Mars
    flybyMars = ob.Flyby(scTarget,scSatellite,cbTarget,closestMars,direction)

    tof_flyby = flybyMars.tof
    tof = tof + tof_flyby
    tof_transferTot = tof_transferTot + tof_flyby

    dv = dv + 0 # no deltaV for flybys

    title = 'mars flyby'
    ob.dispParams(tof_flyby,0,0,title)

    dNuMarsFlyby = MarsOrbit.angularVelocity*tof_flyby


    scMarsFlybyEnd = ob.Orbit(mu=muSun,velocity=flybyMars.Vo,elevationAngle=ob.angle(flybyMars.betaOut,'r'),radius=EarthOrbit.ra)
    '''

    scEarthLaunch = ob.addOrbitPosition(EarthOrbit,trueAnomaly=ob.angle(0,'r'))

    tof_transferTot = 0

    ### oberth from earth to comet 46 P
    scDirectTransfer = ob.Orbit(mu=muSun,radiusApoapsis=scComet1.r,radiusPeriapsis=EarthOrbit.ra)
    scDirectTransfer1 = ob.addOrbitPosition(scDirectTransfer,radius=EarthOrbit.ra)
    scDirectTransfer2 = ob.addOrbitPosition(scDirectTransfer,radius=scComet1.r)

    tof_directTransfer = ob.TOF(scDirectTransfer1,scDirectTransfer2)
    tof = tof + tof_directTransfer
    tof_transferTot = tof_transferTot + tof_directTransfer

    thetaDirectTransfer = abs(scDirectTransfer1.fi.angle-scEarthLaunch.fi.angle)
    dv_directTransfer = ob.lawOfCosines(scDirectTransfer1.v,scEarthLaunch.v,thetaDirectTransfer)
    # dv = dv + dv_directTransfer
    vInfEarth = dv_directTransfer


    offsetDirectTransfer = ob.offsetAngle(scComet2,scDirectTransfer2)
    

    dNuDirectTransfer = scDirectTransfer2.nu.angle-scDirectTransfer1.nu.angle

    # arrive at angle - change in true anomaly for each transfer
    scOrbitStart = NuScience46Ref[0]  - dNuDirectTransfer
    # find where earth is after time of flight
    earthTransfersStart = earthTransfersEnd - (EarthOrbit.angularVelocity*tof_transferTot)

    scOrbitStart =ob.mod2Pi(scOrbitStart)
    earthTransfersStart = ob.mod2Pi(earthTransfersStart)

    # difference between where the spacecraft starts and where the earth starts
    # should be near zero or the solution does not technically converge
    dif = ob.angle(scOrbitStart - earthTransfersStart,'r')
    # dif.toDeg()
    # need to have difference between start and end within one SOI
    # one SOI radius is either 0.354152141601692 rad or
    # 0.5731999569324114 rad.


    distance = 2*EarthOrbit.ra*numpy.sin(dif.angle/2)
    print(distance < Earth.SOI)
    print(dif.angle)
    dif.toDeg()
    if dif.angle > 0.58: # distance within one SOI of earth
        raise ValueError('mission does not return to earth')


    title = 'earth cruise to 46P'
    ob.dispParams(tof_directTransfer,0,0,title) # oberth and hyperbolic capture.



    ### earth hyperbolic arrival
    earthArrival = ob.addOrbitPosition(EarthOrbit,trueAnomaly=ob.angle(0,'r'))
    earthAltitudeArrival = 300*1000 # 300 km
    rp = earthAltitudeArrival + Earth.radius
    earthCapture = ob.Flyby(earthArrival,scReturnCruise2,Earth,rp,'counterClockwise')

    earthPark = ob.Orbit(mu=muEarth,eccentricity=0,radiusPeriapsis=rp)

    # collinear firing at periapsis
    dv_earthArrival = earthCapture.flyby.vp - earthPark.vp
    dv = dv + dv_earthArrival
    
    tof_arrival = earthCapture.tof/2
    tof = tof+tof_arrival

    title = 'Earth Capture'
    ob.dispParams(tof_arrival,0,dv_earthArrival,title)

    dNuEarthCapture = EarthOrbit.angularVelocity*tof_arrival # change in solar-centric 46 ref frame
    earthEnd = earthReturnCruiseEnd + dNuEarthCapture




    ### earth oberth maneuver escape
    parkAltEarth = 300*1000 # 300,000 m, 300 km
    rParkEarth = parkAltEarth + Earth.radius

    parking = ob.Orbit(mu=muEarth,radiusPeriapsis=rParkEarth,eccentricity=0) # force circular

    VOberth = ob.OberthManeuver(muEarth,rParkEarth,vInfEarth)

    dv_oberth = abs(VOberth-parking.vp)
    # one firing to get you all the way to mars into the flyby

    dv = dv + dv_oberth

    oberth1 = ob.Orbit(mu=muEarth,radiusPeriapsis=rParkEarth,radius=rParkEarth,velocity=VOberth,velocityPeriapsis=VOberth,trueAnomaly=ob.angle(0,'r'))
    nu_inf = numpy.arccos(-1/oberth1.e) - 0.001 # to keep from dividing by zero
    oberth2 = ob.changeOrbitPosition(oberth1,radius = Earth.SOI,trueAnomaly= nu_inf)
    
    tof_oberth = ob.TOF(oberth1,oberth2)
    tof = tof + tof_oberth

    title = 'Earth Oberth escape maneuver from '+str(rParkEarth)+' meters from center of earth'
    ob.dispParams(tof_oberth,dv_oberth,0,title)







    ### comet 46P arrival
    scTarget = scComet2 # sun centric target, must both have position at intersection
    scSatellite = scDirectTransfer2 # sun centric transfer
    cbTarget = ob.celestialBody(c46b.mu,c46b.radius,muSun,scDirectTransfer2.r) # generate celestial body object at point of intersection
    direction = "counterClockwise"
    closest46 = cbTarget.SOI/6 # closest altitude, meters
    closest46 = closest46 + c46b.radius # 10 km, 10,000 meters + planet radius.

    print('altitude above 46P ',closest46)
                           
    flyby46 = ob.Flyby(scTarget,scSatellite,cbTarget,closest46,direction)

    undisturbedFalling = flyby46.flyby
    comet46StoppingPoint = ob.changeOrbitPosition(undisturbedFalling,trueAnomaly=ob.angle(0,'r'))



    comet46Park = ob.Orbit(mu=c46b.mu,eccentricity=0,radiusPeriapsis = comet46StoppingPoint.r)
    # stop at periapsis to begin with

    # since velocities are collinear at periapsis, simply subtract
    dv_46Capture = abs(comet46Park.vp-comet46StoppingPoint.vp)
    dv = dv + dv_46Capture

    tof_capture = flyby46.tof/2 # since we're doing half of a hyperbolic flyby, simply divide the tof by two 
    tof = tof+tof_capture

    title = 'Comet 46P capture into circular parking orbit'
    ob.dispParams(tof_capture,0,dv_46Capture,title)





    ### comet 46P departure
    comet46OberthVelocity = ob.OberthManeuver(c46b.mu,comet46Park.ra,vInfComet)
    comet46Oberth = ob.Orbit(mu=c46b.mu,radiusPeriapsis=comet46Park.ra,velocityPeriapsis=comet46OberthVelocity)
    
    dv_46Depart = abs(comet46Oberth.vp-comet46Park.vp)
    dv = dv + dv_46Depart

    nu_inf_46 = numpy.arccos(-1/comet46Oberth.e) - 0.001 
    oberthStart = ob.addOrbitPosition(comet46Oberth,trueAnomaly=ob.angle(0,'r'))
    oberthEnd = ob.addOrbitPosition(comet46Oberth,trueAnomaly=ob.angle(nu_inf_46,'r'))

    tof_46Depart = ob.TOF(oberthStart,oberthEnd)
    tof = tof + tof_46Depart

    title = 'Comet 46P Oberth Maneuver'
    ob.dispParams(tof_46Depart,dv_46Depart,0,title)






    ### science maneuvers
    ### 4 weeks to find landing place, 3 weeks for drilling

    print('orbital period around Comet 46P is ',comet46Park.period/oneDay,' days')
    print('\n')
    # this likely needs some field of view equations, since for the worst case
    # scenario you cannot simply run a whole period and then change planes
    # if you need faster peiods, simply lower the value of "closest46"


    # if you have a 90 degree field of view, you need one period for three angles:
    # starting angle, rotated 45 deg, rotated 90 deg (total)

    title = '45 deg plane change'
    dv_planeChange1 = ob.planeChangeCircle(comet46Park,ob.angle(45,'d'))
    ob.dispPlaneChange(title,dv_planeChange1)
    dv = dv + dv_planeChange1

    title = 'second 45 deg plane change'
    dv_planeChange2 = ob.planeChangeCircle(comet46Park,ob.angle(45,'d'))
    ob.dispPlaneChange(title,dv_planeChange2)
    dv = dv + dv_planeChange1




    


    ### SEND DRILL DOWN
    # fire once to put ourselves on trajectory to hit surface
    # fire once at surface to keep us level
    # firings assumed to be collinear
    # use object solver to specify velocity at surface
    # specify entry condition of elevation angle = 90 to be aligned with the radius
    # descent = ob.Orbit(mu=c46b.mu,radius=c46b.radius,velocity=0.1,radiusApoapsis=comet46Park.ra)
    finalDescentEnd = ob.Orbit(mu=c46b.mu,radiusApoapsis=comet46Park.ra,radius=c46b.radius,radiusPeriapsis=0.5*c46b.radius)
    finalDescentEnd = ob.secondQuadrant(finalDescentEnd)
    finalDescentStart = ob.changeOrbitPosition(finalDescentEnd,radius=finalDescentEnd.ra)
    finalDescentStart = ob.secondQuadrant(finalDescentStart)

    tof_descent = ob.TOF(finalDescentStart,finalDescentEnd)
    # fires collinearly to get on descent traj
    # fires collinearly to end descent traj
    dv_startDescent = abs(comet46Park.va-finalDescentStart.vp)
    dv_endDescent = finalDescentEnd.v - 0

    title = 'sending down the lander'
    ob.dispParams(tof_descent,dv_startDescent,dv_endDescent,title)

    # assume spacecraft mass of 1000 kg
    mSC = 1000
    FSeparate = (mSC*c46b.mu)/comet46Park.ra**2 # mu = G * mass of the planet
    FSurface = (mSC*c46b.mu)/c46b.radius**2
    
    print('at radius of ',comet46Park.ra,' meters, the lander will separate')
    print('the force experienced here will be approximately ', FSeparate,' Newtons')
    print('at the surface of the comet, roughly ',c46b.radius,' meters, the lander will feel ',FSurface,' Newtons')
    print('\n')



    print('\n')
    print('\n')
    print('\n')
    totDeltaV.append(dv)
    totTOF.append(tof)



# plotting

# earth departure

[xe,ye] = ob.plotOrbit(scDirectTransfer1,scDirectTransfer2,show=False,rotate=ob.angle(intersectAt + numpy.pi/2,'r'))
# comet 46 flight, science
rotSci = numpy.pi #numpy.pi
[xs,ys] = ob.plotOrbit(scComet1,scComet2,show=False,rotate=ob.angle(rotSci+ numpy.pi/2,'r'))
# transfer home
returnCruiseOffset.toRad()
rotTsfr = returnCruiseOffset.angle + numpy.pi
test2 = ob.changeOrbitPosition(scReturnCruise2,trueAnomaly=ob.angle(360,'d'))
[xt,yt] = ob.plotOrbit(scReturnCruise1,test2,show=False,rotate=ob.angle(departAt + numpy.pi/2,'r'))


c46orb1 = ob.addOrbitPosition(Comet46P,trueAnomaly=ob.angle(0,'d'))
c46orb2 = ob.addOrbitPosition(Comet46P,trueAnomaly=ob.angle(360,'d'))
[x46,y46] = ob.plotOrbit(c46orb1,c46orb2,show=False,rotate=ob.angle(270,'d'),N=2000)


e1 = ob.addOrbitPosition(EarthOrbit,trueAnomaly=ob.angle(0,'d'))
e2 = ob.addOrbitPosition(EarthOrbit,trueAnomaly=ob.angle(360,'d'))
[xp,yp] = ob.plotOrbit(e1,e2,show=False)

sr = 6.95700E8



#fig = plt.figure()
plt.figure()
plt.plot(xe,ye)
plt.plot(xs,ys)
plt.plot(xt,yt)
plt.plot(x46,y46)
plt.plot(xp,yp)
sr = numpy.linspace(3,sr)
for sr in sr:
    s = ob.Orbit(mu=muSun,eccentricity=0,radiusPeriapsis=sr)
    s1 = ob.addOrbitPosition(s,trueAnomaly=ob.angle(15,'d'))
    s2 = ob.addOrbitPosition(s,trueAnomaly=ob.angle(360+15,'d'))
    [xs,ys] = ob.plotOrbit(s1,s2,show=False)
    plt.plot(xs,ys)
#fig.set_facecolor('k')
plt.rcParams['axes.facecolor'] = 'k'
plt.axis('equal')
plt.show()





print(ob.metersToAU(EarthOrbit.ra))
print(ob.metersToAU(scDirectTransfer1.r))
print(ob.metersToAU(scDirectTransfer2.r))
print(ob.metersToAU(scComet1.r))
print(ob.metersToAU(scComet2.r))
print(ob.metersToAU(scReturnCruise1.r))
print(ob.metersToAU(scComet1.ra))




for i in range(0,3):
    print('total dV for gravitational parameter ', c46b1[i].mu,' is')
    print('m/s: ',totDeltaV[i])
    print('km/s: ',totDeltaV[i]/1000)
    print('total time for same gravitational parameter')
    print('seconds: ',totTOF[i])
    print('days: ',totTOF[i]/oneDay)
    print('years: ',totTOF[i]/(oneDay*365))
    print('\n')

