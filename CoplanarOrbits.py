'''

Script written by Michael Gunnarson for Professor Joshi's EAE 143 B space mission design
class, May-June 2023.

At present, the script is designed to look at two-body mechanics.  It contains several
classes to define orbit types, celestial bodies, and some functions to transfer from one
to another.  The upside to this script is the orbit class allows you to
fully define an orbit based on only a few parameters.  You do not need to know which
parameters should go to which function, you only need to name it correctly

This script has a few limitations.  Primarily, rendezvous are not yet implemented, and will
likely be implemented with another script due to differences in mathematical models.  Additionally,
the math used in this script use the pathced conic assumption, or successive two-body problems,
and assumes instantaneous changes in velocity as well as coplanar orbits.



the Orbit object requires keword arguments (kwargs) of any of the following forms:
all angles in a tuple with value and denotation of 'degrees' or 'radians'

specificEnergy
eccentricity
semilatusRectum
semimajorAxis
angularMomentum
radiusPeriapsis
radiusApoapsis
velocityPeriapsis
velocityApoapsis
radius
trueAnomaly
velocity
turningAngle # for hyperbolic and parabolic trajectories only
elevationAngle
escapeVelocity

the celestialBody object requires the gravitational parameter (called mu in this script) for the planet you want to orbit,
as well as the next largest one.  for example, the earth requires mu_earth as well as mu_sun in order to define the sphere of influence.
If you are orbiting the sun, or the SOI is not of great importance, simply



May 26, 2023: to do

add celestial body lookup table

Center calculations and input arguments that user interfaces with
around that.

figure out flyby implementation, how to end up where you're going

make sure flyby has offset calculations

electric propulsion



'''






### I'm sure it has a name, but it is escaping me at the moment
# if you assume coplanar motion for a more complicated trajectory,
# you need to have some way to keep track of how rotated different ellipses are
# they all rotate around the sun (or solar-system barycenter to be precise)
# so if you have some kind of arbitrary universal radial positioning coordinate
# you can deal with rotated ellipses
# as unintuitive as it is in my head, it would need to point to the periapsis
# simply because open trajectories do not have apoapses

# will have to check how to rotate geometry in arbitrary reference frames,
# though thankflly it is 2-D and no angular velocity will be needed.


# jk it was easy

def rotateGeometry(x,y,theta):
    # wikipedia: rotation matrix:
    # [cos(theta),-sin(theta);sin(theta),cos(theta)]*[x;y]
    import numpy
    xP = x*numpy.cos(theta) - y*numpy.sin(theta)
    yP = x*numpy.sin(theta) + y*numpy.cos(theta)
    return xP,yP


def mod360(d):
    try:
        for i in range(len(d)):
            while d[i] < 0:
                d[i] = d[i] + 360
            while d[i] > 360:
                d[i] = d[i] - 360
        return d
    except:
        while d < 0:
            d = d[i] + 360
        while d > 360:
            d = d - 360
        return d
    

def metersToAU(m):
    oneAU = 1.496E11 # meters
    return m/oneAU

def AUToMeters(AU):
    oneAU = 1.496E11 # meters
    return AU*oneAU


def secondQuadrant(orbit):
    import numpy
    orbit.nu.toRad()
    if orbit.nu.angle < numpy.pi:
        orbit = changeOrbitPosition(orbit,trueAnomaly=angle(orbit.nu.angle + numpy.pi,'r'))
    return orbit

def offsetAngle(c46,spaceCraft):
    '''
    must point to same position, must be same reference frame
    '''
    c46.nu.toRad()
    spaceCraft.nu.toRad()
    offset = angle(c46.nu.angle - spaceCraft.nu.angle,'r')
    return offset

def nu46Ref(orbitStart,orbitEnd,offset):
    orbitStart.nu.toRad()
    orbitEnd.nu.toRad()
    offset.toRad()
    nu46Start = orbitStart.nu.angle + offset.angle
    nu46End = orbitEnd.nu.angle + offset.angle
    return [nu46Start,nu46End]


def lawOfCosines(v1,v2,angleBetween):
    import numpy
    return (v1**2+v2**2-2*v1*v2*numpy.cos(angleBetween))**0.5

def percentDifference(expected,experimental):
    err = abs(expected-experimental)/expected * 100
    return err


def planeChangeCircle(orbit,angle):
    '''
    assumes circular orbit input (va = vp)
    returns deltaV for plane change
    '''
    angle.toRad()
    import numpy
    dv = 2*orbit.va *numpy.sin(angle.angle/2)
    return dv

def planeChange(orbit,angle):
    '''
    obrit must be in position configuration
    need position in orbit where you want plane change to occur
    I recommend firing at the apogee if at all possible, it will lower dV
    returns deltaV for plane change
    '''
    angle.toRad()
    import numpy
    dv = 2*orbit.v * numpy.sin(angle.angle/2)
    return dv

    
def mod2Pi(angle):
    import numpy
    try:
        for i in range(len(angle)):
            while angle[i] < 0:
                angle[i] = angle[i] + numpy.pi*2
            while angle[i] > numpy.pi*2:
                angle[i] = angle[i] - numpy.pi*2
        return angle
    except:
        while angle < 0:
            angle = angle + numpy.pi*2
        while angle > numpy.pi*2:
            angle = angle - numpy.pi*2
        return angle


def dispParams(tof,dv,dv2,title):
    oneDay = 60*60*24 # seconds
    print('time of flight for '+title+' is')
    print('seconds: ',tof)
    print('days: ', tof/oneDay)
    print('deltaV first burning for '+title+' is')
    print('m/s: ',dv)
    print('km/s: ',dv/1000)
    print('deltaV second burning for '+title+' is')
    print('m/s: ',dv2)
    print('km/s: ',dv2/1000)
    print('\n')

def dispPlaneChange(title,dv):
    print('delta v for '+title+' is')
    print('m/s: ',dv)
    print('km/s: ',dv/1000)
    print('\n')



def dictionaryToString(args):
    e = "="
    c = ","
    string = str()
    for key,value in args.items():
        string = string + key + e + str(value) + c
    string = string[:-1]
    return string

def orbitAssertRadians(orbit):
    try:
        orbit.fi.toRad()
        orbit.nu.toRad()
    except:
        pass

def changeOrbitPosition(orbit,**kwargs):
    '''
    kwargs should be new orbital position.  orbit
    must be in position configuration.  If it is not,
    consider using addOrbitPosition(orbit,**kwargs)
    '''
    kwargs = unwrapAngle(kwargs)
    orbit = removeOrbitPosition(orbit)
    orbit = addOrbitPosition(orbit,**kwargs)
    return orbit

def addOrbitPosition(orbit,**kwargs):
    '''
    takes orbit object of standard configuration and new keword arguments to
    add a position.  kwargs must be radius, elevationAngle, trueAnomaly,
    velocity, and/or escapeVelocity
    '''
    kwargs = unwrapAngle(kwargs)
    orbitArgs = stripOrbit(orbit)
    posArgs = dictionaryToString(kwargs)
    evalString = 'Orbit('
    o = eval(evalString + orbitArgs + ',' + posArgs + ')')
    return o


def stripOrbit(orbit):
    orbitArgs = orbit._args
    return orbitArgs

def removeOrbitPosition(orbit):
    '''
    i may have missed something but this should
    remove orbital positioning data, position flag,
    as well as edit kwargs and args.  Returns different object
    than one sent to function
    '''
    orbit = copy(orbit)
    orbit._positionFlag = False
    del orbit.r
    del orbit.v
    del orbit.nu
    del orbit.fi
    del orbit.vesc
    k = orbit._kwargs
    del k['radius']
    del k['velocity']
    del k['trueAnomaly']
    del k['elevationAngle']
    del k['escapeVelocity']
    args = dictionaryToString(k)
    orbit._args = args
    return orbit



class angle:
    def __init__(self,angle,angleType):
        '''
        really simple class to denote angle and type
        expects angle type as 'degrees' or 'radians', 'r', or 'd'
        '''
        if (angleType == 'degrees' or angleType == 'radians' or angleType == 'r' or angleType == 'd')==0:
            raise ValueError("angleType must be one of the following four strings: 'degrees', 'radians', 'r', or 'd'")
        if angleType == 'r':
            angleType = 'radians'
        elif angleType == 'd':
            angleType = 'degrees'
            
        self.angle = angle
        self.angleType = angleType

    def toRad(self):
        import numpy
        if self.angleType == 'radians':
            pass
        else:
            self.angle = numpy.deg2rad(self.angle)
            self.angleType = 'radians'
            
    def toDeg(self):
        import numpy
        if self.angleType == 'degrees':
            pass
        else:
            self.angle = numpy.rad2deg(self.angle)
            self.angleType = 'degrees'
    def print(self):
        print(self.angle,' '+self.angleType)

'''
def checkAngle(a,suppress = False):
    try:
        b = a.angle
        c = a.angleType
    except:
        if suppress == False:
            raise ValueError('angle is not an angle object')
        else:
            return 0
'''

def unwrapAngle(kwargs):
    '''
    for given set of kwargs (in other class or function),
    find which kwarg is an angle object, translate to rads, and remove object
    '''
    for key,value in kwargs.items():
        # literally the jankiest solution i can imagine
        dummyAngle = angle(5,'d')
        if isinstance(value,type(dummyAngle)):
            value.toRad()
            kwargs[str(key)] = value.angle
    return kwargs


def copy(obj):
    import copy
    return copy.deepcopy(obj)

def cylindricalToCartesian(r,theta):
    import numpy
    x = r*numpy.cos(theta)
    y = r*numpy.sin(theta)
    return x,y


def plotOrbit(orbit1,orbit2,show=True,rotate=angle(90,'d'),N=100):
    '''
    orbit1 and orbit2 must be in position configuration for the same orbit
    Trajectory is from 1 to 2
    gives the option to show plot in case you want to chain trajectories together
    '''
    import matplotlib.pyplot as plt
    import numpy
    rotate.toRad()
    if orbit1.E == orbit2.E:

        def _r(e,nu,p):
            return p/(1+e*numpy.cos(nu))
        eccentricity = orbit1.e
        semilatusRectum = orbit1.p
        orbitAssertRadians(orbit1)
        orbitAssertRadians(orbit2)
        nu = numpy.linspace(orbit1.nu.angle,orbit2.nu.angle,N)
        radiusPlot = _r(eccentricity,nu,semilatusRectum)
        #plt.figure()
        x,y = cylindricalToCartesian(_r(eccentricity,nu,semilatusRectum),nu)
        rotate.toRad()
        x,y = rotateGeometry(x,y,rotate.angle)
        #plt.plot(x,y)
        if show == True:
            #plt.show()
            pass
        elif show == False:
            return x,y
        
    else:
        print("printing not possible, make sure orbit1 and orbit2 are positions that lie on the same trajectory")
        

    

class Orbit:
    '''
    the Orbit object requires keword arguments (kwargs) of any of the following forms:

    mu # gravitational parameter
    specificEnergy
    eccentricity
    semilatusRectum
    semimajorAxis
    angularMomentum
    radiusPeriapsis
    radiusApoapsis
    velocityPeriapsis
    velocityApoapsis
    radius
    trueAnomaly
    velocity
    turningAngle # for hyperbolic and parabolic trajectories only
    elevationAngle
    escapeVelocity

    it will automatically solve for the ones you don't input as long as the solution can converge.
    at present, the solver does not account for inclination or geometry beyond semilatus rectum and semimajor axis
    '''

    def __init__(self, **kwargs):
        '''
        defines internal orbital functions as classes with kwarg inputs to allow for polymorphism.  Afterwards, it loops
        through a solution algorithm, appending solutions to the kwarg set of inputs.  This allows for it to iteratively
        solve for all parameters even with many unknown
        '''            
        import numpy
        class specificEnergy:
            def __init__(self,**kwargs):
                if 'specificEnergy' in kwargs:
                    self.specificEnergy = kwargs['specificEnergy']
                elif 'velocity' in kwargs and 'mu' in kwargs and 'radius' in kwargs:
                    self.specificEnergy = kwargs['velocity']**2/2 - kwargs['mu']/kwargs['radius']
                elif 'mu' in kwargs and 'semimajorAxis' in kwargs:
                    self.specificEnergy = -kwargs['mu']/(2*kwargs['semimajorAxis'])
                elif 'velocityPeriapsis' in kwargs and 'mu' in kwargs and 'radiusPeriapsis' in kwargs:
                    self.specificEnergy = kwargs['velocityPeriapsis']**2/2 - kwargs['mu']/kwargs['radiusPeriapsis']
                elif 'velocityApoapsis' in kwargs and 'mu' in kwargs and 'radiusApoapsis' in kwargs:
                    self.specificEnergy = kwargs['velocityApoapsis']**2/2 - kwargs['mu']/kwargs['radiusApoapsis']
                    
              
        class eccentricity:
            def __init__(self,**kwargs):
                if 'eccentricity' in kwargs:
                    self.eccentricity = kwargs['eccentricity']
                elif 'radiusPeriapsis' in kwargs and 'radiusApoapsis' in kwargs:
                    self.eccentricity = (kwargs['radiusApoapsis'] -
                    kwargs['radiusPeriapsis'])/(kwargs['radiusApoapsis'] + kwargs['radiusPeriapsis'])
                elif 'specificEnergy' in kwargs and 'mu' in kwargs and 'angularMomentum' in kwargs:
                    self.eccentricity = (1 + 2*kwargs['specificEnergy']*kwargs['angularMomentum']**2
                                         /kwargs['mu']**2)**0.5

        class semilatusRectum:
            def __init__(self,**kwargs):
                if 'semilatusRectum' in kwargs:
                    self.semilatusRectum = kwargs['semilatusRectum']
                elif 'semimajorAxis' in kwargs and 'eccentricity' in kwargs:
                    self.semilatusRectum = kwargs['semimajorAxis']*(1-kwargs['eccentricity']**2)
                elif 'angularMomentum' in kwargs and 'mu' in kwargs:
                    self.semilatusRectum = kwargs['angularMomentum']**2/kwargs['mu']


        class semimajorAxis:
            def __init__(self,**kwargs):
                if 'semimajorAxis' in kwargs:
                    self.semimajorAxis = kwargs['semimajorAxis']
                elif 'radiusApoapsis' in kwargs and 'radiusPeriapsis' in kwargs:
                    self.semimajorAxis = (kwargs['radiusApoapsis'] + kwargs['radiusPeriapsis'])/2
                elif 'radiusApoapsis' in kwargs and 'eccentricity' in kwargs:
                    self.semimajorAxis = kwargs['radiusApoapsis']/(1 + kwargs['eccentricity'])
                elif 'radiusPeriapsis' in kwargs and 'eccentricity' in kwargs:
                    self.semimajorAxis = kwargs['radiusPeriapsis']/(1 - kwargs['eccentricity'])
                elif 'mu' in kwargs and 'specificEnergy' in kwargs:
                    self.semimajorAxis = -kwargs['mu']/(2*kwargs['specificEnergy'])


        class angularMomentum:
            def __init__(self, **kwargs):
                if 'angularMomentum' in kwargs:
                    self.angularMomentum = kwargs['angularMomentum']
                elif 'mu' in kwargs and 'semilatusRectum' in kwargs:
                    self.angularMomentum = (kwargs['mu']*kwargs['semilatusRectum'])**0.5
                elif 'radiusApoapsis' in kwargs and 'velocityApoapsis' in kwargs:
                    self.angularMomentum = kwargs['radiusApoapsis']*kwargs['velocityApoapsis']
                elif 'radiusPeriapsis' in kwargs and 'velocityPeriapsis' in kwargs:
                    self.angularMomentum = kwargs['radiusPeriapsis']*kwargs['velocityPeriapsis']
                elif 'radius' in kwargs and 'velocity' in kwargs and 'elevationAngle' in kwargs:
                    self.angularMomentum = kwargs['radius']*kwargs['velocity']*numpy.cos(kwargs['elevationAngle'])
         

        class radiusPeriapsis:
            def __init__(self,**kwargs):
                if 'radiusPeriapsis' in kwargs:
                    self.radiusPeriapsis = kwargs['radiusPeriapsis']
                elif 'semimajorAxis' in kwargs and 'eccentricity' in kwargs:
                    self.radiusPeriapsis = kwargs['semimajorAxis']*(1-kwargs['eccentricity'])
                elif 'semimajorAxis' in kwargs and 'radiusApoapsis' in kwargs:
                    self.radiusPeriapsis = 2*kwargs['semimajorAxis']-kwargs['radiusApoapsis']
                elif 'semilatusRectum' in kwargs and 'eccentricity' in kwargs:
                    self.radiusPeriapsis = kwargs['semilatusRectum']/(1+kwargs['eccentricity'])



        class radiusApoapsis:
            def __init__(self,**kwargs):
                if 'radiusApoapsis' in kwargs:
                    self.radiusApoapsis = kwargs['radiusApoapsis']
                elif 'semilatusRectum' in kwargs and 'eccentricity' in kwargs:
                    self.radiusApoapsis = kwargs['semilatusRectum']/(1-kwargs['eccentricity'])
                elif 'semimajorAxis' in kwargs and 'eccentricity' in kwargs:
                    self.radiusApoapsis = kwargs['semimajorAxis']*(1+kwargs['eccentricity'])
                elif 'semimajorAxis' in kwargs and 'radiusPeriapsis' in kwargs:
                    self.radiusApoapsis = 2*kwargs['semimajorAxis']-kwargs['radiusPeriapsis']


        class velocityPeriapsis:
            def __init__(self,**kwargs):
                if 'velocityPeriapsis' in kwargs:
                    self.velocityPeriapsis = kwargs['velocityPeriapsis']
                elif 'angularMomentum' in kwargs and 'radiusPeriapsis' in kwargs:
                    self.velocityPeriapsis = kwargs['angularMomentum']/kwargs['radiusPeriapsis']
                elif 'velocityApoapsis' in kwargs and 'radiusPeriapsis' in kwargs and 'radiusApoapsis' in kwargs:
                    self.velocityPeriapsis = kwargs['radiusApoapsis']*kwargs['velocityApoapsis']/kwargs['radiusPeriapsis']


        class velocityApoapsis:
            def __init__(self,**kwargs):
                if 'velocityApoapsis' in kwargs:
                    self.velocityApoapsis = kwargs['velocityApoapsis']
                elif 'angularMomentum' in kwargs and 'radiusApoapsis' in kwargs:
                    self.velocityApoapsis = kwargs['angularMomentum']/kwargs['radiusApoapsis']
                elif 'velocityPeriapsis' in kwargs and 'radiusPeriapsis' in kwargs and 'radiusApoapsis' in kwargs:
                    self.velocityApoapsis = kwargs['radiusPeriapsis']*kwargs['velocityPeriapsis']/kwargs['radiusApoapsis']


        # optional parameters, requires specific position in orbit
        # used for position configuration of orbit class
        class velocity:
            def __init__(self,**kwargs):
                if 'velocity' in kwargs:
                    self.velocity = kwargs['velocity']
                elif 'specificEnergy' in kwargs and 'mu' in kwargs and 'radius' in kwargs:
                    self.velocity = (2*(kwargs['specificEnergy']+kwargs['mu']/kwargs['radius']))**0.5

        class radius:
            def __init__(self, **kwargs):
                if 'radius' in kwargs:
                    self.radius = kwargs['radius']
                elif 'semilatusRectum' in kwargs and 'eccentricity' in kwargs and 'trueAnomaly' in kwargs:
                    self.radius = kwargs['semilatusRectum']/(1+kwargs['eccentricity']*numpy.cos(kwargs['trueAnomaly']))

        class trueAnomaly:
            def __init__(self,**kwargs):
                if 'trueAnomaly' in kwargs:
                    self.trueAnomaly = kwargs['trueAnomaly']
                elif 'radius' in kwargs and 'semilatusRectum' in kwargs and 'eccentricity' in kwargs:
                    cosu = 1/kwargs['eccentricity'] * (kwargs['semilatusRectum']/kwargs['radius'] - 1)
                    if (abs(cosu)-1) <= .02:
                        cosu = numpy.clip(cosu,-1,1) # round due to peculiarites of numpy's arccos
                    self.trueAnomaly = numpy.arccos(cosu)                    


        class elevationAngle:
            def __init__(self,**kwargs):
                if 'elevationAngle' in kwargs:
                    self.elevationAngle = kwargs['elevationAngle']
                elif 'angularMomentum' in kwargs and 'radius' in kwargs and 'velocity' in kwargs:
                    cosu = kwargs['angularMomentum']/(kwargs['radius']*kwargs['velocity'])
                    #if (abs(cosu)-1) <= .02:
                    #    cosu = numpy.clip(cosu,-1,1) # round due to peculiarites of numpy's arccos
                    self.elevationAngle = numpy.arccos(cosu)
                elif 'eccentricity' in kwargs and 'trueAnomaly' in kwargs:
                    self.elevationAngle = numpy.arctan2(kwargs['eccentricity']*numpy.sin(kwargs['trueAnomaly']),(1+kwargs['eccentricity']*numpy.cos(kwargs['trueAnomaly'])))
 

        class turningAngle:
            def __init__(self,**kwargs):
                if turningAngle in kwargs:
                    self.turningAngle = kwargs['turningAngle']
                elif 'eccentricity' in kwargs:
                    self.turningAngle = 2*numpy.sin(1/kwargs['eccentricity'])


### IMPLEMENT ESCAPE FOR CIRCULAR VELOCITIES
        class escapeVelocity:
            def __init__(self,**kwargs):
                if 'escapeVelocity' in kwargs:
                    self.escapeVelocity = kwargs['escapeVelocity']
                elif 'mu' in kwargs and 'radius' in kwargs:
                    self.escapeVelocity = (2*kwargs['mu']/kwargs['radius'])**0.5


        def checkObjectConfiguration_Position(lst):
            '''
            check if orbit object will be in standard configuration or position configuration
            '''
            try:
                flag = (hasattr(r,'radius') and hasattr(v,'velocity') and hasattr(nu,'trueAnomaly')
                and hasattr(fi,'elevationAngle') and hasattr(vesc,'escapeVelocity'))
            except:
                flag = False
            return flag
                    

        '''
        Iterative solver to find all parameters
        '''
        self._err=False
        loopFlag = False
        # loop through until each orbital attribute is filled
        i = 0
        n = 14*5 # number of attributes * 5, translates to 5 times through the loop
        positionFlag = False # initialize position configuration
            # default configuration is E,e,p,a,H,rp,ra,vp,va
            # may need updating in since ra and a are infinity for parabolic/hyperbolic traj
            # del self.ra or del self.a both work if type is para/hyper and you find this to be necessary to avoid user confusion or avoid errors 
            # position configuration adds r,v,nu,fi,and escape velocity
            ### add table to read output vars
        while loopFlag == False and i < n:
            kwargs = unwrapAngle(kwargs) # INTERNAL CONVERSION TO RADIANS
            args = dictionaryToString(kwargs)
            #print(args)
 
            E = eval('specificEnergy(' + args + ')')
            e = eval('eccentricity(' + args + ')')
            p = eval('semilatusRectum(' + args + ')')
            a = eval('semimajorAxis(' + args + ')')
            H = eval('angularMomentum(' + args + ')')
            rp = eval('radiusPeriapsis(' + args + ')')
            ra = eval('radiusApoapsis(' + args + ')')
            vp = eval('velocityPeriapsis(' + args + ')')
            va = eval('velocityApoapsis(' + args + ')')
            lst = [E,e,p,a,H,rp,ra,vp,va]

            
            # desired architecture for position elements:
            ''' if any(optionalParams):
                    ex = eval('example(' + args + ')')
                    lst.append(optionalParams)
                    '''
            if 'radius' in kwargs or 'velocity' in kwargs or 'trueAnomaly' in kwargs or 'elevationAngle' in kwargs or 'escapeVelocity' in kwargs:
                # update flag
                positionFlag = True
                # evaluate code
                r = eval('radius(' + args + ')')
                v = eval('velocity(' + args + ')')
                nu = eval('trueAnomaly(' + args + ')')
                fi = eval('elevationAngle(' + args + ')')
                vesc = eval('escapeVelocity(' + args + ')')
                # append to lst
                lst2 = [r,v,nu,fi,vesc]
                for element in lst2:
                    lst.append(element)
            
            i += 1
            
            for objct in lst:
                if vars(objct) != {}: # each class has only one attribute, check if it has been filled
                    temp_str = list(vars(objct).keys())[0] # store solved data name
                    if temp_str not in kwargs: # check if in kwargs
                        eval_string = "objct.pp" # dummy string
                        eval_string = eval_string.replace("pp",temp_str) # replace pp with attribute
                        kwargs[temp_str] = eval(eval_string)# append to kwargs



            # check if all orbital data has been filled:
            
            # find if standard orbit AND if position orbit        
            if (hasattr(E,'specificEnergy') and hasattr(e,'eccentricity') and hasattr(p,'semilatusRectum')
                and hasattr(a,'semimajorAxis') and hasattr(H,'angularMomentum') and hasattr(rp,'radiusPeriapsis')
                and hasattr(ra,'radiusApoapsis') and hasattr(vp,'velocityPeriapsis') and hasattr(va,'velocityApoapsis')
                and checkObjectConfiguration_Position(lst)): # hide position data in exception in case you don't have any
                
                loopFlag = True
                break
                

            # check standard orbit and not positin orbit
            elif ((hasattr(E,'specificEnergy') and hasattr(e,'eccentricity') and hasattr(p,'semilatusRectum')
                and hasattr(a,'semimajorAxis') and hasattr(H,'angularMomentum') and hasattr(rp,'radiusPeriapsis')
                and hasattr(ra,'radiusApoapsis') and hasattr(vp,'velocityPeriapsis') and hasattr(va,'velocityApoapsis'))
                and (positionFlag==False)):

                loopFlag = True
                break
            # continue looping until you've hit max number of loops
            if i == n:
                print("solution did not converge")
                self._err = True
                
            
        self.E = E.specificEnergy
        self.e = e.eccentricity
        self.p = p.semilatusRectum
        self.a = a.semimajorAxis
        self.H = H.angularMomentum
        self.rp = rp.radiusPeriapsis
        self.ra = ra.radiusApoapsis
        self.va = va.velocityApoapsis
        self.vp = vp.velocityPeriapsis
        self._mu = kwargs['mu'] # gravitational parameter of body that satellite orbits
        self._positionFlag = positionFlag
        self._args = args
        self._kwargs = kwargs
        

        ## check position flag: add optional implementation
        if positionFlag:
            self.r = r.radius
            self.v = v.velocity
            self.nu = angle(nu.trueAnomaly,'r')
            self.fi = angle(fi.elevationAngle,'r')
            self.vesc = vesc.escapeVelocity
        
        if self.E > 0:
            self.type = 'hyperbolic'
            self.turningAngle = eval('turningAngle(' + args + ')')
        elif numpy.isclose(self.e, 0) :
            self.type = 'circular'
            self.period = 2*numpy.pi*(self.a**3/kwargs['mu'])**0.5 # for circular orbits, a = r, use a so as not to trip position configuration of class
            self.angularVelocity = 2*numpy.pi/self.period
        elif numpy.isclose(self.E,0) or numpy.isinf(self.a) or numpy.isclose(self.e,1):
            self.type = 'parabolic'
            self.turningAngle = eval('turningAngle(' + args + ')')
        elif self.E < 0 and self.e < 1:
            self.type = 'elliptic'
            self.period = 2*numpy.pi*(self.a**3/kwargs['mu'])**0.5
        else:
            print(args)
            raise ValueError('type of orbit could not be determined')




### MUST MAKE HYPERBOLIC ORBIT OBJECT


'''
class inclination:
    # Velocity = (V_planet**2 + V_inf**2 - 2*V_planet*V_infinity*numpy.cos(pi-elevationAngle))**0.5
    pass
class Trajectory:
    pass
class Transfer:
    pass
class Hohmann:
    pass
class OberthManeuver: #pg 41
    # Vr = (2*mu/2 + v_inf**2)**0.5
    pass
class changePlaneAngle: #pg 50
    pass
'''

# chapter 5

### intercept problem:
# calc mean orbital angular velocity, omega: radians/period
# or angle traversed/time of flight if orbital position is known
# calc tof of transfer orbit
# theta = orbital angular velocity * tof
# theta is the angle the planet will move during transfer
# lead angle for hohmann transfer = 180 deg - theta
# my guess is 180 deg comes from trueAnomaly traversed by transfer orbit
# confirmed, 180 deg position is nu_final - nu_initial
# synodic period (time) between 2 objects earth and mars is
# 360 deg/abs(omega_earth-omega_mars)
# in summary:
# intercept problem gives position relationship between
# initial planetary orbit and target planet orbit
# as well as how often this occurs



class Intercept:
    def __init__(self,target,spacecraftTransferInitial):
        '''
        takes orbit of the target body and transfer orbit of the spacecraft
        returns object that holds lead angle, angle traversed by target body
        '''

        def findIntercept():
            '''
            the transfer orbit must intersect the target orbit
            find points of intersection
            flag more than one point
            return intercept object with time of flight, lead angle,
            and planet traverse angle
            return final orbit object if applicables
            '''

            # coarse trueAnomaly for 360 deg, sweep through both orbits
            # find local minima of difference between orbits
            # finer trueAnomaly for those corresponding regions
            # define trueAnomaly of intercept as min distances for each
            # local minima region for the fine true anomaly
            # use to find radius and true anomaly of crossing, no matter how skewed the orbits might be

            
            

        self.orbit = spacecraftPositionFinal
        self.tof = tof
        self.leadAngle = leadAngle
        self.thetaPlanet = thetaPlanet

    def unrwap(self,angleType = 'r'):
        if angleType == 'deg' or angleType == 'degrees' or angleType =='d' or angleType=='degs':
            print('Time of flight of transfer orbit: ',self.tof)
            print('Lead angle (deg) of target body at point of first firing: ',self.leadAngle.toDeg)
            print('Angle (deg) traversed by target body during this time: ',self.thetaPlanet.toDeg)
        elif angleType == 'rad' or angleType == 'radians' or angleType == 'r' or angleType == 'rads':
            print('Time of flight of transfer orbit: ',self.tof)
            print('Lead angle (rads) of target body at point of first firing: ',self.leadAngle.toRad())
            print('Angle (rads) traversed by target body during this time: ',self.thetaPlanet.toRad())


# left to do:
# derive angle difference between velocity vector of target and velcoity
# vector of spacecraft for flyby or intercept or other solution


def escapeVelocity(mu,r):
    return (2*mu/r)**0.5

def OberthManeuver(mu,r,v_inf):
    v_esc = escapeVelocity(mu,r)
    return (v_esc**2 + v_inf**2)**0.5


class infinityOrbit(Orbit):
    def __init__(self,celestialBody,suppress=True,**kwargs):
    ### REMOVE REQUIRED CELESTIAL BODY OR ADD IT TO ORBIT CLASS
        '''
        This class is an odd necessity.  Mu is definitely needed as a kwarg,
        but it has a different set of parameters to find energy.  The assumptn
        is that at the SOI, radius is approximately infinity since the total
        specific energy is almost identical.  This and the impact parameter
        allows for us to sneak our way into defining and orbit.  Impact params
        and d, approach distance, will both be calculatedto make sure you
        don't run into the planet you're trying to orbit.

        I took this class years ago, but I'm pretty sure you need the
        elevation angle in order to define this orbit.  Where it comes from,
        I'm not sure.  Probably geometry along solar centric orbit

        requires celestial body, and kwargs velocitySOI + elevationAngle
        or specifiec approach distance at a minimum
        '''
        import numpy
        # check infinity specific parameters
        kwargs = unwrapAngle(kwargs) # INTERNAL CONVERSION TO RADIANS
        kwargs['mu'] = celestialBody.mu


        if 'velocitySOI' in kwargs:
            kwargs['specificEnergy'] = kwargs['velocitySOI']**2/2 - celestialBody.mu/celestialBody.SOI
            # update positioning data
            # kwargs['velocity'] = kwargs['velocitySOI']
            # kwargs['radius'] = celestialBody.SOI
        if 'approachDistance' in kwargs and 'velocitySOI' in kwargs:
            kwargs['angularMomentum'] = kwargs['velocitySOI']*kwargs['approachDistance']
        if 'specificEnergy' in kwargs and 'angularMomentum' in kwargs:
            kwargs['eccentricity'] = (1 + 2*kwargs['specificEnergy']*kwargs['angularMomentum']**2
                                         /kwargs['mu']**2)**0.5
        if 'eccentricity' in kwargs:
            kwargs['trueAnomaly'] = numpy.arccos(-1/kwargs['eccentricity'])


        # create orbital object with characteristics
        
        super().__init__(**kwargs)
        
        #print(kwargs)
        #super().__init__()

        
        # add infinity specific parameters
        ### impact parameter
        ### approach distance
        if 'approachDistance' not in kwargs:
            self.d = celestialBody.SOI*numpy.cos(self.fi.angle) # approach distance
        else:
            self.d = kwargs['approachDistance']
        self.b = celestialBody.radius*(1+(2*celestialBody.mu/celestialBody.radius)/kwargs['velocitySOI']**2)**0.5 # impact parameter
        if suppress != True:
            if  self.d < self.b:
                print('trajectory will collide with surface of target')
            elif self.d == self.b:
                print('trajectory will graze surface of target')
            elif self.d > self.b:
                print('trajectory will fly by target')
        if self.d < self.b:
            self._collisionFlag = True
            self._flybyFlag = False
        elif self.d == self.b:
            self._collisionFlag = True
            self._flybyFlag = False
        if self.d > self.b:
            self._collisionFlag = False
            self._flybyFlag = True



# implement just the flyby.
# build the tool to find flyby based on multiple reference frames later
#def Flyby(celestialBodyTarget,infinityOrbit,sunCentricSatelliteOrbit,sunCentricTargetOrbit):
#    pass


### flyby pg 113

# V_inf_i = V_inf_o = V_inf


# sun-centric
# beta_i is the angle betweeen planet velocity and spacecraft velocity
# V_i is spacecraft velocity as it reaches planet
# V_planet is the velocity of the planet at point of contact




### GENERATE SOLVER FOR BETA_INITIAL
class Flyby:
    def __init__(self,scTarget,scSatellite,cbTarget,rp,direction,split = False):
        ### OPEN UP FLYBY TRAJECTORY TO ALLOW FOR ENTRANCE CONDITION CALCULATIONS 
        ### BUT NOT ACTUALLY FLY BY
        '''
        scTarget must have position configuration at first SOI entrance
        scSatellite must have same position config
        cbTarget
        tcSatellite must be an infinityOrbit object
        
        betaInitial is the angle (CCW positive) between the velocity of the planet
        and the initial velocity of the satellite as it enters SOI
        # if target is circular, beta_i is also elevationAngle.  I have a hunch
        # that beta is the difference between elevation angles
        # Francis Hale pg 114:
        " In this figure, Beta_in is the intersection angle between the planetary
        velocity V_planet and the incoming heliocentric velocity V_in and Xi_i
        is the angle between the planetary velocity and V_inf_i.  When Beta_i = 0,
        the two velocities are parallel (collinear) as would be the case with
        a Hohmann interplanetary transfer. (With a circular planetary orbit
        as shown, Beta_i is also the heliocentric elevation angle Fi_1.)

        
        naming convention: sc stands for sun or star-centric
        cb stands for celestial body
        ts stands for target cenric
        satellite is your spacecraft, naming convention from Francis J. Hale
        '''
        import numpy
        # trying to solve for betaOut
        # use geometry to figure out what path that puts you on???
        #if tcSatellite._flybyFlag == False:
        #    print('choose a different approach distance periapsis radius for target-centric flyby, your satellite currently crashes into target')
        #else:
        orbitAssertRadians(scSatellite)
        orbitAssertRadians(scTarget)

        betaIn = scSatellite.fi.angle - scTarget.fi.angle
        # target-centric entrance velocity at infinity
        V_inf = (scTarget.v**2 + scSatellite.v**2 - 2*scTarget.v*scSatellite.v*numpy.cos(betaIn))**0.5
        # xi_i is the angle between planetray velocity and V_infinity_i
        xi_i = -numpy.arcsin(scSatellite.v*numpy.sin(betaIn)/V_inf) + numpy.pi


        ### FIND TARGET CENTRIC SATELLITE ORBIT INSTEAD OF HAVE IT GIVEN?

        tcSatellite = infinityOrbit(cbTarget,velocitySOI=V_inf,radius=cbTarget.SOI,velocity=V_inf,radiusPeriapsis = rp)
        orbitAssertRadians(tcSatellite)
        turning_angle = 2*numpy.sin(1/tcSatellite.e)

        if direction == 'clockwise':
            xi_o = xi_i + turning_angle # clockwise rotation
        elif direction == 'counterClockwise':
            xi_o = xi_i - turning_angle # counter clockwise rotation

        V_o = (scTarget.v**2 + V_inf**2 + 2*scTarget.v*V_inf*numpy.cos(xi_o))**0.5

        beta_o = numpy.arcsin(V_inf/V_o * numpy.sin(xi_o))
        beta_i = betaIn
        if beta_o > beta_i:
            dirFlag = 'satellite turned away from the sun'
        elif beta_o < beta_i:
            dirFlag = 'satellite turned towards sun'



        def M(orbit):
            orbitAssertRadians(orbit)
            if orbit.type == 'elliptic':
                    cosu = (orbit.e + numpy.cos(orbit.nu.angle))/(1 + orbit.e*numpy.cos(orbit.nu.angle))
                    if (abs(cosu)-1) <= .1:
                        cosu = numpy.clip(cosu,-1,1) # round due to peculiarites of numpy's arccos
                    u = numpy.arccos(cosu) # eccentric anomaly
                    m = u - orbit.e*numpy.sin(u) # mean anomaly
            elif orbit.type == 'hyperbolic':
                coshF = (orbit.e+numpy.cos(orbit.nu.angle))/(1+orbit.e*numpy.cos(orbit.nu.angle))
                F = numpy.arccosh(coshF)# hypberbolic eccentric anomaly
                m = orbit.e*numpy.sinh(F)-F
            return m




        # tc satellite already defines the point of entrance at SOi
        tof = 2*((-tcSatellite.a)**3/tcSatellite._mu)**0.5 * (M(tcSatellite))
        # orbit defines one half of the hyperbola
        # multiply by two to get the time of flight
                                                              

        self.Vo = V_o
        self.betaOut = beta_o
        self.directionFlag = dirFlag
        self._target = scTarget
        self.flyby = tcSatellite
        self.tof = tof




"""
class postFlyby:
    def __init__(self,flybyObj):
        '''
        figure out how to plot new solar-centric trajectory
        eventually will need to find how to plot the actual flyby

        For separate function:
        assume hyperbolic traj symmetric around periapsis?  use to get
        entry and exit true anomaly to be able to plot the whole thing?
        '''
"""





### DESIGN LAMBERT SOLVER ITERATIVELY:
    # GIVEN TIME OF FLIGHT AND DISTANCE NEEDED TO TRAVERSE,
    # WHAT TRAJECTORY DO YOU NEED?
    # WHY?
    # idea: multiple mars fly bys, one there, one back.
    # if you can get the time at the end of the mission, and use spice
    # to find the pos of the planet, you can use lambert's prob
    # to get where you need to go in the time you need





### escape trajectory: does not need to be parallel with exit velocity
    # how do i calculate the case where it isn't? pg 97



    
### implement/figure out if all periapsis exist on a line
### implement conversion to and from spice


def newtonMethod(g,x_n,tol = 1e-6):
    # x_n+1 = x_n - f(x_n)/f'(x_n) where you want to find f(x) = 0
    # based on code written in ENG 180
    # x_np1 -> x(n+1), x_nm1 -> x(n-1)
    res = 1
    n = 0
    x_nm1 = abs(x_n-.1)

    while res > tol and n < 100001:

        dg = (g(x_n) - g(x_nm1))/(x_n-x_nm1)
        #print(dg)
        x_np1 = x_n - g(x_n)/dg
        

        res = abs(g(x_np1))
        print(x_np1)
        print(g(x_np1))
        print(res)
        

        x_nm1 = x_n
        x_n = x_np1
        n = n+1

    return x_n,res



def TOF(startOrbit,endOrbit):
    if startOrbit.type != endOrbit.type or startOrbit.E != endOrbit.E or startOrbit._mu != endOrbit._mu:
        raise ValueError('both orbit inputs should be same orbit trajectory (same specific energy, eccentricity, type, etc).  Set kewords radius and true anomaly to different values for otherwise same input parameters to denote start and end of time of flight')
    else:
        orbitAssertRadians(startOrbit)
        orbitAssertRadians(endOrbit)
        import numpy
        def M(orbit):
            orbitAssertRadians(orbit)
            if orbit.type == 'elliptic':
                    cosu = (orbit.e + numpy.cos(orbit.nu.angle))/(1 + orbit.e*numpy.cos(orbit.nu.angle))
                    if (abs(cosu)-1) <= .1:
                        cosu = numpy.clip(cosu,-1,1) # round due to peculiarites of numpy's arccos
                    u = numpy.arccos(cosu) # eccentric anomaly
                    if orbit.nu.angle > numpy.pi:
                        u = 2*numpy.pi - u
                    m = u - orbit.e*numpy.sin(u) # mean anomaly
            elif orbit.type == 'hyperbolic':
                coshF = (orbit.e+numpy.cos(orbit.nu.angle))/(1+orbit.e*numpy.cos(orbit.nu.angle))
                F = numpy.arccosh(coshF)# hypberbolic eccentric anomaly
                m = orbit.e*numpy.sinh(F)-F
            return m
        if startOrbit.type == 'elliptic':
            #if startOrbit.nu.angle == 0:
            #   tof = (startOrbit.a**3/startOrbit._mu)**0.5 * (M(endOrbit))
            #else:
            tof = (startOrbit.a**3/startOrbit._mu)**0.5 * (M(endOrbit)-M(startOrbit))
            
        elif startOrbit.type == 'parabolic':
            
            def parabolicTimeFromZero(orbit):
                return orbit.H**3/(2*orbit._mu**2)*(numpy.tan(orbit.nu.angle/2) + 1/3*(numpy.tan(orbit.nu.angle/2))**3)

            tof = parabolicTimeFromZero(endOrbit)-parabolicTimeFromZero(startOrbit)
        elif startOrbit.type == 'circular':
            tof = (startOrbit.a**3/startOrbit._mu)**0.5 * (M(endOrbit)-M(startOrbit))
        elif startOrbit.type == 'hyperbolic':
            if startOrbit.nu.angle ==0:
                tof = ((-startOrbit.a)**3/startOrbit._mu)**0.5 * (M(endOrbit))
            else:
                tof = ((-startOrbit.a)**3/startOrbit._mu)**0.5 * (M(endOrbit) - M(startOrbit))
        else: print('something went wrong')
    return tof



# CHECK IF HYPERBOLIC AND PARABOLIC ALLOW FOR DIFFERENCE IN TIMES AS WELL
# CHECK GEOMETRY OF JACKIE'S TRANSFER ORBIT



class celestialBody: # celestial body to orbit, "celstial body-centered system"
    '''
    the Celestial Body class has two configuratins.  For a planet, it will contain the gravitational
    parameter as well as the sphere of influence.  This requires inputting the gravitational parameter
    of the star it is orbiting as well as the distance between the two centers, averaged as semimajor axis
    for non-circular orbits.  If these values are left as None, it will default to the second configuration
    which is the configuration for a star.  It basically only has mu 
    '''
    def __init__(self,mu_orbitee,radius,mu_star, distance):
        self.mu = mu_orbitee
        def SOI(mu_small,mu_large,dist):
            return dist*(mu_small/mu_large)**0.4
        def SOI2(mu_small,mu_large,dist):
            return dist*(1/3*mu_small/mu_large)**(1/3)
        self.SOI = SOI(mu_orbitee,mu_star,distance)
        self.SOI2 = SOI2(mu_orbitee,mu_star,distance)
        self.radius = radius



# Implement ra, rp, and c into these equations
# Implement collision algorithm with inclination
# Implement circular, elliptical, transfer, hyperbolic, parabolic, escape
# Implement time of flight and flyby angle calculations



"""BALLPARK ESTIMATE BACK OF THE NAPKIN FERMI ESTIMATION"""
"""
muEarth = 3.986e14
aEarth = 1.000*149597870.70*1000
eEarth = 0.0167
muSun = 1.327E20
oneAU = 1.496E11
a_approx_46P = 3.09270540*oneAU
mu46P_density = 0.6*(1/1000)*100**3/1**3 # kg/m^3
mu46P_diameter = 1.4*1000
mu46P_volume = mu46P_diameter**2/4

G = 6.674E-11

mu46P = G*mu46P_volume*mu46P_density

#E = specificEnergy(mu = muEarth,semimajorAxis = aEarth)
#PP = specificEnergy(mu = muEarth)


EarthOrbit = Orbit(mu=muSun,semimajorAxis=aEarth,eccentricity=eEarth)
a = EarthOrbit




'''
Meeting with Gianluca on May 10, 2023
    Find ballpark deltaV for worst case no planetary flyby's.
    Estimate using solar centric transfer deltaV, multiply by 4.5
'''
Earth = celestialBody(muEarth,muSun,oneAU)
Sun = celestialBody(muSun)
Comet46P = celestialBody(mu46P,muSun,a_approx_46P,)

Comet46P_Orbit = Orbit(mu=muSun,semimajorAxis= a_approx_46P,radiusApoapsis = 5.129946*oneAU)

# leave earth at periapsis
# transfer orbit has rp = earth_rp and ra = 46p's apoapsis
transfer_orbit = Orbit(mu=muSun,radiusPeriapsis = EarthOrbit.rp, radiusApoapsis = Comet46P_Orbit.ra) 


dv1_skip_park = transfer_orbit.vp - EarthOrbit.vp
dv2 = abs(transfer_orbit.va - Comet46P_Orbit.va)

dv = 4.5*(dv1_skip_park + dv2) # dv to 46p, dv back, times two ish.  DOES NOT INCLUDE EXTENDED MISSION REQUIREMENTS


# TOF


'''
Meeting with Jackie 5/17/2023
    Find solar-centric transfer orbit
'''

muEarth = 3.986e14
aEarth = 1.000*149597870.70*1000
eEarth = 0.0167
muSun = 1.327E20
oneAU = 1.496E11
a_approx_46P = 3.09270540*oneAU # from wikipedia?  check source!
radiusApoapsis_46P = 5.129946*oneAU # same deal as above


# Solar-Centric Orbit of earth:
EarthOrbit = Orbit(mu=muSun,semimajorAxis=aEarth,eccentricity=eEarth)

# does not include about 10 degree inclination
# require geometry of difference between velocity vectors
Comet46P = Orbit(mu=muSun,semimajorAxis=a_approx_46P,radiusApoapsis=radiusApoapsis_46P)

# Transfer from apoapsis of earth, to apoapsis of 46P
# check alignment of periapsis and apoapsis, not sure if gemotry works.  If it
# doesn't, add about 800 m/s to deltaV1
TransferOrbit = Orbit(mu=muSun,radiusPeriapsis=EarthOrbit.ra,radiusApoapsis=Comet46P.ra)
# TransferOrbit = Orbit(mu=muSun,radiusPeriapsis=EarthOrbit.rp,radiusApoapsis=Comet46P.ra)

# difference between earth apoapsis velocity and transfer trajectory (at periapsis)
deltaV1 = abs(EarthOrbit.va-TransferOrbit.vp) # fire to get onto transfer orbit at periapsis
deltaV2 = abs(TransferOrbit.va-Comet46P.va) # fire to get out of transfer orbit at common apoapsis

print('delta V1 in m/s:',deltaV1)
print('delta V2 in m/s:',deltaV2)

'''
5/17/2023
    Positioning test
'''
trueAnomaly = angle(5,'d')
EarthOrbit = Orbit(mu=muSun,semimajorAxis=aEarth,eccentricity=eEarth)#,trueAnomaly=trueAnomaly)
"""
# heliocentric inertial frame:
# right argument of the ascending node?
