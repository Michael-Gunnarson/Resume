built on a windows machine, I am unsure how to compile on other systems

Author:  Michael Gunnarson

OrbitsTest is an algorithm demonstration.  It takes the orbit of the Earth
around the sun and a satellite in Geostationary orbit and finds their orbital 
parameters.  The strength of this algorithm is that you only need a few well
chosen parameters to find the rest.  by setting a few fields (gravitational 
parameter field mu required) and calling the fill() method, it will automatically
solve for the rest of the orbital parameters.  These are printed to the command line.


compiling on my windows machine:

if you have the make command, cd to the folder the code lives in and type "make"
to get OrbitsTest.exe to compile.  Type ".\OrbitsTest.exe" to run. typing 
"make cleanObj" removes .o files. typing "make clean" removes .o files and OrbitsTest.exe.  

if you do not have the make command, use the following command (assuming you have g++):
"g++ OrbitsTest.cpp CoplanarOrbits.cpp SanityFunctions.cpp -o OrbitsTest.exe; .\OrbitsTest.exe"

if you are like me and had trouble figuring out how to compile things, I used
Cygwin64 to get g++, gcc, f2c, make, and other useful packages on my windwos machine.