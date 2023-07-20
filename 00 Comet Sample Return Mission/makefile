OrbitsTest: OrbitsTest.o CoplanarOrbits.o SanityFunctions.o
	g++ OrbitsTest.o CoplanarOrbits.o SanityFunctions.o -o OrbitsTest

OrbitsTest.o: OrbitsTest.cpp
	g++ -c OrbitsTest.cpp

CoplanarOrbits.o: CoplanarOrbits.cpp CoplanarOrbits.hpp SanityFunctions.hpp
	g++ -c CoplanarOrbits.cpp

SanityFunctions.o: SanityFunctions.cpp SanityFunctions.hpp
	g++ -c SanityFunctions.cpp

cleanObj:
	rm *.o

clean:
	rm *.o OrbitsTest