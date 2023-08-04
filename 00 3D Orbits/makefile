make: OrbitsTest.o Orbit.o SanityFunctions.o VectorMath.o Transfers.o
	g++ OrbitsTest.o Orbit.o SanityFunctions.o VectorMath.o -o OrbitsTest

OrbitsTest.o: OrbitsTest.cpp Orbit.cpp Orbit.hpp OrbitPosition.cpp OrbitPosition.hpp SanityFunctions.cpp SanityFunctions.hpp
	g++ -c OrbitsTest.cpp

Orbit.o: Orbit.cpp Orbit.hpp SanityFunctions.cpp SanityFunctions.hpp VectorMath.cpp VectorMath.hpp
	g++ -c Orbit.cpp

Transfers.o: Transfers.cpp Transfers.hpp Orbit.cpp Orbit.hpp SanityFunctions.cpp SanityFunctions.hpp VectorMath.cpp VectorMath.hpp
	g++ -c Transfers.cpp

SanityFunctions.o: SanityFunctions.cpp SanityFunctions.hpp
	g++ -c SanityFunctions.cpp

VectorMath.o: VectorMath.cpp VectorMath.hpp SanityFunctions.cpp SanityFunctions.hpp
	g++ -c VectorMath.cpp

cleanObj:
	rm *.o

clean:
	rm *.o OrbitsTest