make: OrbitsTest.o Orbits.o SanityFunctions.o VectorMath.o Transfers.o
	g++ OrbitsTest.o Orbits.o SanityFunctions.o VectorMath.o -o OrbitsTest

OrbitsTest.o: OrbitsTest.cpp Orbits.cpp Orbits.hpp SanityFunctions.cpp SanityFunctions.hpp
	g++ -c OrbitsTest.cpp

Orbits.o: Orbits.cpp Orbits.hpp SanityFunctions.cpp SanityFunctions.hpp VectorMath.cpp VectorMath.hpp
	g++ -c Orbits.cpp

Transfers.o: Transfers.cpp Transfers.hpp Orbits.cpp Orbits.hpp SanityFunctions.cpp SanityFunctions.hpp VectorMath.cpp VectorMath.hpp
	g++ -c Transfers.cpp

SanityFunctions.o: SanityFunctions.cpp SanityFunctions.hpp
	g++ -c SanityFunctions.cpp

VectorMath.o: VectorMath.cpp VectorMath.hpp SanityFunctions.cpp SanityFunctions.hpp
	g++ -c VectorMath.cpp

cleanObj:
	rm *.o

clean:
	rm *.o OrbitsTest