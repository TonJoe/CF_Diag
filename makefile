all: LA1.h LA1.cpp CF.h CFN.cpp Monca.h Monca.cpp test.cpp
	g++ -o N14_DiagAll -g -Wall LA1.h LA1.cpp CF.h CFN.cpp Monca.h Monca.cpp test.cpp -I.
clean:
	
