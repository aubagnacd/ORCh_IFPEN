

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include mainDRGEPSpecies.cpp mainDRGEPReactions.cpp mainComputeTrajectories.cpp mainLumping.cpp -c -O3

clean:
	rm *.o
