

all: Main

Main: 
	mpic++ -I$(ORCH_BASE)/Cantera/include computeMultipleInlet.cpp computePremixedFlames.cpp computeAutoIgnition.cpp Analytic_function.cpp particle.cpp flamemodel.cpp -c -O3
clean:
	rm *.o
