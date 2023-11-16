all: conditions mainProgram

conditions : 
	mpic++ -c conditions.cpp

mainProgram:
	mpic++ -pthread -L$(SUNDIALS_LIB) -I$(SUNDIALS_INC) -L$(MPI_LIB) -I$(BOOSTPATH) -I$(ORCH_BASE)/Cantera/include -L$(CT_LIB) -L$(GTCOMB_CT_HOSTTYPE) -o mainProgram ../main/mainProgram.cpp conditions.o  ../tools/tools.o ../tools/gnuplot.o ../tools/fitness_criteria.o ../lumping/lumping.o ../lumping/lumping_species.o ../cantera/computeMultipleInlet.o ../cantera/computePremixedFlames.o ../cantera/computeAutoIgnition.o ../cantera/Analytic_function.o ../cantera/particle.o ../cantera/flamemodel.o ../drgep/drgep.o ../read_write/species.o ../read_write/reaction.o ../read_write/read.o ../read_write/write.o ../read_write/write_QSS.o ../read_write/write_QSS_FORTRAN.o ../read_write/QSSscenario.o ../optimisation/OptimScenario.o ../optimisation/optimisation.o ../optimisation/CreateAnalyticDirectory.o ../main/mainDRGEPSpecies.o ../main/mainDRGEPReactions.o ../main/mainComputeTrajectories.o ../main/mainLumping.o ../tools/outputs.o ../emst/emst_subs.o ../emst/emst.o -O3 -lm -lcantera -luser -lsundials_cvodes -lsundials_nvecserial -lgfortran

clean:
	rm -f mainProgram *.o
