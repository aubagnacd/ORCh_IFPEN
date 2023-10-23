

all: Main

Main: 
	mpifort -w -std=c++11 -L$(MPI_LIB) -I$(BOOSTPATH) -I$(CT_INC) -L$(CT_LIB) -L$(GTCOMB_CT_HOSTTYPE) -I$(HDF5_INC) -L$(HDF5_LIB) emst_subs.f emst.f -c -g -O3 -Wall -lcantera -lm -lhdf5 -lhdf5_hl

clean:
	rm *.o *.mod
