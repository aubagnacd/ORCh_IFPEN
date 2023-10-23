

all: Main

Main: 
	mpicxx -L$(MPI_LIB) -I$(BOOSTPATH) species.cpp reaction.cpp read.cpp write.cpp write_QSS.cpp write_QSS_FORTRAN.cpp QSSscenario.cpp -c -g -O3 -Wall -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -lconverters -lctbase -lctcxx -lcvode -lctf2c -lctf2c -lctlapack -lctblas -lm

clean:
	rm *.o
