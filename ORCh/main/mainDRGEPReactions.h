#ifndef main_drgep_reactions_H
#define main_drgep_reactions_H

//General properties
#include <mpi.h>
#include <algorithm>

//Tools
#include "../tools/tools.h"
#include "../tools/gnuplot.h"
#include "../tools/outputs.h"

//Conditions for the reduction 
#include "conditions.h"

//Compute flames (and DRGEP analysis)
#include "../cantera/computeMultipleInlet.h"
#include "../cantera/computePremixedFlames.h"
#include "../cantera/computeAutoIgnition.h"
#include "../cantera/Analytic_function.h"

//Read and Write mechanisms
#include "../read_write/read.h"
#include "../read_write/write.h"
#include "../read_write/write_QSS.h"
#include "../read_write/write_QSS_FORTRAN.h"

//Cantera properties
#include <Cantera.h>
#include <IdealGasMix.h>
#include <equilibrium.h>
#include <transport.h>
#include <zerodim.h>
#include <user.h>
using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;

   void drgepReactions(ORChInputs inputs,
                       vector<string> listTargets,
                       vector<MultipleInlet*> listInlets, 
                       vector<PremixedFlames*> listFlames,
                       vector<bool> Targets,
                       vector<string> trajectory_ref,
                       int rank);

#endif
