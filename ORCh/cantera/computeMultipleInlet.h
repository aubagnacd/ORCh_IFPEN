#ifndef computeMultipleInlet_H
#define computeMultipleInlet_H

//General properties
#include <mpi.h>
#include <iostream>
#include <complex>
#include <time.h>
#include <ctime>
#include <fstream>
#include "mpi.h"

//ORCh
#include "flamemodel.h"
#include "../drgep/drgep.h"
#include "particle.h"
#include "../read_write/read.h"
#include "../main/conditions.h"

//Cantera properties
#include <Cantera.h>
#include <user.h>

using namespace Cantera;
using namespace Cantera_CXX;
using namespace User;

using std::cout;
using std::endl;
using std::string;

#define PI 3.14159265359

//--------------------
class computeMultipleInlet
{
   public:
   //constructeur
   computeMultipleInlet();

   virtual void getMultipleInlet(ORChInputs inputs, string mech, vector<MultipleInlet*> listInlets, vector<bool> Targets,
                          bool new_mixing, string step, vector<vector<vector<double> > > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target,
                          vector<vector<vector<double> > > &Ym_Trajectories_store, vector<vector<vector<double> > > &Production_Trajectories_ref,
                          vector<vector<vector<double> > > &Consumption_Trajectories_ref, vector<vector<double> > &T_Trajectories_store,
                          vector<bool> &SpeciesIntoReactants);

   virtual void Next_Time_Step_with_drgep(vector<bool> Targets, double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                          vector<vector<double> > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string step, bool print_all_rAB);

   virtual void Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t);

   virtual void Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                    vector<vector<vector<double> > > &Production_Trajectories_ref, vector<vector<vector<double> > > &Consumption_Trajectories_ref, int nInlet, int nLine);

   virtual void getMixedGasesComposition(vector<MultipleInlet*> listInlets, string step);

   virtual void Reacting(vector<Particle*> &listParticles, int nsp, double dt, double Pressure);
   virtual void ReactingParallel(vector<Particle*> &listParticles, int nsp, double dt, double Pressure);
   virtual void ReactingParallelDRGEP(vector<bool> Targets, vector<Particle*> &listParticles, int nsp, double dt, double Pressure, vector<vector<vector<double> > >&R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string  step);

   //destructeur
   virtual ~computeMultipleInlet();

   private:
	int Ifi_rank, Ila_rank, nb_var_loc, nproc;
	int *RecvCounts, *Disp;

	IdealGasMix *mixture;
   drgep *local_drgep;
};


#endif

