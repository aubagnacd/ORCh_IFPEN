#ifndef conditions_H
#define conditions_H

#include "../cantera/flamemodel.h"
#include "../read_write/QSSscenario.h"
#include "../optimisation/OptimScenario.h"

struct ORChInputs {
   int debuglevel;
   int IterationNumber; // number of iterations to run
   string mech; //Current step reference chemical scheme
   string mech_ref; //Reference detailed mechanism
   string mech_desc; //Description of the reference chemical scheme
   string configuration; //Studied combustion regime - OPTIONS: "MultipleInlet"; "PremixedFlames";
   string step; //Step to perform - OPTIONS: "DRGEP_Species"; "DRGEP_Reactions"; "ComputeTrajectories"; "computeQSSCriteria"; "getQSSfile"; "getQSSfileFORTRAN"; "Optimisation"; "Lumping";
   bool new_mixing; //Define if a new mixing of the particles is defined or if the mixing used with the previous step is kept
   bool activateCurl; // mixing model (true --> CURL, false --> emst)
   bool writeTraj; // write inlets trajectories data (only for MultipleInlet)
   bool writeAllPart; // write all particles data (only for MultipleInlet)
   bool drgepTraj; // do DRGEP on trajectories
   bool print_all_rAB; // print species inter-relations diagrams for all inlets/particles (depends on drgepTraj)
   double MixingTime;
   double TimeStep;
   double rMassFlowRate;
};

void conditions(
      ORChInputs &inputs,
      vector<MultipleInlet*> &listInlets, 
      vector<PremixedFlames*> &listFlames, 
      vector<AutoIgnition*> &listIgnitions, 
      vector<string> &listTargets, 
      vector<QSSscenario*> &listQSSscenarios, 
      OptimScenario* &listOptimScenarios,
      vector<string> &trajectory_ref);
#endif
