#include "mainComputeTrajectories.h"

void ComputeTrajectories(ORChInputs inputs, vector<string> listTargets, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, vector<string> trajectory_ref, int rank)
{
   IdealGasMix *mixture  = new IdealGasMix(inputs.mech, inputs.mech_desc);

   int nbLines = 0;
   int nbInlets = 0;
   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size();
      nbLines =  inputs.IterationNumber;
   }

   int nsp_ref = mixture->nSpecies();
   int nreac_ref = mixture->nReactions();

   vector<vector<vector<double> > > Ym_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));
   vector<vector<double> > T_Trajectories_ref (nbInlets, vector<double> (nbLines, 0.0));
   vector<vector<vector<double> > > Production_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));
   vector<vector<vector<double> > > Consumption_Trajectories_ref (nbInlets, vector<vector<double> > (nbLines, vector<double> (nsp_ref, 0.0)));

   vector<vector<vector<double> > > R_AD_Trajectories (nbInlets, vector<vector<double> > (nsp_ref, vector<double>(nsp_ref, 0.0)));
   vector<vector<double> > max_j_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));
   vector<vector<double> > max_jf_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));
   vector<vector<double> > max_jr_on_Target (nsp_ref, vector<double> (nreac_ref,0.0));

   vector<bool> SpeciesIntoReactants (nsp_ref, false);

   vector<vector<double> > R_AD_Premixed (nsp_ref, vector<double>(nsp_ref, 0.0));
   vector<vector<double> > QSS_Criteria (nsp_ref, vector<double> (2, 0.0));

   vector<Species_ORCh*> listSpecies_ref;

   Read *r = new Read();
   r->Read_species(inputs.mech, listSpecies_ref);

   if (inputs.configuration == "MultipleInlet")
   {
      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(inputs, inputs.mech, listInlets, Targets, inputs.new_mixing, "ComputeTrajectories", R_AD_Trajectories, max_j_on_Target, Ym_Trajectories_ref, Production_Trajectories_ref, Consumption_Trajectories_ref, T_Trajectories_ref, SpeciesIntoReactants);
   }

   if (inputs.configuration == "PremixedFlames")
   {
      computePremixedFlames(inputs.mech, "./outputs/Premixed/Premixed_", inputs.mech_desc, listFlames, Targets, "ComputeTrajectories", SpeciesIntoReactants, R_AD_Premixed, max_j_on_Target, max_jf_on_Target, max_jr_on_Target, QSS_Criteria, false);
   }
}