#include "drgep.h"

#include <sstream>
#include <fstream>

//---drgep---
drgep::drgep(IdealGasMix *mixture) //Constructeur
{
   nreac = mixture->nReactions();
   nsp = mixture->nSpecies();

   prodStoichCoeff.resize(nsp);
   reactStoichCoeff.resize(nsp);
   reac_list.resize(nsp);

   for (int k = 0; k < nsp; k++) {
      prodStoichCoeff[k].resize(nreac,0.0);
      reactStoichCoeff[k].resize(nreac,0.0);

      for (int r = 0; r < nreac; r++) {
         prodStoichCoeff[k][r] = mixture->productStoichCoeff(k,r);
         reactStoichCoeff[k][r] = mixture->reactantStoichCoeff(k,r);
         //
         if ((prodStoichCoeff[k][r] != 0.0) || (reactStoichCoeff[k][r]!=0.0)) reac_list[k].push_back(r);
      }
   }
}


//----------------------------------------
//   <drgep_0D_species> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for 0D (AutoIgnition or MultipleInlet) regimes to sort species according to their importance to the creation or destruction of the targets
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//         n: Number of the inlet
//         time: Time for the integration
//   OUT:
//         R_AD_Trajectories: DRGEP matrix of interactions between species
//----------------------------------------
void drgep::drgep_0D_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, int n, double time, bool print) const
{
   vector<double> fwdRates (nreac, 0.0);
   vector<double> revRates (nreac, 0.0);

	mixture->getFwdRatesOfProgress(&fwdRates[0]);
	mixture->getRevRatesOfProgress(&revRates[0]);

   //To store the data files that provide the interactions between species and get DRGEP graphs
   string outputName = "";
   if (print) {
      outputName = "output/rAB_inlet";
      stringstream s_nInlet;
      stringstream s_time;
      s_nInlet << n;
      s_time << time;
      outputName.append(s_nInlet.str()).append("_").append(s_time.str()).append(".dat");
   }
      
   drgep_species(mixture, Targets, R_AD_Trajectories, fwdRates, revRates, outputName); 
}


//----------------------------------------
//   <drgep_0D_reactions> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for 0D (AutoIgnition or MultipleInlet) regimes to sort reactions according to their importance to the creation or destruction of the targets
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//   OUT:
//         rj_for_k: Matrix providing the impact reaction j has on the creation/destruction of species k
//----------------------------------------
void drgep::drgep_0D_reactions(IdealGasMix *mixture, vector<vector<double> > &rj_for_k) const
{
   vector<double> fwdRates (nreac, 0.0);
   vector<double> revRates (nreac, 0.0);

	mixture->getFwdRatesOfProgress(&fwdRates[0]);
	mixture->getRevRatesOfProgress(&revRates[0]);

   vector<vector<double> > rjf_for_k (nsp, vector<double> (nreac,0.0));
   vector<vector<double> > rjr_for_k (nsp, vector<double> (nreac,0.0));

   drgep_reactions(mixture, rj_for_k, rjf_for_k, rjr_for_k, fwdRates, revRates);
}


//----------------------------------------
//   <drgep_1D_species> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for a 1D (Premixed) regime to sort species according to their importance to the creation or destruction of the targets
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         *flow: Pointer to the 1D flame characteristics
//         ino: Node of the 1D flame considered for the DRGEP analysis
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//   OUT:
//         R_AD_Trajectories: DRGEP matrix of interactions between species
//         position_max_wdot: Position of the maximum heat release (m) 
//----------------------------------------
void drgep::drgep_1D_species(IdealGasMix *mixture, StFlow* flow, int ino, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, double position_max_wdot) const
{
   vector<double> fwdRates (nreac, 0.0);
   vector<double> revRates (nreac, 0.0);

   for (int j=0; j<nreac; j++)
   {
      fwdRates[j] = flow->Rates_Progress_Fwd(j,ino);
      revRates[j] = flow->Rates_Progress_Rev(j,ino);
   }

   double position = flow->grid(ino);
   string outputName = "outputs/drgep/rAB_flame";
   stringstream s_nFlame;
   stringstream s_position;
   s_nFlame << 0; //BE CAREFUL: This must be adapted if various flames are computed (the variable nFlame should be added to the inlets)
   s_position << position-position_max_wdot;
   outputName.append(s_nFlame.str()).append("_").append(s_position.str()).append(".dat");

   drgep_species(mixture, Targets, R_AD_Trajectories, fwdRates, revRates, outputName); 
}


//----------------------------------------
//   <drgep_1D_reactions> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for a 1D (Premixed) regime to sort reactions according to their importance to the creation or destruction of the targets
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         *flow: Pointer to the 1D flame characteristics
//         ino: Node of the 1D flame considered for the DRGEP analysis
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//   OUT:
//         rj_for_k: Matrix providing the impact reaction j has on the creation/destruction of species k
//         rjf_for_k: Matrix providing the impact the forward part of reaction j has on the creation/destruction of species k
//         rjr_for_k: Matrix providing the impact the reverse part of reaction j has on the creation/destruction of species k
//----------------------------------------
void drgep::drgep_1D_reactions(IdealGasMix *mixture, StFlow* flow, int ino, vector<vector<double> > &rj_for_k, vector<vector<double> > &rjf_for_k, vector<vector<double> > &rjr_for_k) const
{
   vector<double> fwdRates (nreac, 0.0);
   vector<double> revRates (nreac, 0.0);

   for (int j=0; j<nreac; j++)
   {
      fwdRates[j] = flow->Rates_Progress_Fwd(j,ino);
      revRates[j] = flow->Rates_Progress_Rev(j,ino);
   }

   drgep_reactions(mixture, rj_for_k, rjf_for_k, rjr_for_k, fwdRates, revRates);
}


//----------------------------------------
//   <drgep_species> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for species
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         Targets: List of targets (species), (Targets[k] == true) if target and (Targets[k] == false) if not
//         fwdRates: Reactions forward rates (unit: ...)
//         revRates: Reactions reverse rates (unit: ...)
//         outputName: Name of the file in which the species inter-relations are written in order to plot relation graphs
//   OUT:
//         R_AD_Trajectories: DRGEP matrix of interactions between species
//----------------------------------------
void drgep::drgep_species(IdealGasMix *mixture, vector<bool> Targets, vector<vector<double> > &R_AD_Trajectories, vector<double> fwdRates, vector<double> revRates, string outputName) const
{
   int nel = mixture->nElements(); //Number of elements (typically: C, H, O, N, Ar)
   vector<double> Production_minus_consumption_k (nsp, 0.0);
   vector<double> Production_Atom_A (nel, 0.0);

   for (int k=0; k<nsp; k++)
   {
      double omega_k_prod = 0.0;
      double omega_k_cons = 0.0;
      int j;
      for (int rind=0; rind<reac_list[k].size(); rind++)
      {
         j = reac_list[k][rind];
         omega_k_prod += prodStoichCoeff[k][j]*fwdRates[j]
                        +reactStoichCoeff[k][j]*revRates[j];
         omega_k_cons += reactStoichCoeff[k][j]*fwdRates[j]
                        +prodStoichCoeff[k][j]*revRates[j];
      }

      Production_minus_consumption_k[k] = abs(omega_k_prod-omega_k_cons);

      for (int Atom=0; Atom<nel; Atom++) Production_Atom_A[Atom] += mixture->nAtoms(k,Atom)*(omega_k_prod);
   }

	for (int Atom=0; Atom<nel; Atom++) Production_Atom_A[Atom] = 1.0 / Production_Atom_A[Atom];
	
   vector<vector<double> > scaling_Atom_Target (nel, vector<double> (nsp));
   vector<double> scaling_Target (nsp, 0.0);

   for (int k=0; k<nsp; k++)
   {
      for (int Atom=0; Atom<nel; Atom++)
      {
			scaling_Atom_Target[Atom][k] = mixture->nAtoms(k,Atom)*Production_minus_consumption_k[k]*Production_Atom_A[Atom];

         if (scaling_Atom_Target[Atom][k] > scaling_Target[k]) scaling_Target[k] = scaling_Atom_Target[Atom][k];
      }
   }

   double r_AB_sup,r_AB_inf;
   vector<vector<double> > r_AB (nsp, vector<double>(nsp));
   vector<vector<double> > r_AD_intermediate (nsp, vector<double>(nsp,0.0));

   //r_AB quantifies the direct inter-relation between species 
   for (int ka=0; ka<nsp; ka++)
   {
      int j;
      r_AB_inf = 0.0;
      for (int rind=0; rind<reac_list[ka].size(); rind++)
      {
         j = reac_list[ka][rind];
         r_AB_inf += abs((reactStoichCoeff[ka][j]+prodStoichCoeff[ka][j])*(fwdRates[j]-revRates[j]));
      }

      if (r_AB_inf == 0.0)
      {
         for (int kb=0; kb>nsp;kb++) {
            r_AB[ka][kb] = 0.0;
            r_AD_intermediate[ka][kb] = r_AB[ka][kb];
         }
      } else {
         for (int kb=0; kb<nsp; kb++)
         {
            r_AB_sup = 0.0;

            for (int rind=0; rind<reac_list[kb].size(); rind++)
            {
               j = reac_list[kb][rind];
               r_AB_sup += abs((reactStoichCoeff[ka][j]+prodStoichCoeff[ka][j])*(fwdRates[j]-revRates[j]));
            }

            r_AB[ka][kb] = r_AB_sup/r_AB_inf;
            r_AD_intermediate[ka][kb] = r_AB[ka][kb];
         }
      }
   }

//----------Print the species inter-relations diagrams----------//
   if (outputName != "") {
      ofstream rAB(outputName.c_str());
      for (int ka=0; ka<nsp; ka++)
      {
         rAB << mixture->speciesName(ka) << "  ";
      }
      rAB << endl;

      for (int ka=0; ka<nsp; ka++)
      {
         for (int kb=0; kb<nsp; kb++)
         {
            rAB << r_AB[ka][kb] << "  ";
         }
         rAB << endl;
      }
      rAB.close();
   }
//--------------------------------------------------------------//

   vector<vector<double> > r_AD (nsp, vector<double>(nsp,0.0));

   int nb_interaction_level = 3;

   for (int interaction_level=2; interaction_level<2+nb_interaction_level; interaction_level++)
   {
      if (interaction_level>2) {
         for (int ka=0; ka<nsp; ka++)
         {
            for (int kb=0; kb<nsp; kb++)
            {
               r_AD_intermediate[ka][kb] = r_AD[ka][kb];
            }
         }
      }

      for (int ka=0; ka<nsp; ka++)
      {
         for (int kb=0; kb<nsp; kb++)
         {
            for (int ki=0; ki<nsp; ki++) //i : intermediate species
            {
               if (r_AD[ka][kb] < r_AD_intermediate[ka][ki]*r_AB[ki][kb]) r_AD[ka][kb] = r_AD_intermediate[ka][ki]*r_AB[ki][kb];
            }
         }
      }
   } //end interaction

   for (int ka=0; ka<nsp; ka++)
   {
      if (Targets[ka])
      {
         for (int kb=0; kb<nsp; kb++)
         {
            if (R_AD_Trajectories[ka][kb] < scaling_Target[ka]*r_AD[ka][kb])
            {
               R_AD_Trajectories[ka][kb] = scaling_Target[ka]*r_AD[ka][kb];
            }
         }
      }
   }
}


//----------------------------------------
//   <drgep_reactions> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis for species
//   IN: 
//         *mixture: Pointer to the composition of the mixture to get forward and reverse reaction rates
//         fwdRates: Reactions forward rates (unit: ...)
//         revRates: Reactions reverse rates (unit: ...)
//   OUT:
//         rj_for_k: Matrix providing the impact reaction j has on the creation/destruction of species k
//         rjf_for_k: Matrix providing the impact the forward part of reaction j has on the creation/destruction of species k
//         rjr_for_k: Matrix providing the impact the reverse part of reaction j has on the creation/destruction of species k
//----------------------------------------
void drgep::drgep_reactions(IdealGasMix *mixture, vector<vector<double> > &rj_for_k, vector<vector<double> > &rjf_for_k, vector<vector<double> > &rjr_for_k, vector<double> fwdRates, vector<double> revRates) const
{
   for (int k=0; k<nsp; k++)
   {
      double omega_k_prod = 0.0;
      double omega_k_cons = 0.0;
      for (int j=0; j<nreac; j++)
      {
         omega_k_prod += prodStoichCoeff[k][j]*fwdRates[j]
                        +reactStoichCoeff[k][j]*revRates[j];
         omega_k_cons += reactStoichCoeff[k][j]*fwdRates[j]
                        +prodStoichCoeff[k][j]*revRates[j];
      }

      double max_prod_cons;

      if (omega_k_prod > omega_k_cons)
         max_prod_cons = omega_k_prod;
      else if (omega_k_cons > omega_k_prod)
         max_prod_cons = omega_k_cons;


      for (int j=0; j<nreac; j++)
      {
         double get_rj_for_k;
         double get_rjf_for_k;
         double get_rjr_for_k;

         double omega_k_prod_jf = prodStoichCoeff[k][j]*fwdRates[j];
         double omega_k_prod_jr = reactStoichCoeff[k][j]*revRates[j];
         double omega_k_prod_j = omega_k_prod_jf+omega_k_prod_jr;

         double omega_k_cons_jf = reactStoichCoeff[k][j]*fwdRates[j];
         double omega_k_cons_jr = prodStoichCoeff[k][j]*revRates[j];
         double omega_k_cons_j = omega_k_cons_jf+omega_k_cons_jr;

         if (max_prod_cons > 1e-12)
         {
            if (omega_k_prod_j > omega_k_cons_j)
               get_rj_for_k = (omega_k_prod_j-omega_k_cons_j)/max_prod_cons;
            else if (omega_k_cons_j > omega_k_prod_j)
               get_rj_for_k = (omega_k_cons_j-omega_k_prod_j)/max_prod_cons;
            else
               get_rj_for_k = 0;
   
            if (omega_k_prod_jf > omega_k_cons_jf)
               get_rjf_for_k = (omega_k_prod_jf-omega_k_cons_jf)/max_prod_cons;
            else if (omega_k_cons_jf > omega_k_prod_jf)
               get_rjf_for_k = (omega_k_cons_jf-omega_k_prod_jf)/max_prod_cons;
            else
               get_rjf_for_k = 0;

            if (omega_k_prod_jr > omega_k_cons_jr)
               get_rjr_for_k = (omega_k_prod_jr-omega_k_cons_jr)/max_prod_cons;
            else if (omega_k_cons_jr > omega_k_prod_jr)
               get_rjr_for_k = (omega_k_cons_jr-omega_k_prod_jr)/max_prod_cons;
            else
               get_rjr_for_k = 0;
         }


         if (get_rj_for_k > rj_for_k[k][j])
            rj_for_k[k][j] = get_rj_for_k;

         if (get_rjf_for_k > rjf_for_k[k][j])
            rjf_for_k[k][j] = get_rjf_for_k;

         if (get_rjr_for_k > rjr_for_k[k][j])
            rjr_for_k[k][j] = get_rjr_for_k;
      }
   }
}

drgep::~drgep() //Destructeur
{}