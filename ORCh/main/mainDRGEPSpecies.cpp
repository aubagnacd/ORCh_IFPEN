#include "mainDRGEPSpecies.h"

//----------------------------------------
//   <mainDRGEPSpecies> 
//----------------------------------------
//   DESCRIPTION:
//         Run the DRGEP analysis by
//            1- Computing reference trajectories
//            2- Computing the species DRGEP coefficients (R_AD)
//            3- 
//   IN: 
//         configuration: Either "Premixed" or "MultipleInlet" case
//   OUT:
//
//----------------------------------------
void drgepSpecies(ORChInputs inputs, vector<MultipleInlet*> listInlets, vector<PremixedFlames*> listFlames, vector<bool> Targets, vector<string> trajectory_ref, int rank)
{
   IdealGasMix *mixture  = new IdealGasMix(inputs.mech,inputs.mech_desc);

   int nsp_init = mixture->nSpecies();
   int nreac_init = mixture->nReactions();
   int nbFlame = listFlames.size()-1;   

   stringstream s_nbFlame;
   s_nbFlame << nbFlame;

   vector<Species_ORCh*> listSpecies_init;
   vector<Reaction_ORCh*> listReactions_init;

   Read *r = new Read();
   r->Read_species(inputs.mech, listSpecies_init);
   r->Read_reactions(inputs.mech, listReactions_init);

   int nbIterations = 0;
   int nbInlets = 0;

   if (listInlets.size() > 1)
   {
      nbInlets = listInlets.size();
      nbIterations =  inputs.IterationNumber;
   }

   vector<vector<vector<double> > > R_AD_Trajectories;
   vector<vector<double> > R_AD_Premixed (nsp_init, vector<double>(nsp_init, 0.0));
   vector<vector<double> > max_j_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jf_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<double> > max_jr_on_Target (nsp_init, vector<double> (nreac_init,0.0));
   vector<vector<vector<double> > > Production_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<vector<double> > > Consumption_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations, vector<double> (nsp_init, 0.0)));
   vector<vector<double> > T_Trajectories_ref (nbInlets, vector<double> (nbIterations+1, 0.0));
   vector<vector<vector<double> > > Ym_Trajectories_ref (nbInlets, vector<vector<double> > (nbIterations+1, vector<double> (nsp_init, 0.0)));
   vector<bool> SpeciesIntoReactants (nsp_init, false);
   vector<vector<double> > QSS_Criteria (nsp_init, vector<double> (2, 0.0));

   double sort_R_AD [nsp_init][2]; //sort_R_AD[k][0] -> interaction coefficient //sort_R_AD[k][1] -> species number
   ofstream log("DRGEP_Species.log");

   if (inputs.configuration == "MultipleInlet")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = k;
      }

      computeMultipleInlet *c = new computeMultipleInlet();
      c->getMultipleInlet(inputs,
               inputs.mech,
               listInlets,
               Targets,
               inputs.new_mixing,
               inputs.step,
               R_AD_Trajectories, 
               max_j_on_Target, 
               Ym_Trajectories_ref, 
               Production_Trajectories_ref, 
               Consumption_Trajectories_ref, 
               T_Trajectories_ref,
               SpeciesIntoReactants);

      int nR_AD = (inputs.drgepTraj)?nbInlets:1;
      for (unsigned int n=0; n<nR_AD; n++)
      {
         for (int ka=0; ka<nsp_init; ka++)
         {
            if (Targets[ka])
            {
               for (int kb=0; kb<nsp_init; kb++)
               {
                  if (sort_R_AD[kb][0] < R_AD_Trajectories[n][ka][kb]) sort_R_AD[kb][0] = R_AD_Trajectories[n][ka][kb];
               }
            }
         }
      }
   } //end if (configuration == "MultipleInlet")

   if (inputs.configuration == "PremixedFlames")
   {
      for (int k=0; k<nsp_init; k++)
      {
         sort_R_AD[k][0] = 0.0;
         sort_R_AD[k][1] = k;
      }


      string outputName = "./outputs/Premixed/Ref_DRGEP_Species";
      stringstream s_nbSpeciesInit;
      s_nbSpeciesInit << nsp_init;
      outputName.append(s_nbSpeciesInit.str());

      computePremixedFlames(inputs.mech, 
                            outputName, 
                            inputs.mech_desc, 
                            listFlames, 
                            Targets, 
                            "DRGEP_Species", 
                            SpeciesIntoReactants, 
                            R_AD_Premixed, 
                            max_j_on_Target, 
                            max_jf_on_Target, 
                            max_jr_on_Target, 
                            QSS_Criteria, 
                            true);


      for (int ka=0; ka<nsp_init; ka++) //BE CAREFUL, FOR NOW DRGEP ANALYSIS ON PREMIXED FLAMES IS CODED FOR ONLY ONE FLAME, JUST NEED TO LOOP ON THE NFLAMES
      {
         if (Targets[ka])
         {
            for (int kb=0; kb<nsp_init; kb++)
            {
               if (sort_R_AD[kb][0] < R_AD_Premixed[ka][kb]) sort_R_AD[kb][0] = R_AD_Premixed[ka][kb];
            }
         }
      }
   } //end if (configuration == "PremixedFlames")

   //If the species is within the targets or within the initial reactants its R_AD coefficient is set to one to ensure it is kept while reducing the mechanism  
   for (int k=0; k<nsp_init; k++)
   {
      if (SpeciesIntoReactants[k]) sort_R_AD[k][0] = 1;
      if (Targets[k]) sort_R_AD[k][0] = 1;
   }

   qsort (sort_R_AD, sizeof(sort_R_AD)/sizeof(sort_R_AD[0]), sizeof(sort_R_AD[0]), compare_numbers);

   log << "mech = " << inputs.mech << endl;
   log << "mech_ref = " << inputs.mech_ref << endl;
   log << "trajectory_ref = " << trajectory_ref[0] << endl;
   log << "DRGEP coefficients ---------->" << endl; 


   if (rank ==0)
   {
      cout << endl  << endl << endl << "-------DRGEP coefficients-------" << endl;
      cout << "--------------------------------" << endl;
   }

   for (int k=0; k<nsp_init; k++)
   {
      if (rank ==0) cout << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
      log << sort_R_AD[k][0] << "  " << listSpecies_init[sort_R_AD[k][1]]->m_Name << endl;
   }

   log.close();


   string outputSchemeName;
   int nTarg = 0;
   for (int i = 0; i < Targets.size(); i++) {
      nTarg += Targets[i] ? 1 : 0;
   }
   int Nmin = max(9, nTarg);
   //
   if (rank ==0) {
      string konsole_script;
      ofstream writer_inputs_file;
      string writer_inputs_filename = "inputs.in";
      string yaml_mech = std::regex_replace(inputs.mech,std::regex(".xml"), ".yaml");
      //
      for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>Nmin; nbSpeciesToKeep--)
      {
            outputSchemeName = " ./outputs/mechanisms/drgepSpecies";
            stringstream s_nbSpeciesToKeep; s_nbSpeciesToKeep << nbSpeciesToKeep;
            outputSchemeName.append(s_nbSpeciesToKeep.str());
            writer_inputs_file.open(writer_inputs_filename, ios::out | ios::trunc); 
            writer_inputs_file << yaml_mech << " orch " << nbSpeciesToKeep << outputSchemeName << endl;
            writer_inputs_file.close();
            //
            konsole_script = "python submech_writer.py";
            system(konsole_script.c_str());
            konsole_script = "sed -i \'/dispersion_coefficient/d\' ";
            konsole_script.append(outputSchemeName).append(".xml");
            system(konsole_script.c_str());
            konsole_script = "sed -i \'/quadrupole_polarizability/d\' ";
            konsole_script.append(outputSchemeName).append(".xml");
            system(konsole_script.c_str());
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);

   for (int nbSpeciesToKeep=nsp_init-1; nbSpeciesToKeep>Nmin; nbSpeciesToKeep--)
   {
      outputSchemeName = "./outputs/mechanisms/drgepSpecies";
      stringstream s_nbSpeciesToKeep; s_nbSpeciesToKeep << nbSpeciesToKeep;
      outputSchemeName.append(s_nbSpeciesToKeep.str()).append(".xml");
      
      if (inputs.configuration == "MultipleInlet")
      {
         vector<vector<vector<double> > > Ym_Trajectories (nbInlets, vector<vector<double> > (nbIterations+1, vector<double> (nbSpeciesToKeep, 0.0)));
         vector<vector<double> > T_Trajectories (nbInlets, vector<double> (nbIterations+1, 0.0));
         computeMultipleInlet *c = new computeMultipleInlet();
         c->getMultipleInlet(inputs,
               outputSchemeName,
               listInlets, 
               Targets,
               false,
               "Species", 
               R_AD_Trajectories, 
               max_j_on_Target, 
               Ym_Trajectories, 
               Production_Trajectories_ref, 
               Consumption_Trajectories_ref, 
               T_Trajectories, 
               SpeciesIntoReactants);
      } //end if configuration == "MultipleInlet"

      if (inputs.configuration == "PremixedFlames")
      {
         string outputName = "./outputs/Premixed/Reduced_DRGEP_Species";
         outputName.append(s_nbSpeciesToKeep.str()).append("_").append(s_nbFlame.str());

         computePremixedFlames(outputSchemeName, 
                               outputName, 
                               inputs.mech_desc, 
                               listFlames, 
                               Targets, 
                               "ComputeTrajectories", 
                               SpeciesIntoReactants, 
                               R_AD_Premixed, 
                               max_j_on_Target, 
                               max_jf_on_Target, 
                               max_jr_on_Target, 
                               QSS_Criteria, 
                               false);
      }
  }// end for nbToKeep
}





