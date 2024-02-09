#include "computeMultipleInlet.h"
#include <fstream>	//Huu-Tri: "read from file" library
#include <iostream>	//Huu-Tri
#include <sys/stat.h> 	//Huu-Tri: Check file exist - stat function - 20 Nov 2019

#include <numeric>
#include <iomanip>

#include <fstream>// Huu-Tri NGUYEN 20210212 - To read from text file
#include <iostream>
#include <sstream>
#include <string>

/* EMST model */
extern "C" void emst_(int* mode_emst,int* np_emst,int* nc_emst,double* f_emst,double* state_emst,double* wt_emst,double* omdt_emst,double* fscale_emst,double* cvars_emst,int* info_emst); // EMST mixing model - Edited by Kaidi - Added by Huu-Tri Nguyen 10.12.2019

//---computeMultipleInlet---
computeMultipleInlet::computeMultipleInlet() //Constructeur
{}

void computeMultipleInlet::getMultipleInlet(
   ORChInputs inputs,
   string mech,
   vector<MultipleInlet*> listInlets,
   vector<bool> Targets,
   bool new_mixing,
   string step,
   vector<vector<vector<double> > > &R_AD_Trajectories,
   vector<vector<double> > &max_j_on_Target,
   vector<vector<vector<double> > > &Ym_Trajectories_store,
   vector<vector<vector<double> > > &Production_Trajectories_ref,
   vector<vector<vector<double> > > &Consumption_Trajectories_ref,
   vector<vector<double> > &T_Trajectories_store,
   vector<bool> &SpeciesIntoReactants)
{
    mixture  = new IdealGasMix(mech,inputs.mech_desc);
    int nsp = mixture->nSpecies();
    int nbInlets = listInlets.size();
	int nreac = mixture->nReactions();
	stringstream s_nsp;
	s_nsp << nsp;

    getMixedGasesComposition(listInlets, step);

    if (step == "DRGEP_Species" || step == "DRGEP_Reactions") local_drgep = new drgep(mixture);

    double t = 0.0; //Initial computational time
    int nTot = 0; //Total number of particles

    //Treat parallel stuff
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<int> nbParticles (nbInlets, 0); //Number of particles per inlet
    double Particle_flowRate = 1.0 / inputs.rMassFlowRate;
    for (int n=0; n<nbInlets; n++)
    {
        nbParticles[n] = listInlets[n]->m_flowRate*inputs.rMassFlowRate;
        if (rank == 0)
        {
            cout << "Nb particles  " << n << "  " << nbParticles[n] << endl;
        }
        nTot += nbParticles[n];
    }
    
    if (step != "Optimisation") {
		//
		// DAK : initialize Ifi_rank and Ila_rank (particles to be computed by this proc)
		//
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		
		int Ifi = 0;
		int Ila = 0;
        int proc_offset = 0;
		RecvCounts = new int[nproc];
		Disp = new int[nproc];
		for (int r=0; r<nproc; r++) Disp[r] = 0;
		
        if (inputs.traj_rank0) {
            proc_offset = 1;
            nproc -= 1;
            RecvCounts[0] = 0;
            //
            if (rank == 0) {
                Ila_rank = 0;
                Ifi_rank = 0;
                nb_var_loc = 0;
            }
        }

		for (int r=0; r<nproc; r++)
		{
			if (r < nTot%nproc)
			{
				Ifi = (nTot/nproc)*r + r;
				Ila = (nTot/nproc)*r + (nTot/nproc) + r + 1;
			}
			else
			{
				Ifi = (nTot/nproc)*r + nTot%nproc;
				Ila = (nTot/nproc)*r + (nTot/nproc) + nTot%nproc;
			}
			
			RecvCounts[r+proc_offset] = (Ila-Ifi)*(nsp+1); //for Yks and T
			for (int rb=r+proc_offset; rb<nproc-1+proc_offset; rb++) Disp[rb+1] += (Ila-Ifi)*(nsp+1);
			
			if (rank == r+proc_offset)
			{
				Ifi_rank = Ifi;
				Ila_rank = Ila;
				nb_var_loc = (Ila-Ifi)*(nsp+1);
			}
		}
        //
        if (inputs.drgepTraj) {
            R_AD_Trajectories.resize(nbInlets, vector<vector<double> > (nsp, vector<double>(nsp, 0.0)));
        } else {
            R_AD_Trajectories.resize(Ila_rank-Ifi_rank, vector<vector<double> > (nsp, vector<double>(nsp, 0.0)));
        }
	}
    if(rank==0) cout << "nTot = " << nTot << endl;
    // int AirPart = ceil(nbParticles[0]*0.4); // DAK : should not be hardcoded --> to add in conditions.cpp

    //-----------------------------------------------------------------------------------------------------------------------//
    //   Randomly select the particles that will be mixed or read that into the "Selection.dat" file if new_mixing == false
    //-----------------------------------------------------------------------------------------------------------------------//
    vector<vector<int> > Particle_1;
    vector<vector<int> > Particle_2;

    double delta_t = inputs.TimeStep;
    double tau_t = inputs.MixingTime;
    int nbLines =  inputs.IterationNumber;

    // Huu-Tri Nguyen - 20.01.2020
    if(rank==0) cout << "Time step = " << delta_t << " | Iterations = " << nbLines << " | Mixing time = " << tau_t << endl;

    double Pressure =  dynamic_cast <Characteristics_MultipleInlet*> (listInlets[nbInlets-1])->m_Pressure;
    double F =1.0;

    int Nmix = delta_t*nTot/tau_t;
    if (rank == 0) cout << "Nmix " << Nmix << endl;

    if (Nmix < 0) 
    {
        cout << "Problem with the delta_t and tau_t description, Nmix = " << Nmix << endl;
        getchar();
    }

    // by Kaidi@2019.11 - Huu-Tri Nguyen added - 10.12.2019
    if (Nmix*2 > nTot)
    {
        cout << "Problem with the delta_t and tau_t description, Nmix is too large!" << endl;
        getchar();
    }

    //seed the random number generator
    srand(time(NULL));

    if (new_mixing == false || step == "Optimisation_XML" || step == "Optimisation_Analytical")
    {
        ifstream Selection_read("Selection.dat");

        for (int nb=0; nb<nbLines; nb++)
        {
            Particle_1.push_back(vector<int> (Nmix));
            Particle_2.push_back(vector<int> (Nmix));

            double a;
            Selection_read >> a;
            for (int sp=0; sp<Nmix; sp++)
            {
                Selection_read >> Particle_1[nb][sp] >> Particle_2[nb][sp];
            }
        }
        Selection_read.close();
    }
    else
    {
        if (inputs.activateCurl) {
            ofstream Selection("Selection.dat");

            int select;

            for (int nb=0; nb<nbLines; nb++)
            {
                Particle_1.push_back(vector<int> (Nmix));
                Particle_2.push_back(vector<int> (Nmix));

                Selection << nb << "  ";

                for (int sp=0; sp<Nmix; sp++)
                {
                    //Select first and second particle to be mixed
                    for (int fs=0; fs<2; fs++)
                    {
                        bool particle_found = false;
                        while (!particle_found)
                        {
                            select = rand() % nTot;

                            if (sp == 0 && fs == 0) particle_found = true;
                            if (sp == 0 && fs == 1)
                            {
                                if (select != Particle_1[nb][0]) particle_found = true;
                            }
                            /* Add the corrected code from Kaidi -  Huu-Tri Nguyen 10.12.2019 */
                            bool particle_used = false;
                            for (int sp_test=0; sp_test<sp; sp_test++)
                            {
                                if (select == Particle_1[nb][sp_test] || select == Particle_2[nb][sp_test])
                                {
                                    particle_used = true;
                                }
                            }
                    
                            if (fs == 1)
                            {
                                if (select == Particle_1[nb][sp]) particle_used = true;
                            }
                    
                            if (particle_used) {
                                particle_found = false;
                            }
                            else {
                                particle_found = true;
                            }
                        }

                        if (fs == 0) {
                            Particle_1[nb][sp] = select;
                        }
                        else if (fs == 1) {
                            Particle_2[nb][sp] = select;
                        }
                    }
                    Selection << Particle_1[nb][sp] << " " << Particle_2[nb][sp] << "   ";
                }
                Selection << endl;
            }
            Selection.close();
        }
    }

    //---Trajectories store---
    vector<double> Hm_Trajectories(nbInlets, 0.0);
    vector<double> Zm_Trajectories(nbInlets, 0.0);
    vector<double> Hm_inletIni(nbInlets, 0.0); // Save initial enthalpy of each inlet - Huu-Tri Nguyen - 16.01.2020

    vector<double> Ym(nsp,0.0);
    double Hm = 0.0;
    double Zm = 0.0;
    double Tm = 0.0;
    double density = 0.0;

    //First create the Particles which will transport the gaseous and liquid phases
    vector<Particle*> listParticles;

    double Diameter_init; 
    double tau_vj;
    bool flagEvap = false;
    for (int n=0; n<nbInlets; n++)
    {      
        if (listInlets[n]->m_X_Species != "")
        {
            if (rank == 0) cout << "Set the mole fraction of inlet " << n << endl;
            mixture->setState_TPX(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_X_Species);
        }
        else if (listInlets[n]->m_Y_Species != "")
        {
            if (rank == 0)
            {
                cout << "Set the mass fraction of inlet " << n << endl;
            }
            mixture->setState_TPY(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_Y_Species);
        }

        mixture->getMassFractions(&Ym[0]);
        Hm  = mixture->enthalpy_mass();
        vector<double> Y_transfer (nsp, 0.0);
        for (int k=0; k<nsp; k++) Y_transfer[k] = Ym[k];
        density = mixture->density();

        for (int k=0; k<nsp; k++)
        {
            if (Ym[k] > 0.0)
            {
                SpeciesIntoReactants[k] = true;
            }
        }
	    
        Hm_inletIni[n] = Hm; // Save initial enthalpy of each inlet - Huu-Tri Nguyen - 16.01.2020
	    if(rank ==0) cout << "Hm initial inlet " << n << " = " << Hm_inletIni[n] << endl;

        if (n < nbInlets)  
        {
            //Composition space Lagrangian trajectories
            for (int k=0; k<nsp; k++) Ym_Trajectories_store[n][0][k] = Y_transfer[k];
 
            Hm_Trajectories[n] = Hm;
            T_Trajectories_store[n][0] = listInlets[n]->m_Temperature;
        }
	
        for (int i=0; i<nbParticles[n]; i++)
        {
            if (listInlets[n]->m_liquid)
            {
                flagEvap = true;
                double nbDroplets = Particle_flowRate/(listInlets[n]->m_density_liquid*(PI/6)*pow(listInlets[n]->m_DropletsDiameter,3.0));
                listParticles.push_back(new Particle(
                        Y_transfer, 
                        listInlets[n]->m_Temperature, 
                        Hm, 
                        1, 
                        0.0, 
                        nbDroplets, 
                        listInlets[n]->m_DropletsDiameter, 
                        0.001 /*0.01% gas, 99.9% liquid*/, 
			            Y_transfer, 
                        listInlets[n]->m_density_liquid, 
                        listInlets[n]->m_EvaporationLatentHeat));

                Diameter_init = listInlets[n]->m_DropletsDiameter;
                tau_vj = listInlets[n]->m_Tau_vj;
            }
            else
            {
                listParticles.push_back(new Particle(
                        Y_transfer, 
                        listInlets[n]->m_Temperature, 
                        Hm, 
                        0, 
                        density, 
                        0, 
                        0.0, 
                        1.0 /*100% gas, 0% liquid*/, 
                        vector<double> (nsp, 0.0), 
                        0.0, 
                        0.0));
            }
        }
    }

    //Table with the sensible enthalpy of species i at Tboiling
    vector<double> BoilingEnthalpy (nsp, 0.0);
    if (flagEvap) { // DAK : inlet hardcoded, wrong
        for (int k=0; k<nsp; k++)
        {
            double *Ymb = new double[nsp];
            for (int kbis=0; kbis<nsp; kbis++)
            {
                if (k != kbis) {
                    Ymb[kbis] = 0;
                }
                else 
                {
                    Ymb[kbis] = 1;
                }
            }
            mixture->setState_TPY(listInlets[1]->m_Temperature, listInlets[1]->m_Pressure, Ymb);
            BoilingEnthalpy[k] = mixture->enthalpy_mass();
            delete[] Ymb;
        }
    }

    double dt = inputs.TimeStep;
    double gas_mass_p1;
    double gas_mass_p2;
    double Total_gas_mass;

    vector<double> Mean_Ym(nsp, 0.0);
    double Mean_Hm = 0.0;
    double Mean_Zm = 0.0;
    double Mean_Tm = 0.0;
	
    // Data for scatterplot
    // ofstream store ("outputs/data.dat"); // Store mean values
    // ofstream store_particles ("data_particles.dat");	// Store scatterplot
    // store_particles << "#1:time  2:Particle_number  3:Zfraction  4:T(K) 5:Y_CO2 6:Y_O2 7:particle type  8:ratio  9:Zst " << endl;

    /* Add EMST model initialization - Huu-Tri Nguyen - 10.12.2019 */
    // EMST mixing model initialization -- by Kaidi@2019.12
    int mode_emst = 1;		// 1:-the state variables of all particles are initialized. No mixing is performed (check emst.f)
    int np_emst;		// np_emst: number of particles
    int nc_emst;		// nc_emst: number of compositions variables (+1 to add Enthalpy at the end of array)
    if (!inputs.activateCurl) {
		np_emst = nTot;
		nc_emst = nsp+1;
	} else {
		np_emst = 1;
		nc_emst = 1;
	}
    double *f_emst = new double[np_emst*nc_emst];	// Particle composition
    double *state_emst = new double[np_emst];		// State (or age) variable
    double *wt_emst = new double[np_emst];		// Particle weights. wt(i)>0 is the numerical weight of the i-th particle
    double *fscale_emst = new double[nc_emst];		// Scale factors for the compositions (fscale(j)>0) - see explanation in emst.f
    double cvars_emst[6] = {0,0,0,0,0,0};		// Control variables for expert users. cvars_emst[6] = {0,0,0,0,0,0} is default
    int info_emst = 0;					// = 0 for successful execution; < 0 for error in input
    int i_emst;
    double Cphi = 2; 					// Mixing model constant, see explanation in emst.f
    double omdt_emst = delta_t/tau_t*Cphi;		// Normalized mixing time 
    if(rank ==0) cout << "omdt = " << omdt_emst << endl;
    i_emst = 0;
    if (!inputs.activateCurl) {
        for (int ic=0; ic<nc_emst; ic++)	// ic: running variable for number of composition nc_emst
        {
            for (int ip=0; ip<np_emst; ip++)	// ik: running variable for number of particle np_emst
            {
                if (ic<nsp) {
                    f_emst[i_emst] = listParticles[ip]->m_Yk_gas[ic];
                }
                else
                {
                    f_emst[i_emst] = listParticles[ip]->m_H_gas;		// Ad enthalpy in the end of array to calculate
                }
                i_emst++;
            }
        }

        for (int ip=0; ip<np_emst; ip++)
        {
            state_emst[ip] = 0;
            wt_emst[ip] = 1.0;
        }

        for (int ic=0; ic<nc_emst; ic++)
        {
            if (ic<nsp) {
                fscale_emst[ic] = 0.1;	// Intialization values - recommended by emst.f
            }
            else {
                fscale_emst[ic] = 1e16;
            }
        }

        emst_(&mode_emst,&np_emst,&nc_emst,f_emst,state_emst,wt_emst,&omdt_emst,fscale_emst,cvars_emst,&info_emst);
        if (info_emst != 0)
        {
            cout << "emst initialization failed" << endl;
            getchar();
        }
    }
    /* END Add EMST model initialization - Huu-Tri Nguyen - 10.12.2019 */

    int ndil = 0;
    
    // To count EMST mixed particles
    int modifEMSTParticle = 0; // To count the number of particle modified by EMST

	//
	// DAK - initialization of postprocessing - Export all particles data
	//
	// Check if file names (part_output_name and traj_output_name) exist at the first step
	// If yes, clear the file content
	string part_output_name = "outputs/Particles";
    string traj_output_name = "outputs/Trajectories";
	if ((step == "Reactions")||(step == "DRGEP_Reactions")) {
		stringstream s_nreac;
		s_nreac << nreac;
		part_output_name.append(s_nsp.str()+"_"+s_nreac.str()+".dat");
        traj_output_name.append(s_nsp.str()+"_"+s_nreac.str()+".dat");
	} else if ((step == "Species")||(step == "DRGEP_Species")) {
		part_output_name.append(s_nsp.str()+".dat");
        traj_output_name.append(s_nsp.str()+".dat");
	} else if (step == "Optimisation") {
		stringstream s_nreac;
		s_nreac << nreac;
		stringstream s_rank;
		s_rank << rank;
		part_output_name.append(s_nsp.str()+"_"+s_nreac.str()+"_optim"+s_rank.str()+".dat");
        traj_output_name.append(s_nsp.str()+"_"+s_nreac.str()+"_optim"+s_rank.str()+".dat");
		rank = 0;
	} else {
		part_output_name.append(".dat");
        traj_output_name.append(".dat");
	}
	const int outputs_width = 16;
	ofstream full_dataParticle;
    ofstream full_dataTrajectory;
	if (rank == 0) {
        if (inputs.writeAllPart) {
            if(access(part_output_name.c_str(), F_OK) != -1)
            {
                cout << " -------------- Warning --------------" << endl;
                cout << part_output_name << " exists. Clearing file ... " << endl;
                ofstream full_dataParticle_clear(part_output_name, ios::out | ios::trunc);	//open file in trunc mode to clear the content
                full_dataParticle_clear.close(); //close the file
            }
            full_dataParticle.open(part_output_name); 
            full_dataParticle << std::setw(outputs_width) << "#Time";
            full_dataParticle << std::setw(outputs_width) << "Particle_number";
            for (int k=0; k<nsp; k++) full_dataParticle << std::setw(outputs_width) << mixture->speciesName(k);
            full_dataParticle << std::setw(outputs_width) << "Temperature" << endl;
        }
        //
        //
        //
        if (inputs.writeTraj) {
            if(access(traj_output_name.c_str(), F_OK) != -1)
            {
                cout << " -------------- Warning --------------" << endl;
                cout << traj_output_name << " exists. Clearing file ... " << endl;
                ofstream full_dataTrajectory_clear(traj_output_name, ios::out | ios::trunc);	//open file in trunc mode to clear the content
                full_dataTrajectory_clear.close(); //close the file
            }
            full_dataTrajectory.open(traj_output_name); 
            full_dataTrajectory << std::setw(outputs_width) << "#Time";
            full_dataTrajectory << std::setw(outputs_width) << "Particle_number";
            for (int k=0; k<nsp; k++) full_dataTrajectory << std::setw(outputs_width) << mixture->speciesName(k);
            full_dataTrajectory << std::setw(outputs_width) << "Temperature" << endl;
        }
	}

    /* =============================== */
    /* BEGIN BIG LOOP - Each time step */
    /* =============================== */
    for (int i=0; i<nbLines; i++)
    {
		// =========== FULL DATA PARTICLES - Huu-Tri Nguyen 19.12.2019=========== + DAK corrections
		if(rank==0) //Parallel MPI stuff to prevent print multiple lines, rank = 0 <=> First processor 
		{
			// Store Temperature & Mass fraction for each particle at each time step - Huu-Tri Nguyen 10.12.2019
            if (inputs.writeAllPart) {
                for(int p=0; p<nTot; p++)
                {
                    full_dataParticle << std::setw(outputs_width) << t;
                    full_dataParticle << std::setw(outputs_width) << p;
                    for (int k=0; k<nsp; k++) full_dataParticle << std::setw(outputs_width) << listParticles[p]->m_Yk_gas[k];
                    full_dataParticle << std::setw(outputs_width) << listParticles[p]->m_T_gas << endl;
                }
            }

            if (inputs.writeTraj) {
                for(int p=0; p<nbInlets; p++)
                {
                    full_dataTrajectory << std::setw(outputs_width) << t;
                    full_dataTrajectory << std::setw(outputs_width) << p;
                    for (int k=0; k<nsp; k++) full_dataTrajectory << std::setw(outputs_width) << Ym_Trajectories_store[p][i][k];
                    full_dataTrajectory << std::setw(outputs_width) << T_Trajectories_store[p][i] << endl;
                }
            }
		} //END if(rank==0)
		// =========== END FULL DATA PARTICLES ===========

        if (inputs.drgepTraj || inputs.writeTraj) {
            if (inputs.traj_rank0) {
                if (rank==0) {
                    //Mean values
                    for (int k=0; k<nsp; k++) Mean_Ym[k] = 0.0;
                    Mean_Hm = 0.0;
                    Mean_Tm = 0.0;
                    Total_gas_mass = 0.0;
                    for (int p=ndil; p<nTot; p++)
                    {
                        Total_gas_mass += listParticles[p]->m_P_gas_liquid; // m_P_gas_liquid for 1p = 1
                        for (int k=0; k<nsp; k++) Mean_Ym[k] += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_Yk_gas[k]);
                        Mean_Hm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);
                        Mean_Tm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_T_gas);
                    }
                    for (int k=0; k<nsp; k++) Mean_Ym[k] /= Total_gas_mass;
                    Mean_Hm /= Total_gas_mass;
                    Mean_Tm /= Total_gas_mass; // DAK: should be H/Cp - not important, only postpro

                    cout << endl;		
                    cout << " Mean_Hm at ite " << i << " = " << Mean_Hm <<endl;
                    cout << " Mean_Tm at ite " << i << " = " << Mean_Tm << endl;
                }
            } else {
                //Mean values
                for (int k=0; k<nsp; k++) Mean_Ym[k] = 0.0;
                Mean_Hm = 0.0;
                Mean_Tm = 0.0;
                Total_gas_mass = 0.0;
                for (int p=ndil; p<nTot; p++)
                {
                    Total_gas_mass += listParticles[p]->m_P_gas_liquid; // m_P_gas_liquid for 1p = 1
                    for (int k=0; k<nsp; k++) Mean_Ym[k] += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_Yk_gas[k]);
                    Mean_Hm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_H_gas);
                    Mean_Tm += (listParticles[p]->m_P_gas_liquid*listParticles[p]->m_T_gas);
                }
                for (int k=0; k<nsp; k++) Mean_Ym[k] /= Total_gas_mass;
                Mean_Hm /= Total_gas_mass;
                Mean_Tm /= Total_gas_mass; // DAK: should be H/Cp - not important, only postpro
                if(rank==0) 
                {
                    cout << endl;		
                    cout << " Mean_Hm at ite " << i << " = " << Mean_Hm <<endl;
                    cout << " Mean_Tm at ite " << i << " = " << Mean_Tm << endl;
                }
            }
        }
	
        // ====== LAGRANGIAN TRAJECTORIES - DETERMINISTIC ====== //
        if (inputs.drgepTraj || inputs.writeTraj) {
            if (inputs.traj_rank0) {
                if (rank==0) {
                    for (int n=0; n<nbInlets; n++)
                    {
                        for (int k=0; k<nsp; k++) Ym[k] = (Ym_Trajectories_store[n][i][k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];
                        Hm = (Hm_Trajectories[n]-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;

                        if ((step == "DRGEP_Species" || step == "DRGEP_Reactions") && inputs.drgepTraj)
                        {
                            Next_Time_Step_with_drgep(Targets, Pressure, Ym, Hm, Tm, delta_t, R_AD_Trajectories[n], max_j_on_Target, step, inputs.print_all_rAB);
                        }
                        else if (step == "computeQSSCriteria")
                        {
                            Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t, Production_Trajectories_ref, Consumption_Trajectories_ref, n, i);
                        }
                        else
                        {
                            Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t);
                        }
                        cout << " **** Advance Deterministic inlet " << n << " at " << i << " iterations" << endl;

                        // Huu-Tri - After the reactor
                        for (int k=0; k<nsp; k++) Ym_Trajectories_store[n][i+1][k] = (Ym[k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];
                        Hm_Trajectories[n] = (Hm-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;
                        T_Trajectories_store[n][i+1] = Tm; // DAK : not consistent with H and Yk - not important (postpro only)
                    } // End "for" each inlet
                }
            } else {
                for (int n=0; n<nbInlets; n++)
                {
                    for (int k=0; k<nsp; k++) Ym[k] = (Ym_Trajectories_store[n][i][k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];
                    Hm = (Hm_Trajectories[n]-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;

                    if ((step == "DRGEP_Species" || step == "DRGEP_Reactions") && inputs.drgepTraj)
                    {
                        Next_Time_Step_with_drgep(Targets, Pressure, Ym, Hm, Tm, delta_t, R_AD_Trajectories[n], max_j_on_Target, step, inputs.print_all_rAB);
                    }
                    else if (step == "computeQSSCriteria")
                    {
                        Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t, Production_Trajectories_ref, Consumption_Trajectories_ref, n, i);
                    }
                    else
                    {
                        Next_Time_Step(Pressure, Ym, Hm, Tm, delta_t);
                    }
                    // Huu-Tri TEST - 20210206
                    if(rank==0)
                    {
                        cout << " **** Advance Deterministic inlet " << n << " at " << i << " iterations" << endl;
                    }

                    // Huu-Tri - After the reactor
                    for (int k=0; k<nsp; k++) Ym_Trajectories_store[n][i+1][k] = (Ym[k]-Mean_Ym[k])*exp(-(delta_t/(2*tau_t)))+Mean_Ym[k];
                    Hm_Trajectories[n] = (Hm-Mean_Hm)*exp(-(delta_t/(2*tau_t)))+Mean_Hm;
                    T_Trajectories_store[n][i+1] = Tm; // DAK : not consistent with H and Yk - not important (postpro only)
                } // End "for" each inlet
            }
        }
        // ====== END LAGRANGIAN TRAJECTORIES - DETERMINISTIC ====== //
		//
        //---Particles mixing---
        //
		// DAK : probably hardcoded to add air particles after combustion
        //Add Nmix variation if dilution
        // int Nmix2;
        // if (i < nbLines/2) { // DAK : hardcoded ? probably belongs to conditions.cpp
        //     Nmix2 = delta_t*(nTot-AirPart+1)/tau_t;
        // }
        // else 
        // {
        //     float c = i;
        //     float d = nbLines;
        //     double b = abs(AirPart*(1 - ( 2*( c - d/2 )/d)));
        //     int a = b;
        //     Nmix2 = delta_t*(nTot - a)/tau_t;
        // }

        /* ======== MIXING MODELS ======== */
        /* 2 options: Curl model or EMST model */
	    /* CURL Mixing closure phi_p1 = phi_p2 = (phi_p1 + phi_p2)/2 */
		if(inputs.activateCurl)
		{
			for (int p=0; p<Nmix; p++)
			{
                // Huu-Tri: Only gas case,  gas_mass_p1 = gas_mass_p2 = 0.01 = Particle_flowRate
                gas_mass_p1 = listParticles[Particle_1[i][p]]->m_P_gas_liquid*Particle_flowRate;
                gas_mass_p2 = listParticles[Particle_2[i][p]]->m_P_gas_liquid*Particle_flowRate; 

                if (gas_mass_p1+gas_mass_p2 > 0.0)
                {
                    for (int k=0; k<nsp; k++)
                    {
                        listParticles[Particle_1[i][p]]->m_Yk_gas[k] = F*(gas_mass_p1*listParticles[Particle_1[i][p]]->m_Yk_gas[k]+gas_mass_p2*listParticles[Particle_2[i][p]]->m_Yk_gas[k])/(gas_mass_p1+gas_mass_p2);
                        listParticles[Particle_2[i][p]]->m_Yk_gas[k] = listParticles[Particle_1[i][p]]->m_Yk_gas[k];
                    }
                    listParticles[Particle_1[i][p]]->m_H_gas = (gas_mass_p1*listParticles[Particle_1[i][p]]->m_H_gas+gas_mass_p2*listParticles[Particle_2[i][p]]->m_H_gas)/(gas_mass_p1+gas_mass_p2);//Original
                                
                    //  Huu-Tri NGUYEN - Add a modified enthaply for heat loss
                    //                listParticles[Particle_1[i][p]]->m_H_gas = (gas_mass_p1*listParticles[Particle_1[i][p]]->m_H_gas+gas_mass_p2*listParticles[Particle_2[i][p]]->m_H_gas)/(gas_mass_p1+gas_mass_p2) + varEnthalpyCFD;// Huu-Tri Stochastic heat loss 14 Nov 2019 

                    listParticles[Particle_2[i][p]]->m_H_gas = listParticles[Particle_1[i][p]]->m_H_gas;
                }
			}
		} /* END Commented Curl model to add EMST model - Huu-Tri Nguyen -10.12.2019 */
		else // EMST model
		{
			if(rank ==0) cout << "EMST mixing model" << endl;
			/* Add EMST model - Huu-Tri Nguyen -10.12.2019 */
     		//EMST mixing model -- by Kaidi@2019.12
      		double run_mixing = true;

     		if (run_mixing==true)
      		{
         		mode_emst = 2;	// 2 -mixing is performed and the state variables are incremented.

         		i_emst = 0;
         		for (int ic=0; ic<nc_emst; ic++)
         		{
            			for (int ip=0; ip<np_emst; ip++)
            			{
               				if (ic<nsp)  
						f_emst[i_emst] = listParticles[ip]->m_Yk_gas[ic];
               				else         
						f_emst[i_emst] = listParticles[ip]->m_H_gas;
               			i_emst++;
            			}
         		}
    
				// Mixing with EMST model
				emst_(&mode_emst,&np_emst,&nc_emst,f_emst,state_emst,wt_emst,&omdt_emst,fscale_emst,cvars_emst,&info_emst);
				
				int EMST_mixedParticles =0;			
				for (int ip=0; ip<np_emst; ip++) //HuuTri@20201029: Print number of mixed particles
				{
					if(state_emst[ip] >0)
					{
						EMST_mixedParticles++; // Particle mixed state_emst>0, non-mixed state_emst<0
					}
					// wt_emst[ip];
				}
				if(rank==0) cout << "Mixed particles = " << EMST_mixedParticles << endl;
				if (info_emst != 0)
         		{
            			cout << "emst failed" << endl;
            			// getchar(); // DAK : exit() or MPI Abort ?
         		}

         		i_emst = 0;
         		for (int ic=0; ic<nc_emst; ic++)
         		{
					for (int ip=0; ip<np_emst; ip++)
					{
						if (ic<nsp)
						{
							listParticles[ip]->m_Yk_gas[ic] = f_emst[i_emst];
						}
						else
						{
							listParticles[ip]->m_H_gas = f_emst[i_emst];
						}
						i_emst++;
					}
         		}
			}
		} // End else EMST model
      	// End of mixing model (CURL or EMST)

		//T_o = 353;
		//if (listParticles[ip]->m_T_gas > T_o) listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas * (((rand() / double(RAND_MAX))*2.0-1.0)*1.0/100.0 + 1.0);
		//if (i > nbLines/3)
		// Heatloss Camille - Commented by Huu-Tri Nguyen -10.12.2019 			
		//HT	if (true)
		//HT    {
		//HT		T_o = 353; 
		//HT    	if (listParticles[ip]->m_T_gas > T_o) 
		//HT			listParticles[ip]->m_H_gas = listParticles[ip]->m_H_gas - K*alpha_loss*(listParticles[ip]->m_T_gas - T_o)*delta_t;
		//HT    }
		// Heatloss Camille - Corrected by Huu-Tri Nguyen - 20220107
		// HuuTri@20220107: flagHeatLossStochas == true, use a sink/source term to impose the heat loss (alpha_loss) for Stochastic closure
		// Sink term  = alpha_loss*(T_gas - T_wall)
		// Heat term  = beta_heat*(T_gas - T_wall)
		// The value of alpha_loss, beta_heat also depends on the time step delta_t
		// To calculate enthalpy: It should be multiplied by delta_t: 
		// H_gasAfterLoss = H_gas -  sinkTerm*delta_t =  H_ini - alpha_loss*(T_gas - T_wall)*delta_t
		// H_gasAfterHeated = H_gas -  sourceTerm*delta_t =  H_ini - beta_heat*(T_gas - T_wall)*delta_t
		bool flagHeatLossStochas = false; // DAK --> should be in conditions.cpp
        double T_wall = 1600.0;
        double alpha_loss = 1.5e+05;	//1.5e+06;
        double beta_heat = 7.5e+05;
		if (flagHeatLossStochas)
		{
			for (int p=ndil; p<nTot; p++)
			{
                // HuuTri@20220107: alpha_loss and T_wall should be declared above (find "flagHeatLossStochas" keyword)
                // HuuTri@20220107: If Tgas > Twall, the heat will be lost >> alpha_loss
                if (listParticles[p]->m_T_gas > T_wall)
                {
                    listParticles[p]->m_H_gas = listParticles[p]->m_H_gas - alpha_loss*(listParticles[p]->m_T_gas - T_wall)*delta_t;	
                }
                else // When Tgas < Twall, the gas is heated by walls
                {
                    listParticles[p]->m_H_gas = listParticles[p]->m_H_gas - beta_heat*(listParticles[p]->m_T_gas - T_wall)*delta_t;
                }
            }
		} // end if (flagHeatLossStochas)

	    //Evaporation - DAK : boiling enthalpy probably wrong ...
		if (flagEvap) {
			for (int p=ndil; p<nTot; p++)
			{
				if (listParticles[p]->m_P_gas_liquid < 1.0)
				{
					double Diameter;
					if ((t+dt)/tau_vj < 1) {
						Diameter = Diameter_init*pow((1-((t+dt)/tau_vj)),0.5);
					}
					else {
						Diameter = 0.0;
					}

					double old_Diameter = listParticles[p]->m_droplets_diameter;
					listParticles[p]->m_droplets_diameter = Diameter;

					double gas_mass = listParticles[p]->m_P_gas_liquid*Particle_flowRate;
					double liquid_mass = (1-listParticles[p]->m_P_gas_liquid)*Particle_flowRate;
					double ReleasedVapor = listParticles[p]->m_N_droplets*listParticles[p]->m_density_liquid*(PI/6)*(pow(old_Diameter,3)-pow(Diameter,3));
					listParticles[p]->m_P_gas_liquid = (gas_mass+ReleasedVapor)/Particle_flowRate;

					for (int k=0; k<nsp; k++)
					{
						double Yk_gas_init = listParticles[p]->m_Yk_gas[k];
						listParticles[p]->m_Yk_gas[k] = ((gas_mass*Yk_gas_init)+(ReleasedVapor*listParticles[p]->m_Yk_liquid[k]))/(gas_mass+ReleasedVapor);
					}

					double H_gas_init = listParticles[p]->m_H_gas;
					double ReleasedEnthalpy = 0.0;
					for (int k=0; k<nsp; k++) ReleasedEnthalpy += ReleasedVapor*listParticles[p]->m_Yk_liquid[k]*(BoilingEnthalpy[k]-listParticles[p]->m_EvaporationLatentHeat);
					listParticles[p]->m_H_gas = ((gas_mass*H_gas_init)+(ReleasedEnthalpy))/(gas_mass+ReleasedVapor);
				}
			}
		}

		/* ==================================== */
		/* ============ STOCHASTIC ============ */
		/* ==================================== */
		if (step == "Optimisation") {
			Reacting(listParticles, nsp, dt, Pressure);
		}   
		else
		{
            if ((step == "DRGEP_Species" || step == "DRGEP_Reactions") && !inputs.drgepTraj) {
                ReactingParallelDRGEP(Targets, listParticles, nsp, dt, Pressure, R_AD_Trajectories, max_j_on_Target, inputs.step);
            } else {
			    ReactingParallel(listParticles, nsp, dt, Pressure);
            }
		}	
		/* ==================================== */
	    
		// Store time step	
		t = t + dt;
   	} // END of big loop i

	if (!inputs.drgepTraj) {
        if (step == "DRGEP_Species") {
            int count = nsp*nsp;
            //
            double send_R_AD[count];
            double recv_R_AD[count];
            
            for (int i = 0; i < nsp; i++) {
                for (int j = 0; j < nsp; j++) {
                    send_R_AD[i*nsp+j] = 0.0;
                    for (int p = 0; p < Ila_rank-Ifi_rank; p++) {
                        send_R_AD[i*nsp+j] = (send_R_AD[i*nsp+j] < R_AD_Trajectories[p][i][j])?R_AD_Trajectories[p][i][j]:send_R_AD[i*nsp+j];
                    }
                }
            }
            //
            MPI_Allreduce(send_R_AD, recv_R_AD, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            //
            for (int i = 0; i < nsp; i++) {
                for (int j = 0; j < nsp; j++) {
                    R_AD_Trajectories[0][i][j] = recv_R_AD[i*nsp+j];
                }
            }
        } else if (step == "DRGEP_Reactions") {
            int count = nsp*nreac;
            //
            double send_maxj[count];
            double recv_maxj[count];
            for (int i = 0; i < nsp; i++) {
                for (int j = 0; j < nreac; j++) {
                    send_maxj[i*nreac+j] = max_j_on_Target[i][j];
                }
            }
            //
            MPI_Reduce(send_maxj, recv_maxj, count, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            //
            if (rank == 0) {
                for (int i = 0; i < nsp; i++) {
                    for (int j = 0; j < nreac; j++) {
                        max_j_on_Target[i][j] = recv_maxj[i*nreac+j];
                    }
                }
            }
        }
    } else {
        if (inputs.traj_rank0) {
            if (step == "DRGEP_Species") {
                int count = nbInlets*nsp*nsp;
                //
                if (rank==0) {
                    double send_R_AD[count];
                    //
                    for (int i = 0; i < nsp; i++) {
                        for (int j = 0; j < nsp; j++) {
                            for (int p = 0; p < nbInlets; p++) {
                                send_R_AD[p*nsp*nsp+i*nsp+j] = R_AD_Trajectories[p][i][j];
                            }
                        }
                    }
                    //
                    MPI_Bcast(send_R_AD, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                } else {
                    double recv_R_AD[count];
                    //
                    MPI_Bcast(recv_R_AD, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    //
                    for (int i = 0; i < nsp; i++) {
                        for (int j = 0; j < nsp; j++) {
                            for (int p = 0; p < nbInlets; p++) {
                                R_AD_Trajectories[p][i][j] = recv_R_AD[p*nsp*nsp+i*nsp+j];
                            }
                        }
                    }
                }
            }
        }
    }

    for (int p = 0; p < nTot; p++) {
		delete listParticles[p];
	}

	delete[] f_emst;
	delete[] state_emst;
	delete[] wt_emst;
	delete[] fscale_emst;
	if (step != "Optimisation") {
		delete[] RecvCounts;
		delete[] Disp;
		delete mixture;
    }

    if (step == "DRGEP_Species" || step == "DRGEP_Reactions") delete local_drgep;

} //end of main computeMultipleInlet::getMultipleInlet

void computeMultipleInlet::Reacting(vector<Particle*> &listParticles, int nsp, double dt, double Pressure)
{
    int nTot = listParticles.size();
    vector<double> Ym(nsp,0.0);
    double Hm;
    double Zm;
    double Tm;

    for (int p=0; p<nTot; p++)
    {
        for (int k=0; k<nsp; k++) Ym[k] = listParticles[p]->m_Yk_gas[k];

        Hm = listParticles[p]->m_H_gas;

        Next_Time_Step(Pressure, Ym, Hm, Tm, dt);

        for (int k=0; k<nsp; k++) listParticles[p]->m_Yk_gas[k] = Ym[k];

        listParticles[p]->m_T_gas = Tm;
    }
}

void computeMultipleInlet::ReactingParallel(vector<Particle*> &listParticles, int nsp, double dt, double Pressure)
{
    int nTot = listParticles.size();
    vector<double> Ym(nsp,0.0);
    double Hm;
    double Tm;

    for (int p=Ifi_rank; p<Ila_rank; p++)
    {
        for (int k=0; k<nsp; k++) Ym[k] = listParticles[p]->m_Yk_gas[k];

        Hm = listParticles[p]->m_H_gas;

        Next_Time_Step(Pressure, Ym, Hm, Tm, dt);

        for (int k=0; k<nsp; k++) listParticles[p]->m_Yk_gas[k] = Ym[k];

        listParticles[p]->m_T_gas = Tm;
    }

    double Data_Proc[nb_var_loc];
    double Data_All[nTot*(nsp+1)];

    int count = 0;
    for (int p=Ifi_rank; p<Ila_rank; p++)
    {
        for (int k=0; k<nsp; k++)
        {
            Data_Proc[count] = listParticles[p]->m_Yk_gas[k];
            count += 1;
        }
       
        Data_Proc[count] = listParticles[p]->m_T_gas;
        count += 1;
    }

    MPI_Allgatherv(Data_Proc, nb_var_loc, MPI_DOUBLE, Data_All, RecvCounts, Disp, MPI_DOUBLE, MPI_COMM_WORLD);

    count = 0;
    for (int p=0; p<Ifi_rank; p++)
    {
        for (int k=0; k<nsp; k++)
        {
            listParticles[p]->m_Yk_gas[k] = Data_All[count];
            count += 1;
        }

        listParticles[p]->m_T_gas = Data_All[count];
        count += 1;
    }
    //
	count += nb_var_loc;
	//
	for (int p=Ila_rank; p<nTot; p++)
	{
		for (int k=0; k<nsp; k++)
		{
			listParticles[p]->m_Yk_gas[k] = Data_All[count];
			count += 1;
		}

		listParticles[p]->m_T_gas = Data_All[count];
		count += 1;
	}
}

void computeMultipleInlet::ReactingParallelDRGEP(vector<bool> Targets, vector<Particle*> &listParticles, int nsp, double dt, double Pressure, vector<vector<vector<double> > >&R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string  step)
{
	int nTot = listParticles.size();
	vector<double> Ym(nsp,0.0);
	double Tm, Hm;
	
	for (int p=Ifi_rank; p<Ila_rank; p++)
	{
		Ym = listParticles[p]->m_Yk_gas;
		Hm = listParticles[p]->m_H_gas;
		
		Next_Time_Step_with_drgep(Targets, Pressure, Ym, Hm, Tm, dt, R_AD_Trajectories[p-Ifi_rank], max_j_on_Target, step, false);
		
		listParticles[p]->m_Yk_gas = Ym;
		listParticles[p]->m_T_gas = Tm;
	}
	
	// DAK : if too slow, check shared_memory
	double Data_Proc[nb_var_loc];
	double Data_All[nTot*(nsp+1)];
	int count = 0;
	for (int p=Ifi_rank; p<Ila_rank; p++)
	{
		for (int k=0; k<nsp; k++)
		{
			Data_Proc[count] = listParticles[p]->m_Yk_gas[k];
			count += 1;
		}
		Data_Proc[count] = listParticles[p]->m_T_gas;
		count += 1;
	}
	// cout << "rank " << rank << " particules " << nb_particles <<" time " << time(nullptr)-this_rank  << endl;
	MPI_Allgatherv(Data_Proc, nb_var_loc, MPI_DOUBLE, Data_All, RecvCounts, Disp, MPI_DOUBLE, MPI_COMM_WORLD);
	count = 0;
	for (int p=0; p<Ifi_rank; p++)
	{
		for (int k=0; k<nsp; k++)
		{
			listParticles[p]->m_Yk_gas[k] = Data_All[count];
			count += 1;
		}
		listParticles[p]->m_T_gas = Data_All[count];
		count += 1;
	}
	//
	count += nb_var_loc;
	//
	for (int p=Ila_rank; p<nTot; p++)
	{
		for (int k=0; k<nsp; k++)
		{
			listParticles[p]->m_Yk_gas[k] = Data_All[count];
			count += 1;
		}
		listParticles[p]->m_T_gas = Data_All[count];
		count += 1;
	}
    //
}

void computeMultipleInlet::getMixedGasesComposition(vector<MultipleInlet*> listInlets, string step)
{
	int rank, nproc;
    if (step != "Optimisation")
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	} else {
		rank = -1;
    }

    int nsp = mixture->nSpecies();
    int nbInlets = listInlets.size(); 

	vector<double> Compo_Yk_mixed(nsp,0.0);
    double Compo_H_mixed = 0.0;
    double Total_flowRate = 0.0;

    //Composition of the burned gases
    for (int n=0; n<nbInlets; n++)
    {
        if (listInlets[n]->m_X_Species != "")
        {
			mixture->setState_TPX(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_X_Species);
		} else if (listInlets[n]->m_Y_Species != "") {
			mixture->setState_TPY(listInlets[n]->m_Temperature, listInlets[n]->m_Pressure, listInlets[n]->m_Y_Species);
        }

		vector<double> Ym(nsp,0.0);
        double Hm = 0.0;
		mixture->getMassFractions(&Ym[0]);
		Hm = mixture->enthalpy_mass();

		for (int k=0; k<nsp; k++) Compo_Yk_mixed[k] += listInlets[n]->m_flowRate*Ym[k];
        Compo_H_mixed += listInlets[n]->m_flowRate*Hm;
        Total_flowRate += listInlets[n]->m_flowRate;
    }
   
    for (int k=0; k<nsp; k++)
    {
        Compo_Yk_mixed[k] /= Total_flowRate;
    }
    Compo_H_mixed /= Total_flowRate;
    if (rank == 0)
    {
        cout << endl <<  "Composition to enter for the equilibrium computation to get the Burned gases" << endl;
        cout << "Compo_H_mixed " << Compo_H_mixed << endl;
    }

	mixture->setMassFractions(&Compo_Yk_mixed[0]);
	mixture->setState_HP(Compo_H_mixed, listInlets[0]->m_Pressure);

	vector<double> Xm(nsp);
    double T_mixed = 0.0;
	mixture->getMoleFractions(&Xm[0]);
	T_mixed = mixture->temperature();

    if (rank == 0)
    {
        for (int k=0; k<nsp; k++)
        {
            if (Xm[k] != 0.0) cout << "X_" << mixture->speciesName(k) << ": " << Xm[k] << endl;
        }
        cout << "T_mixed " << T_mixed << endl;
    }
}

void computeMultipleInlet::Next_Time_Step_with_drgep(vector<bool> Targets, double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t, 
   vector<vector<double> > &R_AD_Trajectories, vector<vector<double> > &max_j_on_Target, string step, bool print_all_rAB)
{
	try {
		mixture->setMassFractions(&Ym[0]);
        mixture->setState_HP(Hm, P);

        ConstPressureReactor reac;
        reac.insert(*mixture);
        ReactorNet sim;
        sim.addReactor(reac);

       sim.advance(delta_t);
	} catch (...) {
		return;
	}

    Tm = mixture->temperature();
	mixture->getMassFractions(&Ym[0]);

    local_drgep->drgep_0D_species(mixture, Targets, R_AD_Trajectories, 0, 0.0, print_all_rAB);

    if (step == "DRGEP_Reactions")
    {
        int nsp = mixture->nSpecies();
        int nreac = mixture->nReactions();

        vector<vector<double> > rj_for_k (nsp, vector<double> (nreac,0.0));
        local_drgep->drgep_0D_reactions(mixture, rj_for_k);

        for (int ka=0; ka<nsp; ka++)
        {
            for (int kb=0; kb<nsp; kb++)
            {
                for (int j=0; j<nreac; j++)
                {
                if (max_j_on_Target[ka][j] < R_AD_Trajectories[ka][kb]*rj_for_k[kb][j])
                    max_j_on_Target[ka][j] = R_AD_Trajectories[ka][kb]*rj_for_k[kb][j];
                }
            }
        }
    }
}

//Next_Time_Step without drgep analysis (to use for the computations with optimisation)
void computeMultipleInlet::Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t)
{
	try {
		mixture->setMassFractions(&Ym[0]);
        mixture->setState_HP(Hm, P);

        ConstPressureReactor reac;
        reac.insert(*mixture);
        ReactorNet sim;
        sim.addReactor(reac);

        sim.advance(delta_t);
	} catch (...) {
		return;
	}

    Tm = mixture->temperature();
	mixture->getMassFractions(&Ym[0]);
}

//Next_Time_Step with QSS analysis 
void computeMultipleInlet::Next_Time_Step(double P, vector<double> &Ym, double &Hm, double &Tm, double delta_t,
                    vector<vector<vector<double> > > &Production_Trajectories_ref, vector<vector<vector<double> > > &Consumption_Trajectories_ref, int nInlet, int nLine)
{
    mixture->setMassFractions(&Ym[0]);
    mixture->setState_HP(Hm, P);

    ConstPressureReactor reac;
    reac.insert(*mixture);
    ReactorNet sim;
    sim.addReactor(reac);

    sim.advance(delta_t);

    Tm = mixture->temperature();
    mixture->getMassFractions(&Ym[0]);

    int nreac = mixture->nReactions();
    int nsp = mixture->nSpecies();

    double* fwdRates = new double[nreac];
    double* revRates = new double[nreac];
    mixture->getFwdRatesOfProgress(fwdRates);
    mixture->getRevRatesOfProgress(revRates);

    for (int k=0; k<nsp; k++)
    {
        double omega_k_prod = 0.0;
        double omega_k_cons = 0.0;
        for (int j=0; j<nreac; j++)
        {
            omega_k_prod += mixture->productStoichCoeff(k,j)*fwdRates[j]
                            +mixture->reactantStoichCoeff(k,j)*revRates[j];
            omega_k_cons += mixture->reactantStoichCoeff(k,j)*fwdRates[j]
                            +mixture->productStoichCoeff(k,j)*revRates[j];
        }
        Production_Trajectories_ref[nInlet][nLine][k] = omega_k_prod;
        Consumption_Trajectories_ref[nInlet][nLine][k] = omega_k_cons;
    }
}

computeMultipleInlet::~computeMultipleInlet() //Destructeur
{}

