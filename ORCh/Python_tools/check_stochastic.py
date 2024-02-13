import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cantera as ct
import os, time
# from copy import deepcopy
#
def get_filename(file_type,case,nsp,nreac):
	if (file_type == "Stochastic"):
		if (case == "drgep_species"):
			name = "outputs/Particles"+str(nsp)+".dat"
		else:
			name = "outputs/Particles"+str(nsp)+"_"+str(nreac)+".dat"
	else:
		if (case == "drgep_species"):
			name = "outputs/Trajectories"+str(nsp)+".dat"
		else:
			name = "outputs/Trajectories"+str(nsp)+"_"+str(nreac)+".dat"
	#
	return name
#
class Particle: # functions related to local particle composition/conditions
	def __init__(self, nsp, gas, fuel):
		self.T = np.float64(0.0)
		self.Z = np.float64(-1.0)
		self.Yk = np.zeros(nsp,np.double)
		#
		self.gas = gas
		#
		self.beta_ox = np.float64(0.0)
		self.rbeta = np.float64(0.0) # 1.0/(beta_fuel - beta_ox)
		self.Zst = np.float64(0.0)
		self.Zphi = np.empty((0,1),np.float64)
		self.initZ(nsp, fuel)
	#
	def initZ(self,nsp, fuel):
		#
		# Compute beta_fuel and get surrogate mean composition
		#
		gas = self.gas
		gas.TPX = 300.0, 1.0e5, fuel
		nfuel = fuel.count(":")
		self.wc = gas.atomic_weight("C")
		self.wh = gas.atomic_weight("H")
		self.wo = gas.atomic_weight("O")
		#
		self.cc = 0.0
		self.ch = 0.0
		self.co = 0.0
		X = gas.X
		for k in range(nsp):
			Xloc = X[k]
			spec = gas.species_name(k)
			self.cc += gas.n_atoms(spec,"C") * Xloc
			self.ch += gas.n_atoms(spec,"H") * Xloc
			self.co += gas.n_atoms(spec,"O") * Xloc
		self.cc /= nfuel
		self.ch /= nfuel 
		self.co /= nfuel 
		#print("fuel mean composition (c/h/o) : "+str(self.cc)+"/"+str(self.ch)+"/"+str(self.co))
		#
		zc = gas.elemental_mass_fraction("C")
		zh = gas.elemental_mass_fraction("H")
		zo = gas.elemental_mass_fraction("O")
		beta_fuel = np.float64(zc/(self.cc*self.wc) + zh/(self.ch*self.wh) - 2.0*zo/((self.co+2.0*(self.cc+0.25*self.ch-0.5*self.co))*self.wo))
		#
		# Compute beta_ox and rbeta
		#
		gas.TPX = 300.0, 1.0e5, "O2:0.21,N2:0.79"
		zc = gas.elemental_mass_fraction("C")
		zh = gas.elemental_mass_fraction("H")
		zo = gas.elemental_mass_fraction("O")
		#print("O2 mass fraction in air : "+str(zo))
		self.beta_ox = np.float64(zc/(self.cc*self.wc) + zh/(self.ch*self.wh) - 2.0*zo/((self.co+2.0*(self.cc+0.25*self.ch-0.5*self.co))*self.wo))
		self.rbeta = np.float64(1.0) / (beta_fuel - self.beta_ox)
		self.Zst = -self.beta_ox * self.rbeta
		#
		# Postpro to indicate phi
		#
		for phi in [0.5,2.0,3.0]:
			strfuel = fuel.split(":")
			strfuel = strfuel[0]
			fuel_ox_ratio = gas.n_atoms(strfuel,'C') + 0.25*gas.n_atoms(strfuel,'H') - 0.5*gas.n_atoms(strfuel,'O')
			dil = fuel_ox_ratio*(0.79/0.21)
			ox = fuel_ox_ratio
			strfuel = strfuel + ":"+str(phi)
			reactants = strfuel+', O2:'+str(ox)+', N2:'+str(dil)  # premixed gas composition
			#
			gas.TPX = 300.0, 1.0e5, reactants
			zc = gas.elemental_mass_fraction("C")
			zh = gas.elemental_mass_fraction("H")
			zo = gas.elemental_mass_fraction("O")
			beta = np.float64(zc/(self.cc*self.wc) + zh/(self.ch*self.wh) - 2.0*zo/((self.co+2.0*(self.cc+0.25*self.ch-0.5*self.co))*self.wo))
			Zloc = (beta - self.beta_ox)*self.rbeta
			self.Zphi = np.append(self.Zphi,Zloc)
	#	
	def calcZ(self, Yk):
		gas = self.gas
		gas.TPY = Yk[-1], 1.0e5, Yk[:-1]
		zc = gas.elemental_mass_fraction("C")
		zh = gas.elemental_mass_fraction("H")
		zo = gas.elemental_mass_fraction("O")
		beta = np.float64(zc/(self.cc*self.wc) + zh/(self.ch*self.wh) - 2.0*zo/((self.co+2.0*(self.cc+0.25*self.ch-0.5*self.co))*self.wo))
		self.Z = (beta - self.beta_ox)*self.rbeta

		return self.Z
#
start = time.time()
# -----------------------------------------------------------------------------
# Inputs start here
# -----------------------------------------------------------------------------
mech = "mechanisms/gri12.xml" # user input - initial mech for this step
ppro_type = "Stochastic"
case = "drgep_species" # user input - drgep species / drgep reaction
fuel = "CH4:1.0" # user input - could be read from input.ini
eps = np.float64(1.0e-6)
species_to_plot = ["CO2","H2O","CO"]
# PAH_to_plot = []
# Tcut = 1400.0
# Zcut = 0.09
# -------------------------------------------------------------------------------
# Inputs end here
# -------------------------------------------------------------------------------
gas = ct.Solution(mech)
specs = gas.species_names
n2plot = len(species_to_plot)
# nPAH2plot = len(PAH_to_plot)
#
nsp = gas.n_species
nreac = gas.n_reactions
# Open file
filename = get_filename(ppro_type,case,nsp,nreac)
df = pd.read_csv(filename, delim_whitespace=True)
df = df.where(df > 1.0e-20, 0.0)
#
# Find nparts
#
nparts = df.loc[:,"Particle_number"].max()+1
ntimesteps = df.loc[:,"#Time"].nunique()
dt = df[df["#Time"] > 0.0].loc[:,"#Time"].min()
#
# Read reference case
#
Yref = np.zeros((ntimesteps,nparts,n2plot))
# YPAHref = np.zeros((ntimesteps,nparts,nPAH2plot))
Tref = np.zeros((ntimesteps,nparts))
Zref = np.zeros((ntimesteps,nparts))
Zlist = np.empty((0,1),np.float64)
Ymax =  np.empty((0,n2plot),np.float64)
# YPAHmax = np.empty((0,nPAH2plot),np.float64)
Zind = np.zeros((ntimesteps,nparts),dtype = np.int32)
local_part = Particle(nsp,gas,fuel)
col_name_array = gas.species_names
col_name_array.append("Temperature")
for i in range(ntimesteps):
	local_df = df.loc[i*nparts:(i+1)*nparts-1,:]
	#
	# Read particles composition
	# Save reference values to compute error due to reduction
	#
	# Tref[i,:] = deepcopy(local_df["Temperature"].values)
	Tref[i,:] = local_df["Temperature"].values
	Yref[i,:,:] = local_df[species_to_plot].values
	# YPAHref[i,:,:] = local_df[PAH_to_plot].values
	Zref[i,:] = np.apply_along_axis(local_part.calcZ, 1, local_df[col_name_array].values)
	for j in range(nparts):
		#
		# Fill Zlist, Zind, Ymax and YPAHmax
		#
		if (Zref[i,j] not in Zlist):
			Zlist = np.append(Zlist,Zref[i,j])
			Ymax = np.concatenate((Ymax,np.zeros((1,n2plot))))
			# YPAHmax = np.concatenate((YPAHmax,np.zeros((1,nPAH2plot))))
		locZind = np.where(Zlist == Zref[i,j])
		Zind[i,j] = locZind[0][0]
		# Ymax --> maximum value of each species for a given mixture fraction
		Ymax[locZind,:] = np.maximum(Ymax[locZind,:],Yref[i,j,:])
		# YPAHmax[locZind,:] = np.maximum(YPAHmax[locZind,:],YPAHref[i,j,:])
#
# Loop on mechanisms versions (species or reactions)
#
mech_range = nsp if (case == "drgep_species") else nreac
Terror = np.empty((0,1),np.float64)
Tmean = np.empty((0,1),np.float64)
T99 = np.empty((0,1),np.float64)
Zerror = np.empty((0,1),np.float64)
Yerror = np.empty((0,n2plot),np.float64)
# YPAHerror = np.empty((0,nPAH2plot),np.float64)
#
tmp_error = np.zeros(n2plot)
# tmp_PAHerror = np.zeros(nPAH2plot)
#
nred_list = np.empty((0,1),np.int32)
for ired in range(mech_range):
	nred = mech_range - ired - 1
	errY = np.zeros((ntimesteps,nparts,n2plot))
	# errYPAH = np.zeros((ntimesteps,nparts,nPAH2plot))
	errT = np.zeros((ntimesteps,nparts))
	errZ = np.zeros((ntimesteps,nparts))
	T = np.zeros((ntimesteps,nparts))
	Z = np.zeros((ntimesteps,nparts))
	Y = np.zeros((ntimesteps,nparts,n2plot))
	# YPAH = np.zeros((ntimesteps,nparts,nPAH2plot))
	#
	if (case == "drgep_species"):
		filename = get_filename(ppro_type,case,nred,nreac)
		nsp = nred
		mech = "outputs/mechanisms/drgepSpecies"+str(nred)+".xml"
	else:
		filename = get_filename(ppro_type,case,nsp,nred)
	#
	if (os.path.exists(filename)):
		print(f"reading {filename:s}")
		gas = ct.Solution(mech)
		col_name_array = gas.species_names
		col_name_array.append("Temperature")
		#
		if (set(species_to_plot).issubset(col_name_array)): # and (set(PAH_to_plot).issubset(col_name_array)):			
			df = pd.read_csv(filename, delim_whitespace=True)
			isfull = df.loc[:,"#Time"].nunique()
			if (isfull==ntimesteps):
				local_part = Particle(nsp,gas,fuel)
				nred_list = np.append(nred_list, np.int32(nred))
				#
				for i in range(ntimesteps):
					local_df = df.loc[i*nparts:(i+1)*nparts-1,:]
					#
					# Read trajectories particles composition
					#
					T[i,:] = local_df["Temperature"].values
					Y[i,:,:] = local_df[species_to_plot].values
					# YPAH[i,:,:] = local_df[PAH_to_plot].values
					Z[i,:] = np.apply_along_axis(local_part.calcZ, 1, local_df[col_name_array].values)
					#
					# T --> simple relative error
					errT[i,:] = np.abs(T[i,:]-Tref[i,:]) / Tref[i,:]
					# Y --> relative error if ref values superior than threshold (eps) and current value superior than 1% of max at this mixture fraction (Z)
					errY[i,:,:] = np.divide(np.abs(Y[i,:,:]-Yref[i,:,:]), Yref[i,:,:], out=np.zeros_like(Y[i,:,:]), where=((Yref[i,:,:] > eps) & (Y[i,:,:] > Ymax[Zind[i,:],:]*1.0e-2)))
					# PAH --> to be done (like Y + Z/T constraint)
					# errYPAH[i,:,:] = np.divide(np.abs(YPAH[i,:,:]-YPAHref[i,:,:]), YPAHref[i,:,:], out=np.zeros_like(YPAH[i,:,:]), where=(Z[i,:] > Zcut) & (T[i,:] > Tcut) & (YPAHref[i,:,:] > eps) & (YPAH[i,:,:] > YPAHmax[Zind[i,:],:]*1.0e-2))
					# Z --> simple relative error (Z can be null), only sanity check (Z error should be close to 0 at any time)
					errZ[i,:] = np.divide(np.abs(Z[i,:]-Zref[i,:]), Zref[i,:], np.zeros_like(Zref[i,:]), where = Zref[i,:] > 0.0)
				#
				Terror = np.append(Terror, np.max(errT))
				Tmean = np.append(Tmean, errT.mean())
				T99 = np.append(T99, np.percentile(errT,99))
				Zerror = np.append(Zerror, np.max(errZ))
				#
				for k in range(n2plot):
					tmp_error[k] = np.max(errY[:,:,k])
				# for k in range(nPAH2plot):
				# 	tmp_PAHerror[k] = np.max(errYPAH[:,:,k])
				Yerror = np.concatenate((Yerror, [tmp_error]),axis=0)
				# YPAHerror = np.concatenate((YPAHerror, [tmp_PAHerror]),axis=0)
			else:
				if (case == "drgep_species"):
					mech = "outputs/mechanisms/drgepSpecies"+str(nred+1)+".xml"
					gas = ct.Solution(mech)
				break	
		else:
			if (case == "drgep_species"):
				mech = "outputs/mechanisms/drgepSpecies"+str(nred+1)+".xml"
				gas = ct.Solution(mech)
			break
	else:
		break
#
plt.plot(nred_list,Terror,marker='o', color = 'r')
plt.grid(True)
plt.savefig("T_max.png",format="png")
plt.clf()
#
plt.plot(nred_list,Tmean,marker='o', color = 'r')
plt.grid(True)
plt.savefig("T_mean.png",format="png")
plt.clf()
#
plt.plot(nred_list,T99,marker='o', color = 'r')
plt.grid(True)
plt.savefig("T_99.png",format="png")
plt.clf()
#
plt.plot(nred_list,Zerror,marker='o', color = 'r')
plt.grid(True)
plt.savefig("Z_max.png",format="png")
plt.clf()
#
ind_list = []
for spec in species_to_plot:
	ind_list.append(gas.species_index(spec))
#
# PAHind_list = []
# for spec in PAH_to_plot:
# 	PAHind_list.append(gas.species_index(spec))
#
for k in range(n2plot):
	specname = gas.species_name(ind_list[k])
	plt.plot(nred_list,Yerror[:,k],marker='o', color = 'r')
	plt.grid(True)
	plt.savefig(specname+"_max.png",format="png")
	plt.clf()

# for k in range(nPAH2plot):
# 	specname = gas.species_name(PAHind_list[k])
# 	plt.plot(nred_list,YPAHerror[:,k],marker='o', color = 'r')
# 	plt.grid(True)
# 	plt.savefig(specname+"_max.png",format="png")
# 	plt.clf()

end = time.time()
print(f"duration : {end-start:f}")