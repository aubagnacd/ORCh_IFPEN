import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cantera as ct
import pandas as pd
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
		self.beta_ox = np.float64(0.0)
		self.rbeta = np.float64(0.0) # 1.0/(beta_fuel - beta_ox)
		self.Zst = np.float64(0.0)
		self.Zphi = np.empty((0,1),np.float64)
		self.initZ(nsp, gas, fuel)
	#
	def initZ(self,nsp, gas, fuel):
		#
		# Compute beta_fuel and get surrogate mean composition
		#
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
	def calcZ(self, gas):
		gas.TPY = self.T, 1.0e5, self.Yk
		zc = gas.elemental_mass_fraction("C")
		zh = gas.elemental_mass_fraction("H")
		zo = gas.elemental_mass_fraction("O")
		beta = np.float64(zc/(self.cc*self.wc) + zh/(self.ch*self.wh) - 2.0*zo/((self.co+2.0*(self.cc+0.25*self.ch-0.5*self.co))*self.wo))
		self.Z = (beta - self.beta_ox)*self.rbeta
# -------------------------------------------------------------------------
# Inputs start
# -------------------------------------------------------------------------
ppro_type = "Stochastic"
mech = 'mechanisms/gri12.xml' # user input - initial mech for this step
case = "drgep_species" # user input - drgep species / drgep reaction
fuel = "CH4:1.0" # user input - could be read from input.ini
Tmin, Tmax = 300.0, 2600.0 # temperature range to plot
# --------------------------------------------------------------------------
# Inputs end
# --------------------------------------------------------------------------
gas = ct.Solution(mech)
nsp = gas.n_species
nreac = gas.n_reactions
indPAH = gas.species_index("C2H2")
#
# Get filenames and read files
#
local_part = Particle(nsp,gas,fuel)
# Stochastic
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
Zhist = np.empty((0,1),np.float64)
Thist = np.empty((0,1),np.float64)
for i in range(ntimesteps):
	local_df = df.loc[i*nparts:(i+1)*nparts-1,:]
	#
	Z = np.empty((0,1),np.float64)
	PAH = np.empty((0,1),np.float64)
	T = local_df["Temperature"].values
	#
	# Read stochastic particles composition
	#
	for j in range(nparts):
		local_part.Yk = local_df.loc[local_df["Particle_number"] == j,gas.species_names].values
		local_part.T = T[j]
		local_part.calcZ(gas)
		#
		# Add particle to matplotlib scatterplot
		#
		Z = np.append(Z,local_part.Z)
		PAH = np.append(PAH, np.max((np.float64(1.0e-10),local_part.Yk[0][indPAH])))
	#
	# Save Z and T
	#
	Zhist = np.append(Zhist,Z)
	Thist = np.append(Thist,T)
	#
	# Write current instant scatterplot
	#
	fig, ax = plt.subplots()
	ax.scatter(Zhist,Thist, c = '0.7', s = 25)
	sc = ax.scatter(Z, T, c=PAH, s=15, norm=colors.LogNorm(vmin=1.0e-6, vmax=1.0e-1), cmap='inferno')
	ax.set_ylabel(r"$T$ $[K]$")
	ax.set_xlabel(r"$Z$ $[-]$")
	ax.vlines(local_part.Zst,300.0,3000.0)
	ax.vlines(local_part.Zphi[0],300.0,3000.0)
	ax.vlines(local_part.Zphi[1],300.0,3000.0)
	ax.vlines(local_part.Zphi[2],300.0,3000.0)
	plt.xlim(-0.01,1.01)
	plt.ylim(Tmin,Tmax)
	plt.colorbar(sc)
	fig.tight_layout()
	fig.savefig("./TZ_plot_"+str(i)+".png")
	plt.close()
	print("plotted time_steps "+str(i))
