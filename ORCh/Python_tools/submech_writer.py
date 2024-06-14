import cantera as ct
import pandas as pd
import numpy as np
import sys, os

ct.suppress_thermo_warnings()

def WriteReducedMech(full_mech, reduced_mech, to_keep):
    import yaml2ck_orch, ck2cti_orch, ctml_writer_orch
    all_species = ct.Species.list_from_file(full_mech)
    species = []

    # Filter species
    for S in all_species:
        if S.name in to_keep:
            species.append(S)

    # Filter reactions, keeping only those that only involve the selected species
    ref_phase = ct.Solution(thermo='ideal-gas', kinetics='gas', species=all_species)
    all_reactions = ct.Reaction.list_from_file(full_mech, ref_phase)
    reactions = []

    for R in all_reactions:
        if not all(reactant in to_keep for reactant in R.reactants):
            continue

        if not all(product in to_keep for product in R.products):
            continue
        
        if R.third_body != None:
            tbs = list(R.third_body.efficiencies.keys())
            if not all (spec in to_keep for spec in tbs):
                dummy = {k:v for k,v in R.third_body.efficiencies.items() if k in to_keep}
                R.third_body.efficiencies = dummy
                #
                removed_spec = list(set(tbs) - set(list(dummy.keys())))
                if any(sname in R.equation for sname in removed_spec):
                    print(f"removing {R.equation:s} (Troe without collider)")
                    continue
        
        reactions.append(R)

    gas = ct.Solution(name="reduced_mech",
                    thermo="ideal-gas", kinetics="gas",
                    transport_model="mixture-averaged",
                    species=species, reactions=reactions)

    # Save the resulting mechanism for later use
    gas.update_user_header({"description": "reduced mech"})
    gas.write_yaml("tmp_mech.yaml", header=True)
    yaml2ck_orch.write_ck("tmp_mech.yaml", "mech.dat", "therm.dat", "tran.dat")
    ck2cti_orch.call_me("mech.dat", "therm.dat", "tran.dat", f"{reduced_mech:s}.cti")
    ctml_writer_orch.call_me(f"{reduced_mech:s}.cti")
    os.remove(f"{reduced_mech:s}.cti")
                    

def GetSpeciesToKeep(filename, threshold):
    data = pd.read_csv(filename, sep ='\s+')
    data = data.drop(data[data["Coef"]<threshold].index)
    return data

def GetOrderedSpeciesList(filename, nred):
    data = pd.read_csv(filename, sep ='\s+',skiprows=4, header=None,names=["Coef","Species"])
    #
    return data.tail(nred)

def main():
    f = open("inputs.in","r")
    lines = f.readlines()
    line = lines[0].split()
    basemechfile = line[0]
    if (line[1] == "orch"):
        nred = int(line[2])
        tokeep = GetOrderedSpeciesList("DRGEP_Species.log", nred)
        WriteReducedMech(basemechfile, f"{line[3]:s}", tokeep["Species"].to_list())
    else:
        target = float(line[2])
        tokeep = GetSpeciesToKeep("DRGEP_Species.log", target)
        WriteReducedMech(basemechfile, f"{line[3]:s}", tokeep["Species"].to_list())

if __name__ == "__main__":
    main()