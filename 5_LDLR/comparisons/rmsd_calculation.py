from modeller import *
from modeller.scripts import complete_pdb
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

# log.verbose()    # request verbose output
# env = Environ()
# env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
# env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# # read model file
# mdl = complete_pdb(env, 'LDLR.B99990005.pdb')

# # Assess all atoms with DOPE:
# s = Selection(mdl)
# model_dope=s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='LDLR-model5.profile',
#               normalize_profile=True, smoothing_window=15)



# #List of protein structure files (PDB format)
# protein_files = ['1n7d.pdb']  # Add your files

# #Loop through the protein structure files
# for pdb_file in protein_files:
#     # read model file
#     mdl = complete_pdb(env, pdb_file, model_segment=('FIRST:A', 'LAST:A'))

# # Assess all atoms with DOPE:
# s = Selection(mdl)
# output_file = pdb_file.replace('.pdb', '_profile.profile')

# template_dope=s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='LDLR.profile',
#               normalize_profile=True, smoothing_window=15)


# print("Dope score of Model: ",model_dope)
# print("Dope score of Template: ",template_dope)

def calculate_rmsd(template_pdb, model_pdb):
    # Load structures using BioPython's PDB parser
    parser = PDBParser(QUIET=True)
    template_structure = parser.get_structure("template", template_pdb)
    model_structure = parser.get_structure("model", model_pdb)

    # Extract atomic coordinates
    template_coords = []
    model_coords = []
    for template_model, model_model in zip(template_structure, model_structure):
        for template_chain, model_chain in zip(template_model, model_model):
            for template_res, model_res in zip(template_chain, model_chain):
                for template_atom, model_atom in zip(template_res, model_res):
                    if template_atom.get_id() == model_atom.get_id():
                        template_coords.append(template_atom.get_coord())
                        model_coords.append(model_atom.get_coord())

    # Calculate RMSD using the Kabsch algorithm
    template_coords = np.array(template_coords)
    model_coords = np.array(model_coords)

    # Center the coordinates
    template_center = np.mean(template_coords, axis=0)
    model_center = np.mean(model_coords, axis=0)
    template_coords -= template_center
    model_coords -= model_center

    # Calculate the rotation matrix using SVD
    h = np.dot(template_coords.T, model_coords)
    u, _, v_t = np.linalg.svd(h)
    rotation_matrix = np.dot(v_t.T, u.T)

    # Apply the rotation and calculate RMSD
    model_coords_rotated = np.dot(model_coords, rotation_matrix)
    rmsd = np.sqrt(np.mean(np.sum((template_coords - model_coords_rotated) ** 2, axis=1)))
    return rmsd

template_pdb=['1n7d','1ijq','3p5c','3p5b','3m0c']
best_model_pdb={'Basic':'LDLR.B99990005.pdb', 'Advanced':'LDLR.B99990004.pdb'}

for (key,val) in best_model_pdb.items():
    print(key + ' modelling -')
    print()
    for temp in template_pdb:
        print('RMSD Value for '+temp+": ",calculate_rmsd(temp+'.pdb', val))
    print()