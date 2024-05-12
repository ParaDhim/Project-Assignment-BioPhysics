from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

#List of protein structure files (PDB format)
protein_files = ['1n7d.pdb']  # Add your files

#Loop through the protein structure files
for pdb_file in protein_files:
    # read model file
    mdl = complete_pdb(env, pdb_file, model_segment=('FIRST:A', 'LAST:A'))

# Assess all atoms with DOPE:
s = Selection(mdl)
output_file = pdb_file.replace('.pdb', '_profile.profile')
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='LDLR.profile',
              normalize_profile=True, smoothing_window=15)


#model 5 is best as it has minimum dope score and maximum GA341 score