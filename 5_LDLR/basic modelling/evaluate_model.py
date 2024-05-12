from modeller import *
from modeller.scripts import complete_pdb

log.verbose()    # request verbose output
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, 'LDLR.B99990005.pdb')

# Assess all atoms with DOPE:
s = Selection(mdl)
dope_score=s.assess_dope( file='LDLR-model5.profile',
              normalize_profile=True, smoothing_window=15)


print(dope_score)