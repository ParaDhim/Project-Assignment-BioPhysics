from modeller import *
from modeller.automodel import *
# from modeller import soap_protein_od

env = Environ()
a = AutoModel(env, alnfile='Target-template.ali',
              knowns='1n7d', sequence='LDLR',
              assess_methods=(assess.DOPE,
                            #   soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 10
a.make()