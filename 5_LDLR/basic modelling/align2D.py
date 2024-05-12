from modeller import *

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='1n7d', model_segment=('FIRST:A','LAST:A'))

aln.append_model(mdl, align_codes='1n7d', atom_files='1n7d.pdb')
aln.append(file='seq.pir', align_codes='LDLR')
aln.align2d()
aln.write(file='Target-template.ali', alignment_format='PIR')
aln.write(file='Target-template.pap', alignment_format='PAP')