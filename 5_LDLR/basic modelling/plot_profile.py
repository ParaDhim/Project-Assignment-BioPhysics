import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend, which does not require a graphical display

import matplotlib.pyplot as plt  # Import pyplot from matplotlib

import modeller

def get_profile(profile_file, seq):
    '''Read `profile_file` into a Python array, and add gaps corresponding to
       the Alignment sequence `seq.'''
    # Read all non-comment and non-blank lines from the file:
    with open(profile_file, 'r') as f:
        vals = []
        for line in f:
            if not line.startswith('#') and len(line) > 10:
                spl = line.split()
                vals.append(float(spl[-1]))
    # Insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    # Add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

e = modeller.Environ()
a = modeller.Alignment(e, file='LDLR-mult.ali')  # Use 'Alignment' instead of 'alignment'

model10 = get_profile('LDLR-model10.profile', a['MDH'])
template_1bdm = get_profile('LDLR.profile', a['3M0C'])

# Plot the template and model profiles in the same plot for comparison:
plt.figure(figsize=(20, 12))  # Create a figure using Matplotlib
plt.xlabel('Residue index')
plt.ylabel('DOPE per-residue score')
plt.plot(model10, color='blue', linewidth=2, label='Model10')
plt.plot(template_1bdm, color='red', linewidth=2, label='1bdm')
plt.legend()
plt.savefig('model10-1bdm.png', dpi=96)  # Save the plot as an image
