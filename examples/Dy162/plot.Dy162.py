import numpy as np
import matplotlib.pyplot as plt
import matplotlib


font = { 'size'   : 18, 'family' : 'serif'}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)


dat = np.loadtxt('Dy162.out')
plt.plot(dat[1:,3], dat[1:,7], 'ks--', label='HF')
plt.tick_params(axis='both', direction='in')

plt.xlabel(r'$\beta$ (MeV$^{-1}$)')
plt.ylabel(r'$\langle Q_{20} \rangle_T$ (fm$^2$)')
plt.savefig('Dy162.transition.eps', bbox_inches='tight')
plt.close()
