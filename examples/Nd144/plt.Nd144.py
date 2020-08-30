import numpy as np
import matplotlib.pyplot as plt
import matplotlib

font = { 'size'   : 18, 'family' : 'serif'}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

dat = np.loadtxt('Nd144.out')

plt.plot(dat[1:,3], dat[1:,19], 'bo-', label='neutron')
plt.plot(dat[1:,3], dat[1:,18], 'rs--', label='proton')

plt.tick_params(axis='both', direction='in')

plt.xlim(1, 6.0)
plt.ylabel(r'$\langle uv \Delta \rangle_q$ (MeV)')
plt.xlabel(r'$\beta$ (MeV$^{-1}$)')
plt.legend()


plt.savefig('Nd144.transition.eps', bbox_inches='tight')

