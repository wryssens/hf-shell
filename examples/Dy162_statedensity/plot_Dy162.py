import numpy as np
import matplotlib.pyplot as plt
import matplotlib
################################################################################
# Some libraries that are necessary for the inset axis. Feel free to disable 
# if you do not have the right libraries installed.
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.legend_handler import HandlerLineCollection
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
################################################################################
# Since we have no knowledge of the configuration of your system, we manually
# bring the level_densities.py auxiliary script here.
import os
os.system('cp ../../scripts/level_densities.py .')
from level_densities import *
################################################################################

font = {'size'   : 18, 'family' : 'serif'}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

A = 162
fname = 'Dy162.out'

# We generate the relevant quantities from the HF-SHELL output
(E_c, S_c, rho_s, rho_l) = generate_level_densities(fname)
E_c = E_c - E_c[0]

fig, ax = plt.subplots(1)
axins = zoomed_inset_axes(ax, 4, loc=4) 

dat = np.genfromtxt(fname)
ind = np.where(dat[:,3] > 0.815)[0]
ind2= np.where(dat[:,3] < 0.852)[0] 
# Plot the dashed line for the shape transition
axins.plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], 'k--', lw=1.0)
ax.plot([E_c[ind2[0]], E_c[ind[-1]]], [rho_s[ind2[0]], rho_s[ind[-1]]], 'k--', lw=1.0)
# Remove numerical artefacts from the shape transition, where the canonical 
# quantities are discontinous.
for k in ind:
  if k in ind2:
    E_c[k] = np.nan
# Remove numerical artefacts from the pairing transition
ind = np.where(dat[:,3] > 2.9)[0]
ind2= np.where(dat[:,3] < 4.0)[0] 
# Plot the dashed line for the pairing transition
axins.plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], 'k--', lw=1.0)
ax.plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], 'k--', lw=1.0)


# Remove numerical artefacts from the shape transition, where the canonical 
# quantities are discontinous.
for k in ind:
  if k in ind2:
    E_c[k] = np.nan

# Plot the remaining data
ax.semilogy(E_c - E_c[0], rho_s, 'k', ms=0.2)
axins.semilogy(E_c - E_c[0], rho_s, 'k', ms=0.2)

# All the data is plotted, now we just add all relevant plotting.
axins.set_xlim(0.0,2.3)
axins.set_ylim(0.5,300)

ax.tick_params(axis='both', direction='in', which='both')
axins.tick_params(axis='both', direction='in', which='both')


axins.set_xticks(np.arange(0,2.5,1.0))

axins.xaxis.tick_top()
axins.set_xlabel(r'$E_x$ (MeV)')
axins.set_ylabel(r'$\rho$ (MeV$^{-1}$)')

axins.xaxis.set_label_position('top')

ax.set_xlim(0.0,35)
ax.set_ylim(0.1,10e18)

alpha = "0.7"
mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec=alpha, ls='--')

ax.set_xlabel(r'$E_x$ (MeV)')
ax.set_ylabel(r'$\rho$ (MeV$^{-1}$)')

plt.savefig('Dy162.eps', bbox_inches='tight')
plt.close()
os.system('rm level_densities.py')
os.system('rm level_densities.pyc')

#-------------------------------------------------------------------------------
