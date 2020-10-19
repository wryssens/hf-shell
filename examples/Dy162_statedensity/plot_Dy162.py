#-------------------------------------------------------------------------------
# This script plots the state density of 162Dy, as obtained with HF-SHELL and
# as obtained through SMMC calculations. The SMMC results were taken from the 
# following reference: 
#
#  Y. Alhassid, G. F. Bertsch, C. N. Gilbreth and H. Nakada
#  Benchmarking mean-field approximations to level densities,
#  Phys. Rev. C 93, 044320 (2016).
#
#-------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

################################################################################
# Some libraries that are necessary to draw the inset axis. Feel free to disable 
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
font = {'size'   : 23, 'family' : 'serif'}
matplotlib.rc('font', **font)
matplotlib.rc('text', usetex=True)
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

A = 162
fname = 'Dy162.out'

# We generate the relevant quantities from the HF-SHELL output
(E_c, S_c, rho_s, rho_l) = generate_level_densities(fname)
E_c = E_c - E_c[0]

fig, ax = plt.subplots(2, figsize=(7,10.5))
axins      = zoomed_inset_axes(ax[0], 3.5, loc=4)

plt.subplots_adjust(hspace=0.01)

dat = np.genfromtxt(fname)
ind = np.where(dat[:,3] > 0.815)[0]
ind2= np.where(dat[:,3] < 0.852)[0] 
# Plot the dashed line for the shape transition
axins.plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], \
                                                                  'k--', lw=1.0)
ax[0].plot(   [E_c[ind2[0]], E_c[ind[-1]]], [rho_s[ind2[0]], rho_s[ind[-1]]], \
                                                                  'k--', lw=1.0)

# Remove numerical artefacts from the shape transition, where the canonical 
# quantities are discontinous.
for k in ind:
  if k in ind2:
    E_c[k] = np.nan

# Remove numerical artefacts from the pairing transition
ind = np.where(dat[:,3] > 2.9)[0]
ind2= np.where(dat[:,3] < 4.0)[0] 

# Plot the dashed line for the pairing transition
axins.plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], \
                                                                  'k--', lw=1.0)
ax[0].plot([E_c[ind[-1]], E_c[ind2[0]]], [rho_s[ind[-1]], rho_s[ind2[0]]], \
                                                                  'k--', lw=1.0)
# Remove numerical artefacts from the shape transition, where the canonical 
# quantities are discontinous.
for k in ind:
  if k in ind2:
    E_c[k] = np.nan

# Plot the remaining data
ax[0].semilogy(E_c - E_c[0], rho_s, 'k', ms=0.2, label='FTHFB')
axins.semilogy(E_c - E_c[0], rho_s, 'k', ms=0.2, label='FTHFB')


# Plot the SMMC data
dat = np.loadtxt('SMMC.dat')
ax[0].errorbar(dat[:,0], np.exp(dat[:,1]), yerr = np.exp(dat[:,1])*dat[:,2], 
           marker='o', color='b', linestyle='', ms=4, mfc='white', label='SMMC')
axins.errorbar(dat[:,0], np.exp(dat[:,1]), yerr = np.exp(dat[:,1])*dat[:,2], 
               marker='o', color='b', linestyle='', ms=4, mfc='white', zorder=0)

# COnstruct an interpolator for the mean-field results
f = np.interp(dat[::-1,0],E_c - E_c[0], rho_s)
# Plot the ratio of state densities
ax[1].errorbar(dat[:,0], np.exp(dat[:,1])/f[::-1],  \
               yerr = np.exp(dat[:,1])/f[::-1] * dat[:,2], color='b', \
               marker='o', linestyle=' ', mfc='white')
# Plot the reference line
ax[1].plot([0,55], [1.0, 1.0], color='0.7')
#-------------------------------------------------------------------------------
# All the data is plotted, now we just add all relevant plotting details.
axins.set_xlim(-0.05,3.5)
axins.set_ylim(  0.5,8000)

ax[0].tick_params(axis='both', direction='in', which='both')
ax[1].tick_params(axis='both', direction='in', which='both')
axins.tick_params(axis='both', direction='in', which='both')

axins.set_xticks(np.arange(0,3.5,1.0))
ax[0].set_xticklabels([])
axins.xaxis.tick_top()
axins.set_xlabel(r'$E_x$ (MeV)')
axins.set_ylabel(r'$\rho$ (MeV$^{-1}$)')

axins.xaxis.set_label_position('top')

ax[0].set_xlim(0.0,55)
ax[0].set_ylim(0.21,10e20)

ax[1].set_xlim(0, 55)
ax[1].set_ylim(0.5, 120)
ax[1].set_yscale('log')

alpha = "0.7"
mark_inset(ax[0], axins, loc1=2, loc2=3, fc="none", ec=alpha, ls='--')

ax[0].set_ylabel(r'$\rho$ (MeV$^{-1}$)')

ax[1].set_xlabel(r'$E_x$ (MeV)')
ax[1].set_ylabel(r'$\rho_{\rm SMMC} / \rho_{\rm FTHFB}$')

ax[0].legend(loc='upper left', fontsize=22)
plt.savefig('Dy162.eps', bbox_inches='tight')
plt.close()
os.system('rm level_densities.py')
os.system('rm level_densities.pyc')
#-------------------------------------------------------------------------------
