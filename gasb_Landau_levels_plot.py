#!/usr/bin/env python3

# Python
import sys
from itertools import product

# Anaconda
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Kwant
import kwant
import tinyarray
import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes

def main():

    # Define the system
    parameters  = gasb.params_97
    hamiltonian = gasb.hamiltonian_97_k_plus()
    Width       = shapes.W_STD
    a_lattice   = shapes.A_STD
    syst_inf  = gasb.just_lead_builder(hamiltonian, Width, a_lattice, symmetry=1)

    # Calculate the energies at kx == 0
    min_field     = 0
    max_field     = 5000
    N_pts         = 5001
    efield_values = np.linspace(min_field, max_field, N_pts)
    energies_conf = trans.collected_energies(syst_inf, parameters, efield_values)
    energies_free = trans.continuous_levels_eF(0, 0, hamiltonian, parameters, efield_values)

    # Save the results
    path_data  = './data/band_structure/'
    name       = '97_landau_levels_' + 'eF_min_' + str(min_field) + \
                '_max_' + str(max_field) + '_N_' + str(N_pts)
    np.savez(path_data + name, efields = efield_values,
             Energies_conf = energies_conf,
             Energies_free = energies_free)

    # Plot the results
    "Formatação para os gráficos:"
    FONT_LABELS = 18
    FONT_TITLES = 20
    font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
    matplotlib.rc('font', **font)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    "Actual plot:"
    fig, ax = plt.subplots(figsize=(7,8))
    ax.plot(efield_values, energies_conf, color='k', linestyle='', marker=',', linewidth=0.7)
    ax.plot(efield_values, energies_free, color='b', linewidth=2)
    ax.set_xlabel(r'eF [meV]',fontsize=FONT_TITLES)
    ax.set_ylabel(r'$\varepsilon$ [meV]',fontsize=FONT_TITLES)
    # ax.set_ylim(410, 550)
    ax.set_xlim(min_field, max_field)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
