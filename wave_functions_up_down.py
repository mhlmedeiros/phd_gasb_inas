#!/usr/bin/env python
'''
Script for analysis of wavefunctions on GaSb/InAs/GaSb simmetric quantum wells.

This piece code is part of the project "phd_gasb_inas", which comprises the work
related to the Phd. Dissertation named: "Quantum transport of charge and spin in
topological insulators 2D".

Author: Marcos Medeiros
email: mhlmedeiros@gmail.com
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Kwant related stuff
import kwant
import tinyarray

# local application imports
from hamiltonians import gasb_hamiltonian as gasb
from system_geometry import shapes
from transport_tools import bands_and_currents as tools


def map_density(ax, syst, psi_sqrd, colormap = "Reds"):
    # Plot the results:
    # fig, ax = plt.subplots(1,1)
    kwant.plotter.map(syst, psi_sqrd, ax = ax, fig_size = (7,3), cmap=colormap)
    tools.edit_axis(ax,'dens')
    # plt.tight_layout()
    # plt.show()
    return 0

def density_in_line(syst, states, Op = np.eye(3)):
    y_stack = []
    def line(site):
        (x, y) = site.pos
        # half   = shapes.L_STD/2
        half   = 0
        delta  = shapes.A_STD
        ans = abs(x - half) < delta
        if ans == True : y_stack.append(y)
        return ans

    rho_line = kwant.operator.Density(syst, Op, where = line, sum = False)
    dos_line = np.array([rho_line(p) for p in states])
    return dos_line, np.array(y_stack)

def plot_dos_in_line(dos_line):
    fig, ax = plt.subplots(1, 2, figsize = (10,5))
    ax[0].plot(dos_line[0], color = 'red')
    ax[1].plot(dos_line[1], color = 'blue')
    plt.tight_layout()
    plt.show()

def normalize(dos_in_line):
    # return sum(dos_in_line)/max(sum(dos_in_line))
    return sum(dos_in_line)

def print_info_dos_line(y_values, dos_in_line):
    print(80*"=")
    print("Size of dos_both: ", dos_in_line.shape)
    print("Size of y_both: ", y_values.shape)
    print("y_both:\n", y_values)

# def calcula_density():
#     hamiltonian_up = gasb.hamiltonian_103_up()
#     centralShape = shapes.Rect()
#     syst_up = gasb.system_builder(hamiltonian_up, centralShape)

def main():
    # Define the system:
    hamiltonian_up = gasb.hamiltonian_97_up()
    hamiltonian_dn = gasb.hamiltonian_97_up_y_inv()
    lead_ham = gasb.hamiltonian_97_up(-100)
    centralShape = shapes.Rect()
    syst_up = gasb.system_builder(hamiltonian_up, lead_ham, centralShape)
    syst_dn = gasb.system_builder(hamiltonian_dn, lead_ham, centralShape)

    # Calculate the wave_function:
    energia = 442
    parametros = gasb.params_97
    parametros['eF'] = 60
    # parametros = dict(GammaLead = parametros["GammaC"], V = 100, **parametros )
    wf_up = kwant.wave_function(syst_up, energy=energia, params=parametros)
    wf_dn = kwant.wave_function(syst_dn, energy=energia, params=parametros)
    modes_up = wf_up(0) # from left lead
    modes_dn = wf_dn(0) # from left lead


    # Calculate the density:
    rho_up = kwant.operator.Density(syst_up)
    rho_dn = kwant.operator.Density(syst_dn)
    psi_up = sum(rho_up(p) for p in modes_up)
    psi_dn = sum(rho_dn(p) for p in modes_dn)


    # Calculate dos in a line
    dos_in_line_up, y_values_up = density_in_line(syst_up, modes_up)
    dos_in_line_dn, y_values_dn = density_in_line(syst_dn, modes_dn)


    # Plot the results:
    fig, ax = plt.subplots(2, 2, figsize = (14,6))
    y_values_up = y_values_up * (shapes.A0 / 10) # conversion to nm^{-1}
    y_values_dn = y_values_dn * (shapes.A0 / 10) # conversion to nm^{-1}
    min_line, max_line = -0.7 * shapes.L_STD, 0.7 * shapes.L_STD

    map_density(ax[0][0], syst_up, psi_up, colormap = "Reds")
    ax[0][0].vlines(0, min_line, max_line, linestyle = "--")
    ax[0][0].set_title("spin up")
    map_density(ax[1][0], syst_dn, psi_dn, colormap = "Blues")
    ax[1][0].vlines(0, min_line, max_line, linestyle = "--")
    ax[1][0].set_title("spin up inv.")

    ax[0][1].plot(y_values_up, normalize(dos_in_line_up),
                    marker = ".", markersize = 2.5, linestyle = "-", color = "red")
    ax[1][1].plot(y_values_dn, normalize(dos_in_line_dn),
                    marker = ".", markersize = 2.5, linestyle = "-" )
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
