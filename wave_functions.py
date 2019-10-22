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
    kwant.plotter.map(syst, psi_sqrd, ax = ax, fig_size = (7,3), cmap=colormap)
    tools.edit_axis(ax,'dens')
    return 0

def density_in_line(syst, states, Op = np.eye(6)):
    y_stack = []
    def line(site):
        (x, y) = site.pos
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
    # return sum(dos_in_line)
    return sum(dos_in_line)/max(sum(dos_in_line))

def print_info_dos_line(y_values, dos_in_line):
    print(80*"=")
    print("Size of dos_both: ", dos_in_line.shape)
    print("Size of y_both: ", y_values.shape)
    print("y_both:\n", y_values)



def main():
    # Define the system:
    # hamiltonian = gasb.hamiltonian_97_k_plus()
    hamiltonian = gasb.free_ham(6)
    lead_ham = gasb.free_ham(6)
    centralShape = shapes.Rect()
    syst = gasb.system_builder(hamiltonian, lead_ham, centralShape)


    # Calculate the wave_function:
    energia = 440
    parametros = gasb.params_97
    parametros['Eta3'] = 0
    parametros['Eta2'] = 0
    parametros['eF']   = 50
    parametros = dict(GammaLead = 1.5 * parametros["GammaC"],
                      ShiftLead = -0.00, **parametros )
    wf = kwant.wave_function(syst, energy=energia, params=parametros)
    modes_left  = wf(0)
    modes_right = wf(1)
    # modes_total = np.vstack((wf(0), wf(1)))


    # Calculate the density:
    sigma_z = tinyarray.array([[1,0],[0,-1]])
    spin_proj= np.kron(sigma_z, np.eye(3))
    identity = np.eye(6)
    rho = kwant.operator.Density(syst, identity)
    psi_left = sum(rho(p) for p in modes_left)
    psi_right = sum(rho(p) for p in modes_right)


    # Calculate dos in a line
    dos_in_line_from_left, y_values_left  = density_in_line(syst, modes_left, identity)
    dos_in_line_from_right, y_values_right = density_in_line(syst, modes_right, identity)


    # Plot the results:
    fig, ax = plt.subplots(2,2,figsize=(14,6))
    y_values_left = y_values_left * (shapes.A0 / 10) # conversion to nm^{-1}
    y_values_right = y_values_right * (shapes.A0 / 10) # conversion to nm^{-1}
    min_line, max_line = -0.7 * shapes.L_STD, 0.7 * shapes.L_STD

    map_density(ax[0][0], syst, psi_left, colormap = "Oranges")
    ax[0][0].vlines(0, min_line, max_line, linestyle = "--")
    ax[0][0].set_title("left lead")
    map_density(ax[1][0], syst, psi_right, colormap = "Oranges")
    ax[1][0].vlines(0, min_line, max_line, linestyle = "--")
    ax[1][0].set_title("right lead")

    ax[0][1].plot(y_values_left, normalize(dos_in_line_from_left),
                    marker = ".", markersize = 2.5, linestyle = "-" )
    ax[1][1].plot(y_values_right, normalize(dos_in_line_from_right),
                    marker = ".", markersize = 2.5, linestyle = "-" )
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
