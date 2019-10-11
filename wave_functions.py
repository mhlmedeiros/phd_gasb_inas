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

# local application imports
from hamiltonians import gasb_hamiltonian as gasb
from system_geometry import shapes
from transport_tools import bands_and_currents as tools


def plot_density(syst, psi_sqrd):
    # Plot the results:
    fig, ax1 = plt.subplots(1,1)
    kwant.plotter.map(syst, psi_sqrd, ax=ax1)
    tools.edit_axis(ax1,'dens')
    plt.tight_layout()
    plt.show()
    return 0

def density_in_line(syst, states):

    def line(site):
        (x, y) = site.pos
        half   = shapes.L_STD/2
        delta  = shapes.A_STD
        return abs(x - half) <= delta

    rho_line = kwant.operator.Density(syst, where = line, sum = False)
    dos_line = np.array([rho_line(p) for p in states])
    return dos_line

def plot_dos_in_line(dos_line):
    fig, ax = plt.subplots(1, 2, figsize = (10,5))
    ax[0].plot(dos_line[0], color = 'red')
    ax[1].plot(dos_line[1], color = 'blue')
    plt.tight_layout()
    plt.show()


def main():
    # Define the system:
    hamiltonian = gasb.hamiltonian_103_k_plus()
    centralShape = shapes.Rect()
    syst = gasb.system_builder(hamiltonian, centralShape)

    # Calculate the wave_function:
    energia = 435
    parametros = gasb.params_103
    parametros['eF'] = 25
    wf = kwant.wave_function(syst, energy=energia, params=parametros)
    psi = wf(0)

    # # Calculate the density:
    # rho = kwant.operator.Density(syst)
    # psi_sqrd = sum(rho(p) for p in psi)
    #
    # # Plot the results:
    # plot_density(syst, psi_sqrd)

    # Calculate dos in a line
    dos_in_line_from_left  = density_in_line(syst, psi)
    dos_in_line_from_both  = density_in_line(syst, np.vstack((wf(0),wf(1))))
    plt.plot(sum(dos_in_line_from_both))
    plt.show()
    # print(dos_in_line.shape)
    # print(dos_in_line)
    # plot_dos_in_line(dos_in_line_from_left)
    # plot_dos_in_line(dos_in_line_from_both)
    print(sum(dos_in_line_from_both).shape)



if __name__ == '__main__':
    main()
