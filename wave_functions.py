'''
Script for analysis of wavefunctions on GaSb/InAs/GaSb simmetric quantum wells.

This piece code is part of the project "phd_gasb_inas", which comprises the work
related to the Phd. Dissertation named: "Quantum transport of charge and spin in
topological insulators 2D".

Author: Marcos Medeiros
email: mhlmedeiros@gmail.com
'''
#!/usr/bin/env python

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

    # Calculate the density:
    rho = kwant.operator.Density(syst)
    psi_sqrd = sum(rho(p) for p in psi)

    # Plot the results:
    plot_density(syst, psi_sqrd)


if __name__ == '__main__':
    main()
