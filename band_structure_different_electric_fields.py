#!/usr/bin/env python3

"""
Here we have a script that print out the bands for the spin-UP and the spin-Down
obtained by adopting a uncoupled Hamiltonian (Eta3 == Eta2 == 0).

Each of those band structures is presented in separated panels. A third panel is
also presented with both band structures (Up and Down) and the structure  for a
coupled system for comparison.

"""
# Anaconda
import numpy as np
import matplotlib.pyplot as plt

# kwant
import kwant

# local application imports
from hamiltonians import gasb_hamiltonian as gasb
from system_geometry import shapes
from transport_tools import bands_and_currents as tools

def main():
    # Define Hamiltonians
    hamiltonian       = gasb.hamiltonian_97_k_plus()


    # Define the system
    symm_dir = 1 # symmetry direction
    syst_lead_up = gasb.just_lead_builder(hamiltonian, symmetry = symm_dir)


    # Set params
    eF0 = 0
    eF1 = 20
    eF2 = 60
    parametros = gasb.params_97

    # Define momentum array
    Nkx = 200
    porcent = 0.15
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent

    # Define bands
    parametros["eF"] = eF0
    bands_left   = kwant.physics.Bands(syst_lead_up, params = parametros)
    parametros["eF"] = eF1
    bands_center = kwant.physics.Bands(syst_lead_up, params = parametros)
    parametros["eF"] = eF2
    bands_right  = kwant.physics.Bands(syst_lead_up, params = parametros)

    # Calculate Bands for the stripped system
    energies_left   = [bands_left(k) for k in momenta]
    energies_center = [bands_center(k) for k in momenta]
    energies_right  = [bands_right(k) for k in momenta]

    # Calculate the Bands for the bulk system
    bulk_left   = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, hamiltonian, parametros, eF_value = eF0)
    bulk_center = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, hamiltonian, parametros, eF_value = eF1)
    bulk_right  = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, hamiltonian, parametros, eF_value = eF2)

    # Plot bands for strip
    E_min_plot = 420
    E_max_plot = 460
    max_kx = 0.15
    fig, ax = plt.subplots(1, 3, figsize=(14,7))
    tools.band_with_line_gasb(ax[0], momenta, energies_left,
                            kx_max = max_kx, E_max = E_max_plot, E_min = E_min_plot)
    tools.band_with_line_gasb(ax[1], momenta, energies_center,
                            kx_max = max_kx, E_max = E_max_plot, E_min = E_min_plot)
    tools.band_with_line_gasb(ax[2], momenta, energies_right,
                            kx_max = max_kx, E_max = E_max_plot, E_min = E_min_plot)

    # Plot bands for bulk
    momenta_bulk = tools.trans_momenta(momenta)
    ax[0].plot(momenta_bulk, bulk_left, linewidth=1.5, color = "red")
    ax[0].set_title("x = %.1f meV" % eF0)

    ax[1].plot(momenta_bulk, bulk_center, linewidth=1.5, color = "red")
    ax[1].set_title("x = %.1f meV" % eF1)

    ax[2].plot(momenta_bulk, bulk_right, linewidth=1.5, color = "red")
    ax[2].set_title("x = %.1f meV" % eF2)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
