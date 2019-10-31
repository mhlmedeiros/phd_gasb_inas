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
    ham_up       = gasb.hamiltonian_97_up()
    ham_dn       = gasb.hamiltonian_97_down()
    ham_total  = gasb.hamiltonian_97_k_plus()

    # Define the system
    symm_dir = 1 # symmetry direction
    syst_lead_up = gasb.just_lead_builder(ham_up, symmetry = symm_dir)
    syst_lead_dn = gasb.just_lead_builder(ham_dn, symmetry = symm_dir)
    syst_total = gasb.just_lead_builder(ham_total)

    # Set params
    eF = 60
    parametros = gasb.params_97
    parametros["eF"] = eF

    # Define momentum array
    Nkx = 200
    porcent = 0.15
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent

    # Define bands
    bands_coupled = kwant.physics.Bands(syst_total, params = parametros)
    bands_up = kwant.physics.Bands(syst_lead_up, params = parametros)
    bands_dn = kwant.physics.Bands(syst_lead_dn, params = parametros)

    # Calculate Bands for the stripped system
    energies_up = [bands_up(k) for k in momenta]
    energies_dn = [bands_dn(k) for k in momenta]
    energies_coupled = [bands_coupled(k) for k in momenta]

    # Calculate the Bands for the bulk system
    bulk_up = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_up, parametros, eF_value = eF)
    bulk_dn = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_dn, parametros, eF_value = eF)
    bulk_total = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_total,parametros, eF_value = eF)

    # Plot bands for strip
    E_min_plot = 420
    E_max_plot = 460

    fig, ax = plt.subplots(1, 3, figsize=(16,8))
    tools.band_with_line_gasb(ax[0], momenta, energies_up,
                            kx_max = 0.15, E_max = E_max_plot, E_min = E_min_plot)
    tools.band_with_line_gasb(ax[1], momenta, energies_dn,
                            kx_max = 0.15, E_max = E_max_plot, E_min = E_min_plot)

    tools.band_with_line_gasb(ax[2], momenta, energies_up,
                            kx_max = 0.15, E_max = E_max_plot, E_min = E_min_plot,
                            linestyle_plot = "-", color_plot = "red",
                            label_plot = "Up")
    tools.band_with_line_gasb(ax[2], momenta, energies_dn,
                            kx_max = 0.15, E_max = E_max_plot, E_min = E_min_plot,
                            linestyle_plot = "-", color_plot = "blue",
                            label_plot = "Down.")
    tools.band_with_line_gasb(ax[2], momenta, energies_coupled,
                            kx_max = 0.15, E_max = E_max_plot, E_min = E_min_plot,
                            linestyle_plot = "--", color_plot = "black",
                            marker_plot = None, label_plot = "Coupled")

    # Plot bands for bulk
    momenta_bulk = tools.trans_momenta(momenta)
    ax[0].plot(momenta_bulk, bulk_up, linewidth=2.5)
    # ax[1].plot([0.03733],[440],marker = "o", markersize = 7, color = "blue")
    ax[0].set_title("spin-up")

    ax[1].plot(momenta_bulk, bulk_dn, linewidth=2.5)
    # ax[0].plot([0.02784],[440], marker = "o", markersize = 7, color = "red")
    ax[1].set_title("spin-down")

    ax[2].set_title("both spins")
    ax[2].legend(fontsize = 12, bbox_to_anchor = (0.68, 0.64))

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
