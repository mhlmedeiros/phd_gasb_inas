#!/usr/bin/env python3

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
    ham_dn       = gasb.hamiltonian_97_up_y_inv()
    # ham_total  = gasb.hamiltonian_97_k_plus()

    # Define the system
    syst_lead_up = gasb.just_lead_builder(ham_up, symmetry = -1)
    syst_lead_dn = gasb.just_lead_builder(ham_dn, symmetry = -1)
    # syst_total = gasb.just_lead_builder(ham_total)

    # Set params
    eF = 50
    parametros = gasb.params_97
    parametros["Eta3"] = 0
    parametros["Eta2"] = 0
    parametros["eF"] = eF
    params_orig = gasb.params_97


    # Define momentum array
    Nkx = 200
    porcent = 0.15
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent

    # Define bands
    bands_up = kwant.physics.Bands(syst_lead_up, params = parametros)
    bands_dn = kwant.physics.Bands(syst_lead_dn, params = parametros)
    # bands_coupled = kwant.physics.Bands(syst_total, params = params_orig)

    # Calculate Bands for the stripped system
    energies_up = [bands_up(k) for k in momenta]
    energies_dn = [bands_dn(k) for k in momenta]
    # energies_coupled = [bands_coupled(k) for k in momenta]

    # Calculate the Bands for the bulk system
    # bulk_up = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_up, parametros, eF_value = eF)
    # bulk_dn = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_dn, parametros, eF_value = eF)
    # bulk_total = tools.continuous_bands_2D(momenta/shapes.A_STD, 0, ham_total, params_orig, eF_value = eF)

    # Plot bands for strip
    fig, ax = plt.subplots(1, 3, figsize=(16,8))
    tools.band_with_line_gasb(ax[0], momenta, energies_up,
                            kx_max = 0.15, E_max = 470, E_min = 400)
    tools.band_with_line_gasb(ax[1], momenta, energies_dn,
                            kx_max = 0.15, E_max = 470, E_min = 400)

    tools.band_with_line_gasb(ax[2], momenta, energies_up,
                            kx_max = 0.15, E_max = 470, E_min = 400,
                            linestyle_plot = "-", color_plot = "red",
                            label_plot = "Up")
    tools.band_with_line_gasb(ax[2], momenta, energies_dn,
                            kx_max = 0.15, E_max = 470, E_min = 400,
                            linestyle_plot = "--", color_plot = "blue",
                            label_plot = "Up-inv.")
    # tools.band_with_line_gasb(ax[2], momenta, energies_coupled,
    #                         kx_max = 0.15, E_max = 470, E_min = 400,
    #                         linestyle_plot = "--", color_plot = "black",
    #                         marker_plot = None, label_plot = "Coupled")

    # Plot bands for bulk
    momenta_bulk = tools.trans_momenta(momenta)
    # ax[0].plot(momenta_bulk, bulk_up, linewidth=2.5)
    # ax[1].plot([0.03733],[440],marker = "o", markersize = 7, color = "blue")
    ax[0].set_title("spin-up orig.")

    # ax[1].plot(momenta_bulk, bulk_dn, linewidth=2.5)
    # ax[0].plot([0.02784],[440], marker = "o", markersize = 7, color = "red")
    ax[1].set_title("spin-up inv.")

    ax[2].set_title("orig. and inv.")
    ax[2].legend(fontsize = 12, bbox_to_anchor = (0.69, 0.612))

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
