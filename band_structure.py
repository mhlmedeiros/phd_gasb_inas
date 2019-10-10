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
    # Define the system
    ham_up       = gasb.hamiltonian_97_up()
    ham_dn       = gasb.hamiltonian_97_down()
    syst_lead_up = gasb.just_lead_builder(ham_up)
    syst_lead_dn = gasb.just_lead_builder(ham_dn)

    ham_total  = gasb.hamiltonian_97_k_plus()
    syst_total = gasb.just_lead_builder(ham_total)

    # Set params
    parametros = gasb.params_97
    parametros["Eta3"] = 0
    parametros["Eta2"] = 0
    parametros["eF"] = 50

    params_orig = gasb.params_97

    # Calculate bands
    Nkx = 200
    porcent = 0.15
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent
    # momenta_dn = np.linspace(-1, 1, Nkx) * np.pi * porcent


    bands_up = kwant.physics.Bands(syst_lead_up, params = parametros)
    bands_dn = kwant.physics.Bands(syst_lead_dn, params = parametros)

    bands_coupled = kwant.physics.Bands(syst_total, params = params_orig)

    energies_up = [bands_up(k) for k in momenta]
    energies_dn = [bands_dn(k) for k in momenta]

    energies_coupled = [bands_coupled(k) for k in momenta]

    print(np.array(energies_up).shape)
    # print(len(energies_up[0]))

    # Plot bands
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
                            linestyle_plot = "-", color_plot = "blue",
                            label_plot = "Down")
    tools.band_with_line_gasb(ax[2], momenta, energies_coupled,
                            kx_max = 0.15, E_max = 470, E_min = 400,
                            linestyle_plot = "--", color_plot = "black",
                            marker_plot = None, label_plot = "Coupled")
    ax[0].set_title("spin-up")
    ax[1].set_title("spin-down")
    ax[2].set_title("both and coupled")
    ax[2].legend(fontsize = 12, bbox_to_anchor = (0.69, 0.612))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
