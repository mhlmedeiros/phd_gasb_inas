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

    # Set params
    parametros = gasb.params_97
    parametros["Eta3"] = 0
    parametros["Eta2"] = 0
    parametros["eF"] = 50

    # Calculate bands
    Nkx = 200
    porcent = 0.15
    momenta_up = np.linspace(-1, 1, Nkx) * np.pi * porcent
    momenta_dn = np.linspace(-1, 1, Nkx) * np.pi * porcent


    bands_up = kwant.physics.Bands(syst_lead_up, params = parametros)
    bands_dn = kwant.physics.Bands(syst_lead_dn, params = parametros)

    energies_up = [bands_up(k) for k in momenta_up]
    energies_dn = [bands_dn(k) for k in momenta_dn]

    # Plot bands
    fig, ax = plt.subplots(1, 2, figsize=(16,8))
    tools.band_with_line_gasb(ax[0], momenta_up, energies_up,
                            kx_max = 0.15, E_max = 470, E_min = 400)
    tools.band_with_line_gasb(ax[1], momenta_dn, energies_dn,
                            kx_max = 0.15, E_max = 470, E_min = 400)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
