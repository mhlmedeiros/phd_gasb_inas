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
    ham = gasb.hamiltonian_97_k_plus()
    syst_lead = gasb.just_lead_builder(ham)

    # Set params
    parametros = gasb.params_97
    parametros["Eta3"] = 0
    parametros["Eta2"] = 0
    parametros["eF"] = 50

    # Calculate bands
    Nkx = 200
    porcent = 0.15
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent
    bands = kwant.physics.Bands(syst_lead, params = parametros)
    energies = [bands(k) for k in momenta]


    # Plot bands
    fig, ax = plt.subplots(1, 1, figsize=(8,8))
    tools.band_with_line_gasb(ax, momenta, energies,
                            kx_max = 0.15, E_max = 450, E_min = 400)


if __name__ == '__main__':
    main()
