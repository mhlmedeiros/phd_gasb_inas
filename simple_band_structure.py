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
    # Define the Hamiltonian
    ham_total  = gasb.free_ham(norbs = 3)

    # Define the system
    symm_dir = 1 # symmetry direction
    syst_total = gasb.just_lead_builder(ham_total, symmetry = symm_dir)
    print(syst_total)


    # Set params
    V_shift = 0
    parametros = gasb.params_110
    params_dict = dict(GammaLead = parametros["GammaC"], V = V_shift, **parametros)
    print()
    for k in params_dict.keys():
        print(k," = ", params_dict[k])

    # Define momentum array
    Nkx = 201
    porcent = 1.
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent

    # Define bands
    bands_lead = kwant.physics.Bands(syst_total, params = params_dict)
    print("ok 1")

    # Calculate Bands for the stripped system
    energies_lead = [bands_lead(k) for k in momenta]
    print("ok 2")

    # Plot bands for strip
    E_max_plot = 1000.
    E_min_plot = 0.
    kmax       = 1.
    eFermi     = 448.
    fig, ax = plt.subplots(1, 1, figsize=(8,8))
    tools.band_with_line_gasb(ax, momenta, energies_lead, E_line = eFermi,
                            kx_max = kmax, E_max = E_max_plot, E_min = E_min_plot)


    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
