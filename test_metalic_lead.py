import numpy as np
import matplotlib.pyplot as plt

# # Kwant
# import kwant
# import tinyarray
# import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes

def main():
    syst_format = shapes.Rect()
    syst_hamiltonian = gasb.hamiltonian_97_k_plus()
    syst = gasb.system_builder(syst_hamiltonian, syst_format, norbs = 6)

    Nkx = 200
    porcent = 0.25
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent
    params = gasb.params_97
    params["Eta3"] = 0
    params["Eta2"] = 0

    new_params = dict(GammaLead = params["GammaC"], ShiftLead = -0.00, **params)


    # bands = trans.band_values(syst, momenta, new_params, eF_value = 0)


    fig, axis = plt.subplots(1, 1, figsize = (10,8))

    trans.current_density(axis, syst, new_params, eF_value = 60, energy = 442, lead_index = 0, colormap="Oranges")
    plt.show()


    # plt.plot(momenta, bands)
    # plt.grid()
    # plt.xlim(-.2,.2)
    # plt.ylim(400,500)
    # plt.show()

if __name__ == "__main__":
    main()
