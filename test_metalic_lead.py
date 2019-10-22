import numpy as np
import matplotlib.pyplot as plt

# # Kwant
import kwant
# import tinyarray
# import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes

def main():
    syst_format = shapes.Rect()
    syst_hamiltonian = gasb.hamiltonian_97_up()
    lead_ham = gasb.free_ham(norbs=3)
    syst = gasb.system_builder(syst_hamiltonian, lead_ham, syst_format)
    # kwant.plot(syst)


    Nkx = 200
    porcent = 0.25
    momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent
    params = gasb.params_97
    params["Eta3"] = 0
    params["Eta2"] = 0

    new_params = dict(GammaLead = 1.5 * params["GammaC"],
                      ShiftLead = -0.00,
                      **params)


    # bands = trans.band_values(syst, momenta, new_params, eF_value = 0)


    fig, axis = plt.subplots(1, 1, figsize = (10,8))

    trans.current_density(axis, syst, new_params, eF_value = 60, energy = 442, lead_index = 0, colormap="Reds")
    plt.show()


    plt.plot(momenta, bands)
    plt.grid()
    plt.xlim(-.2,.2)
    plt.ylim(400,500)
    plt.show()

if __name__ == "__main__":
    main()
