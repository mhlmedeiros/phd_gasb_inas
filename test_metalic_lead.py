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

    ## Define system:
    syst_format = shapes.Rect()
    syst_hamiltonian = gasb.hamiltonian_97_up()
    lead_ham = gasb.hamiltonian_97_up(100)
    syst = gasb.system_builder(syst_hamiltonian, lead_ham, syst_format)
    # kwant.plot(syst)

    ## Define the parameters of the Hamiltonians (for system and leads)
    params = gasb.params_97
    params["Eta3"] = 0
    params["Eta2"] = 0
    new_params = dict(GammaLead = 1.5 * params["GammaC"], V = 0., **params)

    # ## Calculate density
    # wf = kwant.wave_function(syst, energy = energia, params = new_params)
    # modes_left = wf(0)
    # # sigma_z = tinyarray.array([[1,0],[0,-1]])
    # # spin_proj= np.kron(sigma_z, np.eye(3))
    # identity = np.eye(3)
    # rho = kwant.operator.Density(syst, identity)
    # psi_left = sum(rho(p) for p in modes_left)
    #
    # fig_dens, axis_dens = plt.subplots(1, 1, figsize = (10,8))
    # kwant.plotter.map(syst, psi_left, ax = axis_dens, cmap = "Oranges")
    energia = 442
    eF = 60
    fig_current, axis = plt.subplots(1, 1, figsize = (8,5))
    trans.current_density(axis, syst, new_params, eF_value = eF, energy = energia, lead_index = 0, colormap="Reds")
    plt.show()

    # Nkx = 200
    # porcent = 1.
    # momenta = np.linspace(-1, 1, Nkx) * np.pi * porcent
    # bands = trans.band_values(syst, momenta, new_params, eF_value = 0)
    # fig_bands, axis_band = plt.subplots(1, 1, figsize = (10,8))
    # axis_band.plot(momenta, bands)
    # # plt.xlim(-.2,.2)
    # # plt.ylim(400,500)
    # plt.grid()
    # plt.show()

if __name__ == "__main__":
    main()
