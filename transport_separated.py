#!/usr/bin/env python3

# Python
import sys
from itertools import product

# Anaconda
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Kwant
import kwant
import tinyarray
import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes


def calc_transp_Fermi(system, Fermi_energies, parameters):

    transport_up = np.empty(len(Fermi_energies))
    transport_dn = np.empty(len(Fermi_energies))
    transport_total = np.empty(len(Fermi_energies))

    for i in range(len(Fermi_energies)):

        # compute the scattering matrix at a given energy
        smatrix = kwant.smatrix(system, Fermi_energies[i], params = parameters)

        # compute the transmission probability from lead 0 to lead 1
        # transport_up[i] = smatrix.transmission((1, 0), (0, 0))
        # transport_dn[i] = smatrix.transmission((1, 1), (0, 1))
        transport_total[i] = smatrix.transmission(1, 0)

    return transport_up, transport_dn, transport_total

def calc_transp_Field(system, e_field_values, parameters, Fermi_energy=440):

    transport_up = np.empty(len(e_field_values))
    transport_dn = np.empty(len(e_field_values))
    transport_total = np.empty(len(e_field_values))

    for i in range(len(e_field_values)):
        # update electric field value
        params_dict['eF'] = e_field_values[i]

        # compute the scattering matrix at a given Fermi energy
        smatrix = kwant.smatrix(system, Fermi_energy, params = parameters)

        # compute the transmission probability from lead 0 to lead 1
        transport_up[i] = smatrix.transmission((1, 0), (0, 0))
        transport_dn[i] = smatrix.transmission((1, 1), (0, 1))
        transport_total[i] = smatrix.transmission(1, 0)

    return transport_up, transport_dn, transport_total

def main():

    shapeScattering = shapes.Rect(shapes.W_STD, shapes.L_STD)

    params_raw       = gasb.params_97
    gammaLead =  36.917
    V_shift = 100
    params_dict      = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
    hamiltonian_syst = gasb.hamiltonian_97_k_plus()

    hamiltonian_lead = gasb.free_ham(norbs = 6)
    sistema       = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, shapeScattering)
    # sistema_dn       = gasb.system_builder(hamiltonian_syst_dn, hamiltonian_lead, shapeScattering)

    Fermi_initial, Fermi_final, N_values = 435, 440, 2001
    Fermi_values = np.linspace(Fermi_initial, Fermi_final, N_values)

    params_dict['eF'] = 62.

    # " Names for the files: "
    path_data = "./data/transport/"
    name_preffixe = "data_" + str(Fermi_initial) + "_" + str(Fermi_final) + "_meV_Fermi_"
    name_Fermi = name_preffixe + str(N_values)
    name_transport_up = name_preffixe + "Transport_UP_" + str(N_values)
    name_transport_dn = name_preffixe + "Transport_DN_" + str(N_values)
    name_transport_total = name_preffixe + "Transport_Total_" + str(N_values)
    np.save(path_data + name_Fermi + ".npy", Fermi_values)
    np.savetxt(path_data + name_Fermi + ".txt", Fermi_values)

    # " Transport: "
    transport_up, transport_dn, transport_total = calc_transp_Fermi(sistema, Fermi_values, params_dict)

    # Save the results in "numpy" format and "txt"
    # np.save(path_data + name_transport_up + ".npy", transport_up)
    # np.savetxt(path_data + name_transport_up + ".txt", transport_up)
    #
    # np.save(path_data + name_transport_dn + ".npy", transport_dn)
    # np.savetxt(path_data + name_transport_dn + ".txt", transport_dn)

    np.save(path_data + name_transport_total + ".npy", transport_total)
    np.savetxt(path_data + name_transport_total + ".txt", transport_total)

if __name__ == '__main__':
    main()
