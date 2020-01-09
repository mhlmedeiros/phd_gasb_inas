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

# This following comand turns off the warnings for
# be opening many figures (default = 20)
# matplotlib.rcParams['figure.max_open_warning'] = 50
#
# path_geral = "/home/marcos/Desktop/projetos_trabalho/"
# path_fig   = path_geral + "images/" # place to figures
# path_data  = path_geral + "data_hdf5_phd_gasb_inas/" # place to figures
# folder_suf = "_Angstroms_GaSb_InAs/metallic_leads/fermi_fixed_electric_field_var/"

#
# def names_Fermi_fixo(esp, eF, energia,V_shift, lead):
#     name_geral = esp + "_total_transport_eF_" \
#         + str(energia) + "meV"\
#         + "_VShift_" + str(V_shift)\
#         + "_from_lead_" + str(lead) + ".png"
#
#     name_fig      = name_geral + ".png"
#     name_dataset  = name_geral
#
#     return name_fig, name_geral
#
#
# def read_file_and_prepare():
#
#     try:
#         infile = sys.argv[1]
#     except:
#         print("Usage: " + sys.argv[0] + " input_file_name"); sys.exit(1)
#
#     file = open(infile,'r')
#
#     numero_de_entradas = int(file.readline()) # first line will be an lonely integer
#
#     sistemas = []
#     eFs = []
#     mus = []
#
#     for i in range(numero_de_entradas):
#         espessura = file.readline()
#
#         N_eF_val  = int(file.readline())
#         eF_ends   = [float(x) for x in list(file.readline().split())]
#         eF_values = np.linspace(eF_ends[0], eF_ends[1], N_eF_val)
#
#         N_m
#
#         us     = int(file.readline())
#         mu_ends   = [float(x) for x in list(file.readline().split())]
#         mu_values = np.linspace(mu_ends[0], mu_ends[1], N_mus)
#
#         sistemas.append(espessura.rstrip('\n'))
#         eFs.append(eF_values)
#         mus.append(mu_values)
#
#     file.close();
#
#     return numero_de_entradas, sistemas, eFs, mus
#
#
# def plot_transport_total(esp, eF, energia, centralShape, V_shift = 100, porcent = 0.25, Nkx = 501, lead = 0, gammaLead =  36.917):
#
#     """
#     esp = "97", "103" ou "110"
#     eF = campo elétrico aplicado na direção-z
#     energia = nível de Fermi (para calcular correntes)
#     """
#
#     folder_fig = esp + folder_suf
#
#     params_raw = eval("gasb.params_" + esp)
#     params_dict = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
#     hamiltonian_syst = eval("gasb.hamiltonian_" + esp + "_k_plus()")
#
#     hamiltonian_lead = gasb.free_ham(norbs = 6)
#     sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, centralShape)
#
#     #
#     fig1 = plt.figure(figsize=(10,10))
#     ax1 = fig1.add_subplot(111)
#
#
#     "Transport: "
#     energies = []
#     transport = []
#     for ie in range(100):
#         energy = ie * 0.01
#
#         # compute the scattering matrix at a given energy
#         smatrix = kwant.smatrix(syst, energy)
#
#         # compute the transmission probability from lead 0 to
#         # lead 1
#         energies.append(energy)
#         transport.append(smatrix.transmission(1, 0))
#
#     "Saving the plot"
#     name_fig,_ = names(esp, eF, energia, V_shift, lead)
#     plt.savefig(path_fig + folder_fig + name_fig)
#     # plt.show()
#     return 0
#

def main():

    shapeScattering = shapes.Rect(shapes.W_STD, shapes.L_STD)

    params_raw       = gasb.params_97
    gammaLead =  36.917
    V_shift = 100
    params_dict      = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
    hamiltonian_syst = gasb.hamiltonian_97_k_plus()

    hamiltonian_lead = gasb.free_ham(norbs = 6)
    sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, shapeScattering)

    eF_initial, eF_final, N_values = 62, 64, 501
    eF_values = np.linspace(eF_initial, eF_final, N_values)

    "Transport: "
    energy = 440
    transport = np.empty(len(eF_values))
    for i in range(len(eF_values)):
        params_dict['eF'] = eF_values[i]

        # compute the scattering matrix at a given energy
        smatrix = kwant.smatrix(sistema, energy, params = params_dict)

        # compute the transmission probability from lead 0 to lead 1
        transport[i] = smatrix.transmission(1, 0)



    name_preffixe = "data_" + str(eF_initial) + "_" + str(eF_final) + "_"
    name_eF = name_preffixe + "eF_" + str(N_values)
    name_transport = name_preffixe + "Transport_" + str(N_values)

    np.save(name_eF + ".npy", eF_values)
    np.savetxt(name_eF + ".txt", eF_values)

    np.save(name_transport + ".npy", transport)
    np.savetxt(name_transport + ".txt", transport)

if __name__ == '__main__':
    main()
