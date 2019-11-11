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
matplotlib.rcParams['figure.max_open_warning'] = 50;

# def plot_bands_with_transport(esp, eF, energia, centralShape, V_shift = 100, porcent = 0.25, Nkx = 501, lead = 0, gammaLead =  36.917):
#
#     """
#     esp = "97", "103" ou "110"
#     eF = campo elétrico aplicado na direção-z
#     energia = nível de Fermi (para calcular correntes)
#     """
#
#     path_fig = "/home/marcos/Desktop/projetos_trabalho/images/";
#     folder = esp + "_Angstroms_GaSb_InAs/metallic_leads/";
#     name_fig = esp + "_bands_transport_eF_" \
#             + str(eF) + "meV_mu_" \
#             + str(energia) + "meV"\
#             + "_VShift_" + str(V_shift)\
#             + "_lead_" + str(lead) + ".png";
#
#     params_raw = eval("gasb.params_" + esp)
#     params_dict = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
#     hamiltonian_syst = eval("gasb.hamiltonian_" + esp + "_k_plus()")
#     # hamiltonian_lead = eval("gasb.hamiltonian_" + esp + "_k_plus(" + str(V_shift) + ")")
#     hamiltonian_lead = gasb.free_ham(norbs = 6)
#     sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, centralShape)
#
#     vec_k_limited_disc = np.linspace(-1, 1, Nkx) * np.pi * porcent;
#     vec_k_limited_cont = np.linspace(-1, 1, Nkx) * np.pi/shapes.A_STD * porcent;
#
#     fig1 = plt.figure(figsize=(10,10))
#     ax1 = fig1.add_subplot(121)
#     ax2 = fig1.add_subplot(322)
#     ax3 = fig1.add_subplot(324)
#     ax4 = fig1.add_subplot(326)
#
#     "Bands: "
#     cont_energies = trans.continuous_bands_2D(kx_array = vec_k_limited_cont,
#                                     ky_value = 0,
#                                     hamiltonian = hamiltonian_syst,
#                                     params = params_dict,
#                                     eF_value = eF)
#     disc_energies = trans.band_values(ham_syst = hamiltonian_syst,
#                                     momenta = vec_k_limited_disc,
#                                     params = params_dict,
#                                     eF_value = eF)
#     trans.bands_cont2D_and_discr(axis = ax1,
#                                 cont_energies = cont_energies,
#                                 disc_energies = disc_energies,
#                                 momenta = vec_k_limited_disc,
#                                 kx_max = 0.20,
#                                 E_min = 405,
#                                 E_max = 460,
#                                 E_line = energia)
#
#     "Transport: "
#     trans.current_spin(ax2, sistema, params_dict,
#                     eF_value = eF,
#                     energy = energia,
#                     lead_index=lead,
#                     spin='total',
#                     colormap="Oranges")
#     trans.current_spin(ax3, sistema, params_dict,
#                     eF_value = eF,
#                     energy = energia,
#                     lead_index=lead,
#                     spin='up',
#                     colormap="Reds")
#     trans.current_spin(ax4, sistema, params_dict,
#                     eF_value = eF,
#                     energy = energia,
#                     lead_index=lead,
#                     spin='down',
#                     colormap="Blues")
#     plt.tight_layout()
#     plt.savefig(path_fig + folder + name_fig)
#     # plt.show()
#     return 0
#
#
# def read_file_and_prepare_many_eFs():
#     """
#     This function reads the file that carries the input information
#     and return lists that may be iterated to automatically perform
#     calculations in a unch of systems at once.
#
#     If the file was not passed the "try-except" block will handle of
#     this exception.
#
#     The file has to have the following format (WITHOUT THE COMENTS):
#
#     3               # number of system's configurations (: different widths or different Fermi levels): positive integer
#     103             # width of the system: 97, 103 or 110 are accepted
#     0 10 20 30      # values of "eF" for the first system
#     430             # Fermi level for the first system
#     103             # Start of the second system: width
#     0 2.5 5 7.5     # values of "eF" for the second system
#     415             # Fermi level for the second system
#     97              # Start of the third system: width
#     0 10 20 30 40   # values of "eF" for the third system
#     435             # Fermi level for the third system
#     """
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
#
#     sistemas = []
#     eFs = []
#     mus = []
#
#     for i in range(numero_de_entradas):
#         espessura = file.readline()
#         eF_values = [float(x) for x in list(file.readline().split())]
#         mu_values = float(file.readline())
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
# def read_file_and_prepare_many_Fermi():
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
#         eF_values = float(file.readline())
#         N_mus = int(file.readline())
#         mu_ends = [float(x) for x in list(file.readline().split())]
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
# def calculateForManyeFs(shape):
#     """
#
#     Use this function when the INPUT_FILE.txt follows the format discribed
#     by (WITHOUT THE COMENTS):
#
#     3               # number of system's configurations : positive integer
#     103             # width of the system: 97, 103 or 110 are accepted
#     0 10 20 30      # values of "eF" for the first system
#     430             # Fermi level for the first system
#     103             # Start of the second system: width
#     0 2.5 5 7.5     # values of "eF" for the second system
#     415             # Fermi level for the second system
#     97              # Start of the third system: width
#     0 10 20 30 40   # values of "eF" for the third system
#     435             # Fermi level for the third system
#
#     """
#
#     n_sist, sists, eFields, eFermis = read_file_and_prepare_many_eFs()
#
#     for i in range(n_sist):
#         espessura = sists[i]
#         energia = eFermis[i]
#
#         for eF in eFields[i]:
#             plot_bands_with_transport(
#                 esp          = espessura,  # 97, 103 or 110
#                 eF           = eF,         # electric field normal to the system
#                 energia      = energia,    # Fermi level
#                 centralShape = shape,      # geometry of scattering reagion
#                 V_shift      = 100,       # Potential on-site on the leads
#                 porcent      = 0.25,       # Brillouin zone fraction shown
#                 Nkx          = 501,        # number of points to form the band
#                 lead         = 0           # lead where the incident wave come from
#             )
#
#
# def calculateForManyFermiEnergies(shape):
#     """
#
#     Use this function when the INPUT_FILE.txt follows the format discribed
#     by (WITHOUT THE COMENTS):
#
#     3               # number of system's configurations: positive integer
#     103             # width of the system: 97, 103 or 110 are accepted
#     30              # value of "eF" for the first system
#     11              # number of Fermi energies to use in first system
#     436 438         # min. and max. of Fermi energies for 1st system
#     97              # width of the 2nd system
#     60              # value of "eF" for the 2nd system
#     11              # number of Fermi energies to use in 2nd system
#     437 440         # min. and max. of Fermi energies for 2nd system
#     110             # width of the 3rd system
#     30              # value of "eF" for the 3rd system
#     11              # number of Fermi energies to use in 3rd system
#     428 430         # min. and max. of Fermi energies for 3rd system
#
#     """
#
#     n_sist, sists, eFields, eFermis = read_file_and_prepare_many_Fermi()
#
#     for i in range(n_sist):
#         espessura = sists[i]
#         eF = eFields[i]
#
#         for e_Fermi in eFermis[i]:
#             plot_bands_with_transport(
#                 esp          = espessura, # 97, 103 or 110
#                 eF           = eF,        # electric field normal to the system
#                 energia      = e_Fermi,   # Fermi level
#                 centralShape = shape,     # geometry of scattering reagion
#                 V_shift      = 100,      # Potential on-site on the leads
#                 porcent      = 0.25,      # Brillouin zone fraction shown
#                 Nkx          = 501,       # number of points to form the band
#                 lead         = 0          # lead where the incident wave come from
#             )

# Formatação para os gráficos:
FONT_LABELS = 18
FONT_TITLES = 20
font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
matplotlib.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



def formatter_current_axis(value, tick_number):
    #Copiado para phd
    "Function for changing units of axis"
    return int(round(value * shapes.A0/10))



def edit_axis(axis, spin):
    # Copiado para phd
    "Editing axis"
    axis.set_ylim(-shapes.W_STD/2, shapes.W_STD/2)
    axis.xaxis.set_major_locator(plt.MultipleLocator(shapes.L_STD/2))
    axis.yaxis.set_major_locator(plt.MultipleLocator(shapes.W_STD/2))

    if spin.lower() == 'total':
        axis.set_title("(b)", fontsize=FONT_TITLES)
    elif spin.lower() == 'up':
        axis.set_title("(c)", fontsize=FONT_TITLES)
    elif spin.lower() == 'down':
        axis.set_title("(d)", fontsize=FONT_TITLES)

    if (spin.lower() == 'total') or (spin.lower() == 'up'):
        axis.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
        axis.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
    else:
        axis.xaxis.set_major_formatter(plt.FuncFormatter(formatter_current_axis))
        axis.yaxis.set_major_formatter(plt.FuncFormatter(formatter_current_axis))
        axis.set_xlabel(r'$x$ [nm]',fontsize=FONT_LABELS)
        axis.set_ylabel(r'$y$ [nm]',fontsize=FONT_LABELS)



def main():
    """
    This code is for generate maps of current, without bands structures and
    with custom operators. This will allow a better exploratory process.
    """
    # Define the system
    hamiltonian = gasb.hamiltonian_97_k_plus()
    lead_ham = gasb.free_ham(norbs = 6)
    centralShape = shapes.Rect()
    syst = gasb.system_builder(hamiltonian, lead_ham, centralShape)

    # Calculate the wave function:
    energia = 442
    parametros = gasb.params_97
    parametros['eF'] = 60
    parametros['Eta2'] = 0
    parametros['Eta3'] = 0
    parametros = dict(GammaLead = parametros["GammaC"], V = 100, **parametros)
    wf = kwant.wave_function(syst, energy=energia, params=parametros)
    # modes = wf_dn(0) # from left lead

    # Define the operator
    σz = tinyarray.array([[1,0],[0,-1]])
    Mz = np.kron(σz, np.eye(3))

    # Current
    colormap = "Reds"
    J_spin = kwant.operator.Current(syst, Mz)
    current_spin = sum(J_spin(psi, params = parametros) for psi in wf(0))


    # Plot and/or save the data
    fig, axis =  plt.subplots(1,1, figsize=(8,8))
    kwant.plotter.current(syst, current_spin,
                            cmap = colormap,
                            colorbar = False,
                            show = False,
                            ax = axis)
    edit_axis(axis)
    plt.show()






if __name__ == '__main__':
    main()