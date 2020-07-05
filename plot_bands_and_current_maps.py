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
matplotlib.rcParams['figure.max_open_warning'] = 50

path_geral = "/home/marcos/Desktop/projetos_trabalho/"
path_fig   = path_geral + "images/" # place to figures
path_data  = path_geral + "data_hdf5_phd_gasb_inas/" # place to figures
folder_suf = "_Angstroms_GaSb_InAs/metallic_leads/fermi_fixed_electric_field_var/"

# For remaking the system:
def make_system(esp="97", gammaLead=36.917, V_shift=100, width = shapes.W_STD, length = shapes.L_STD):
    shapeScattering = shapes.Rect(width, length)

    # folder_fig = esp + folder_suf

    params_raw = eval("gasb.params_" + esp)
    params_dict = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
    hamiltonian_syst = eval("gasb.hamiltonian_" + esp + "_k_plus()")

    hamiltonian_lead = gasb.free_ham(norbs = 6)
    sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, shapeScattering)

    return sistema


# For loading the data from the file:
def load_currents(path):
    data = np.load(path)
    total, up, down = data["Total"], data["Up"], data["Dn"]
    return total, up, down

# For the plotting the current map:
def plot_map(path,spin="Up",colormap="Reds",axis=None):
    data = np.load(path)

    syst = make_system()
    if axis == None:
        fig, axis = plt.subplots()

    kwant.plotter.current(syst, data[spin], cmap = colormap, colorbar = False, show = False, ax=axis)
    tools.edit_axis(axis,'total')
    axis.set_title(" ", fontsize=tools.FONT_TITLES)
    # plt.show()


def plot_bands_with_currents_from_files(path_data_bands, path_data_currents,fermi_energy):

    """

    """
    fig1 = plt.figure(figsize=(10,10))
    ax1 = fig1.add_subplot(121)
    ax2 = fig1.add_subplot(322)
    ax3 = fig1.add_subplot(324)
    ax4 = fig1.add_subplot(326)

    "Bands: "

    bands_total = np.load(path_data_bands)
    kx, E_free, E_conf = bands_total['kx'], bands_total['E_free'], bands_total['E_conf']

    trans.bands_cont2D_and_discr(axis = ax1,
                                free_elec_energies = E_free,
                                confined_elec_energies = E_conf,
                                momenta = kx,
                                kx_max = 0.12,
                                E_min = 435,
                                E_max = 455,
                                E_line = fermi_energy)

    "Currents: "
    # (Re)Make the system
    syst = make_system()

    # Load the data
    total_current, up_current, down_current = load_currents(path_data_currents)
    all_currents = [total_current, up_current, down_current]
    color_maps   = ["Oranges", "Reds", "Blues"]
    current_type = ["total","up","down"]
    axes         = [ax2, ax3, ax4]

    # Plot the maps
    for current, colormap, name, axis in zip(all_currents, color_maps, current_type, axes):
        kwant.plotter.current(syst, current, cmap = colormap, colorbar = False, show = False, ax=axis)
        trans.edit_axis(axis, name)
        # plt.savefig("../../"+name+"_current_"+str(Fermi_energy)+"L_600.png")

    plt.tight_layout()

    # "Saving the plot"
    # name_fig,_ = names(esp, eF, energia, V_shift, lead)
    # plt.savefig(path_fig + folder_fig + name_fig)
    plt.show()
    return


def main():
    path_bands    = "./data/band_structure/97_band_eF_62.0meV_0.25BZ_Nkx_501.npz"
    path_currents = "./data/local_currents/97_currents_eF_62.0meV_Fermi_439.961meV_VShift_100_lead_0_gammaLead_36.917.npz"
    energy_Fermi  = float(path_currents.split('_')[-7][:-3])

    plot_bands_with_currents_from_files(path_bands,
                                        path_currents,
                                        energy_Fermi)





if __name__ == '__main__':
    main()
