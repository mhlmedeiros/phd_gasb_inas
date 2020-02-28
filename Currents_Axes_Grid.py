#!/usr/bin/env python3

# Anaconda
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

# Kwant
import kwant
import tinyarray
import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes


# In[44]:


path_data_currents = "./data/local_currents/97_currents_eF_62.0meV_Fermi_447.5meV_VShift_100_lead_0_gammaLead_36.917.npz"
path_data_bands     = "./data/band_structure/97_band_eF_62.0meV_0.25BZ_Nkx_501.npz"

list_curent_name = list(path_data_currents.split('_'))
print("Fermi energy = ",list_curent_name[-7][:-3], " meV")


# In[45]:


# Load the data

def load_currents(path):
    data = np.load(path)
    total, up, down = data["Total"], data["Up"], data["Dn"]
    return total, up, down

total_current, up_current, down_current = load_currents(path_data_currents)
all_currents = [total_current, up_current, down_current]
color_maps   = ["Oranges", "Reds", "Blues"]
current_type = ["total","up","down"]

bands_total = np.load(path_data_bands)
kx, E_free, E_conf = bands_total['kx'], bands_total['E_free'], bands_total['E_conf']


# In[46]:


def make_system(esp="97", gammaLead=36.917, V_shift=100, width = shapes.W_STD, length = shapes.L_STD):
    shapeScattering = shapes.Rect(width, length)

    # folder_fig = esp + folder_suf

    params_raw = eval("gasb.params_" + esp)
    params_dict = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
    hamiltonian_syst = eval("gasb.hamiltonian_" + esp + "_k_plus()")

    hamiltonian_lead = gasb.free_ham(norbs = 6)
    sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, shapeScattering)

    return sistema


# In[47]:



# Formatação para os gráficos:
FONT_LABELS = 26
FONT_TITLES = 28
font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
matplotlib.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def formatter_current_axis(value, tick_number):
    #
    "Function for changing units of axis"
    return int(round(value * shapes.A0/10))

def edit_axis(axis, spin):
    "Editing axis"
    axis.set_ylim(-shapes.W_STD/2, shapes.W_STD/2)
    axis.xaxis.set_major_locator(plt.MultipleLocator(shapes.L_STD/2))
    axis.yaxis.set_major_locator(plt.MultipleLocator(shapes.W_STD/2))
    if spin.lower() == 'down':
        axis.xaxis.set_major_formatter(plt.FuncFormatter(formatter_current_axis))
        axis.yaxis.set_major_formatter(plt.FuncFormatter(formatter_current_axis))
        axis.set_xlabel(r'$x$ [nm]',fontsize=FONT_TITLES)
        axis.set_ylabel(r'$y$ [nm]',fontsize=FONT_TITLES, rotation=90)



# In[10]:


# (Re)Make the system
syst = make_system()


# In[79]:


L = 12
W = .98*L
F = plt.figure(1, (L, W))

ax_bands = plt.subplot(121)

grid_current = ImageGrid(F, 122,  # similar to subplot(122) row-column-num
                 nrows_ncols=(3, 1),
                 axes_pad=0.1,
                 add_all=True,
                 label_mode="1",
                 )

# Plot the bandstructure
fermi_energy = float(list_curent_name[-7][:-3])
trans.bands_cont2D_and_discr(axis = ax_bands,
                            free_elec_energies = E_free,
                            confined_elec_energies = E_conf,
                            momenta = kx,
                            kx_max = 0.12,
                            E_min = 435,
                            E_max = 455,
                            E_line = fermi_energy)


# Plot the maps for currents

axes_currents = [grid_current[0], grid_current[1], grid_current[2]]
labels = [r'(b)',r'(c)',r'(d)']

for current, colormap, name, label, axis in zip(all_currents, color_maps, current_type, labels, axes_currents):
    axis.yaxis.tick_right()
    axis.yaxis.set_label_position("right")
    kwant.plotter.current(syst, current, cmap = colormap, colorbar = False, show = False, ax=axis)

    axis.text(-0.01, 0.5, label, horizontalalignment='right',
         verticalalignment='center', transform=axis.transAxes, fontsize=FONT_TITLES)

    edit_axis(axis, name)
plt.subplots_adjust(wspace = 0.2)
# plt.tight_layout()
plt.show()
