import kwant

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from system_geometry import shapes
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as tools


# For plotting:
font = {'family' : 'serif', 'weight' : 'bold', 'size': tools.FONT_LABELS}
matplotlib.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



# Remake the system:
def make_system(esp="97", gammaLead=36.917, V_shift=100, width = shapes.W_STD, length = shapes.L_STD):
    shapeScattering = shapes.Rect(width, length)

    # folder_fig = esp + folder_suf

    params_raw = eval("gasb.params_" + esp)
    params_dict = dict(GammaLead =  gammaLead, V = V_shift, **params_raw)
    hamiltonian_syst = eval("gasb.hamiltonian_" + esp + "_k_plus()")

    hamiltonian_lead = gasb.free_ham(norbs = 6)
    sistema          = gasb.system_builder(hamiltonian_syst, hamiltonian_lead, shapeScattering)

    return sistema


# Load the data from the file

def load_currents(path):
    data = np.load(path)
    total, up, down = data["Total"], data["Up"], data["Dn"] 
    return total, up, down

# Plot the current mapping
def plot_map(path,spin="Up",colormap="Reds",axis=None):
    data = np.load(path)

    syst = make_system()
    if axis == None:
        fig, axis = plt.subplots()
    
    kwant.plotter.current(syst, data[spin], cmap = colormap, colorbar = False, show = False, ax=axis)
    tools.edit_axis(axis,'total')
    axis.set_title(" ", fontsize=tools.FONT_TITLES)
    # plt.show()
    

 
def main():
    Fermi_energy = 439.8905
    # path_data = "./data/local_currents/97_currents_eF_62.0meV_Fermi_"+ str(Fermi_energy) +"meV_VShift_100_lead_0_gammaLead_36.917.npz"
    path_data = "./data/local_currents/97_currents_eF_62.0meV_Fermi_439.8905meV_VShift_100_lead_0_gammaLead_36.917_L_600.npz"
    total_current, up_current, down_current = load_currents(path_data)
    all_currents = [total_current, up_current, down_current]

    syst = make_system(length = 2*shapes.W_STD)
    
    for current, colormap, name in zip(all_currents, ["Oranges", "Reds", "Blues"],["total","up","down"]):
        fig, axis = plt.subplots()
        kwant.plotter.current(syst, current, cmap = colormap, colorbar = False, show = False, ax=axis)
        tools.edit_axis(axis,'total')
        axis.set_title(r"$\varepsilon = $ "+str(Fermi_energy)+" meV", fontsize=tools.FONT_TITLES)
        plt.savefig("../../"+name+"_current_"+str(Fermi_energy)+"L_600.png")


if __name__ =='__main__':
    main()