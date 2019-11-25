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
