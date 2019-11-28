#!/usr/bin/env python3

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






def main():
    """
    This code is for generate maps of current, without bands structures and
    with custom operators. This will allow a better exploratory process.
    """
    # Define the system
    hamiltonian = gasb.hamiltonian_97_k_plus()
    lead_ham = gasb.free_ham(norbs = 6)
    centralShape = shapes.Rect(Wmax = 300, Lmax=500)
    syst = gasb.system_builder(hamiltonian, lead_ham, centralShape, a_lattice=20)
    
    fig, axis = plt.subplots(1,1,figsize=(10,5))
    kwant.plot(syst, ax = axis, show = False)
    axis.axis('off')
    axis.set_aspect('equal')
    plt.show()


if __name__ == '__main__':
    main()
