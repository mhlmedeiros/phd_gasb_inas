#!/usr/bin/env python3

# Python
import sys
from itertools import product

# Anaconda
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# For eigenvalue computation
import scipy.sparse.linalg as sla

# Kwant
import kwant
import tinyarray
import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes


def plot_spectrum(syst, Efields, parameters):

    energies = []
    for eF in Efields:
        # Set the value of the electrical field
        parameters['eF'] = eF

        # Obtain the Hamiltonian as a sparse matrix
        ham_mat = syst.hamiltonian_submatrix(params=parameters, sparse=True)

        # we only calculate the 15 lowest eigenvalues
        ev = sla.eigsh(ham_mat.tocsc(), k=15, sigma=0,
                       return_eigenvectors=False)

        energies.append(ev)

    plt.figure()
    plt.plot(Efields, energies)
    plt.xlabel("electric field [meV]")
    plt.ylabel("energy [meV]")
    plt.show()


def calcula_energies(syst, parameters):

    energies = []

    # Obtain the Hamiltonian as a sparse matrix
    ham_mat = syst.hamiltonian_submatrix(params=parameters, sparse=True)
    N = ham_mat.shape[0]
    # we only calculate the 15 lowest eigenvalues
    ev = sla.eigsh(ham_mat.tocsc(), k=N//1000, sigma=440,
                    # which="LA",
                    return_eigenvectors=False)
    # TO DO: ordenar valores


    

    plt.figure()
    plt.plot(ev)
    plt.xlabel("N ")
    plt.ylabel("Energy [meV]")
    plt.show()

    np.savetxt("test_closed_system_BE.txt", ev)

def main():

    parameters  = gasb.params_97
    hamiltonian = gasb.hamiltonian_97_k_minus()
    syst_shape  = shapes.Rect(shapes.W_STD, shapes.L_STD)
    
    system = gasb.just_system_builder(hamiltonian, syst_shape)
    sub_matrix = system.hamiltonian_submatrix(params=parameters, sparse=True) 
    print(type(sub_matrix))
    print(sub_matrix.shape)


    eF_values = np.linspace(0,200,201)
    # plot_spectrum(system, eF_values, parameters)

    parameters["eF"] = 62.0
    calcula_energies(system, parameters)


if __name__ == '__main__':
    main()