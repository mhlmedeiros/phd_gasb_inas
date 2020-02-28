#!/usr/bin/env python3

# Python
import sys
from itertools import product

# Anaconda
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sla

# Kwant
import kwant
import tinyarray
import kwant.continuum

# This repository
from hamiltonians import gasb_hamiltonian as gasb
from transport_tools import bands_and_currents as trans
from system_geometry import shapes

def calc_spectrum(syst, Efields, parameters):

    energies = []
    n_calculado = 0
    for Ef in Efields:
        # Obtain the Hamiltonian as a sparse matrix
        n_calculado += 1
        parameters['eF'] = Ef
        ham_mat = syst.hamiltonian_submatrix(params=parameters, sparse=True)
        if n_calculado == 1:  print("Dimenção da matriz: ", ham_mat.shape)
        # print("Calculando: ", n_calculado,"/", len(Efields))


        # we only calculate the 15 lowest eigenvalues
        central_energy = 435
        ev = sla.eigsh(ham_mat.tocsc(),
                    k = 200,
                    sigma=central_energy,
                    return_eigenvectors = False)

        energies.append(ev)



    return energies

def main():
    # Define the system
    raio_simples = 10 * shapes.A_STD
    # raio_nm  = 150
    # raio_std = raio_nm * 10/shapes.A0
    print("QD radius: ", raio_simples, " Bohr radia.")

    # lattice_const_Bohr_radia = 200
    lattice_const_Bohr_radia = shapes.A_STD

    formato = shapes.Circular(raio_simples)
    H = gasb.hamiltonian_97_k_plus()
    quantum_dot = gasb.just_system_builder(H, formato, a_lattice=lattice_const_Bohr_radia)

    # Check that the system looks as intended.
    kwant.plot(quantum_dot)

    # Define parameters
    min_field     = 0
    max_field     = 1000
    N_pts         = 1000
    eF_values  = np.linspace(min_field, max_field, N_pts)
    parametros = gasb.params_97

    # Calculate and save energy levels of the dot
    energies = calc_spectrum(quantum_dot, eF_values, parametros)
    np.savez("circular_quantum_dot_energies_around_445_meV_field_0_62_and_200_eigenstates_R_150nm_a_200A0.npz", Ef_values=eF_values, Energies=energies)

    # and plot the results
    plt.figure()
    plt.plot(eF_values, energies,linestyle=' ', marker=',', color='k')
    plt.xlabel("E. field [meV]")
    plt.ylabel(r"$\varepsilon$ [meV]")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
