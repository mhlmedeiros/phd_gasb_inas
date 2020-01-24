import kwant, kwant.continuum
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Hamiltonian
def hamiltonian_3D_mag_field():
    
    H_subs = {"K": "gamma * (k_x**2 + k_y**2 + (k_z - alpha * (y*k_y + x*k_x)*E)**2)",
            "V" : "V",
            "Z" : "g_eff * E * k_y" # For k_x == 0
            }
    
    sympify = kwant.continuum.sympify
    H = sympify("K + V", locals = H_subs)

    return H

# System
def lead_syst(hamiltonian, W = 100, T = 10, a_lattice = 1):

    template = kwant.continuum.discretize(hamiltonian, grid = a_lattice)

    def lead_shape(site):
        (x, y, z) = site.pos
        return (-W/2 < y < W/2) and (-T/2 < z < T/2)

    lead = kwant.Builder(kwant.TranslationalSymmetry([-a_lattice,0,0]))
    lead.fill(template, lead_shape, (-a_lattice,0,0))

    return lead


def system_builder(hamiltonian, L_syst = 100, W_syst = 50, T_syst = 10, a_syst = 1):

    template = kwant.continuum.discretize(hamiltonian, grid = a_syst)

    def central_shape(site):
        (x, y, z) = site.pos
        return (-L_syst/2 < x < L_syst/2) and (-W_syst/2 < y < W_syst/2) and (-T_syst/2 < z < T_syst/2)

    syst = kwant.Builder()
    syst.fill(template, central_shape, (0, 0, 0))

    lead = lead_syst(hamiltonian, W = W_syst, a_lattice= a_syst)

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    syst = syst.finalized()

    return syst

# Energy levels



def main():
    H = hamiltonian_3D_mag_field()
    
    lead = lead_syst(H, W=50, T=10, a_lattice=1)

    # syst = system_builder(H)
    # kwant.plot(syst)
    
    # print(type(lead.finalized()))

    parameters = dict(
        gamma = 1,
        alpha = 1,
        E     = 0,
        V     = 0,
        g_eff = 1
    )

    # bands = kwant.physics.Bands(lead.finalized(), params=parameters)
    # momenta = np.linspace(-np.pi,np.pi,101)
    # energies = [bands(k) for k in momenta]
    # print(type(energies[0]))
    # print(len(energies[0]))


    E_field_values = np.linspace(0,100,2001)
    energies = []
    for E_value in E_field_values:
        parameters['E'] = E_value
        bands = kwant.physics.Bands(lead.finalized(), params=parameters)
        # momenta = np.linspace(-np.pi,np.pi,101)
        # energies = [bands(k) for k in momenta]
        energies.append(bands(0))

    name_file = "./data/3D_effective_system/levels_at_kx_zero/gamma_1_alpha_1_V_0_g_eff_1_E_0_100_2001_pts.npz"
    np.savez(name_file, E_field=E_field_values, energy_levels=energies)
    
    plt.plot(E_field_values, energies)
    plt.show()

if __name__ == "__main__":
    main()