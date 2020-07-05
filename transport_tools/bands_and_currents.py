import kwant
import scipy
from sympy import init_printing
from sympy import *
import tinyarray
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


"""
    # The following import statement is crucial
    # without it we'd lost all the important
    # constant definitions
"""

from system_geometry import shapes
from hamiltonians.gasb_hamiltonian import *


# Formatação para os gráficos:
FONT_LABELS = 30
FONT_TITLES = 32
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
        axis.set_xlabel(r'$x$ [nm]',fontsize=FONT_TITLES)
        axis.set_ylabel(r'$y$ [nm]',fontsize=FONT_TITLES)

def trans_momenta(k_x):
    return k_x * (shapes.A_STD*shapes.A0*10**(-1))**(-1) # conversion from nm^{-1}

def current_density(axis, syst, parameters, eF_value = 0, energy = 428, lead_index=0, colormap="Reds"):
    parameters["eF"] = eF_value

    wf = kwant.wave_function(syst, energy=energy, params=parameters)
    J_spin = kwant.operator.Current(syst)

    current_spin = sum(J_spin(psi, params = parameters) for psi in wf(lead_index))
    # kwant.plotter.current(syst, current_spin, cmap = colormap, colorbar = False, show = False, ax=axis, density=1/9)
    kwant.plotter.current(syst, current_spin, cmap = colormap, colorbar = False, show = False, ax=axis)
    edit_axis(axis, "none") # change units to nm
    axis.set_title(" ")

    return 0

def current_spin(syst, parameters, eF_value = 0, energy = 428, lead_index=0, spin="up", colormap="Reds", axis=None):

    """
    Gera mapa da densidade de corrente para elétrons incidentes
    da lead esquerda. Se o axis é dado a função gera uma figura
    caso contrário somente a array resultante é retornada. Para
    se gerar uma figura posterioriormente, o systema deve ser
    construído de modo identico àquele usado para calcular a
    corrente (syst).

    "syst": Note que o sistema é TOTAL com a matriz Hamiltoniana 6x6;
    "parameters": dicionario -> define const. para a Hamiltoniana;
    "eF_value" : valor da magnitude do campo elétrico em meV;
    "energy": potencial eletroquímico dos eletrons incidentes;
    "lead_index": "0" para elétrons incidentes da esquerda e "1" para elétrons da direita;
    "spin": projeção do pseudo spin
    "colormap" = mapa de cores usado para 'plotting'
    "axis" = eixo usado para plot, caso seja 'None' a função retorna a array com os valores de corrente mas não gera figura

    """

    parameters["eF"] = eF_value
    if spin.lower() == "total" : up, down = 1, 1
    elif spin.lower() == "up"  : up, down = 1, 0
    elif spin.lower() == "down": up, down = 0, 1

    sz_sub_matrix   = np.array([[up, 0],[0, down]])
    Sz_total_matrix = np.kron(sz_sub_matrix, np.eye(3))
    wf = kwant.wave_function(syst, energy=energy, params=parameters)
    J_spin = kwant.operator.Current(syst, Sz_total_matrix)

    current_spin = sum(J_spin(psi, params = parameters) for psi in wf(lead_index))
    # kwant.plotter.current(syst, current_spin, cmap = colormap, colorbar = False, show = False, ax=axis, density=1/9)

    if axis != None:
        kwant.plotter.current(syst, current_spin, cmap = colormap, colorbar = False, show = False, ax=axis)
        edit_axis(axis, spin) # change units to nm

    return current_spin

def continuous_bands_2D(kx_array, ky_value, hamiltonian, params, eF_value = 0):
    params["eF"] = eF_value
    H_k = hamiltonian.subs(params)
    H_k = kwant.continuum.lambdify(H_k)

    if len(params) > 25:
        Bands_cont_2D = np.array([scipy.linalg.eigvalsh(H_k(k_x = kx, k_y = ky_value,
            AlphaC = αCGeral, AlphaV = αVGeral,
            DeltaGamma = ΔγGeral, DeltaR = ΔE_ΔR,
            GammaC = γCGeral, GammaV = γVGeral,
            Px = PxGeral)) for kx in kx_array])
    else:
        Bands_cont_2D = np.array([scipy.linalg.eigvalsh(H_k(k_x = kx, k_y = ky_value,
            AlphaC = αCGeral, AlphaV = αVGeral,
            DeltaGamma = ΔγGeral, DeltaE = ΔE_ΔR,
            Px = PxGeral)) for kx in kx_array])
    return Bands_cont_2D

def plot_bands_cont2D(axis, H_func_symb, kx_array, ky_value, params):
    '''
    This function generates a simple 2D plot for the
    six energy bands of the passed hamiltonian. Such hamiltonian
    must be k_x and k_y dependent, following the kwant convention.

    Input:
        * Function that generates symbolic Hamiltonian:
            - "hamiltonian_InAs"; or
            - "hamiltonian_GaSb"
        * values for "kx": np.array
        * one single value for "ky": float or int
        * parameters for the hamiltonian: dict
    '''

    H_symb = H_func_symb()
    Bands = continuous_bands_2D(kx_array, ky_value,
                             hamiltonian = H_symb,
                             params = params)

    # fig = plt.figure(figsize=(6,6))
    # ax1 = fig.add_subplot(111)
    ax1.plot(kx_array, Bands, linewidth=2.,)
    plt.tight_layout()
    plt.show()

def band_values(ham_syst, momenta, params, eF_value = 0):

    """
    Calculate the energy bands and return it as a matrix.

        * syst      : kwant finalized system
        * momenta   : array with values for kx in units of shapes.A_STD^{-1}
        * params    : dictionary with parameters for Hamiltonian
        * eF        : value for electric field perpendicualr to the QW
    """
    lead = lead_metalica(ham_syst).finalized()
    params['eF'] = eF_value
    bandas = kwant.physics.Bands(lead, params = params)

    energies = [bandas(k) for k in momenta]

    return energies

def band_with_line_gasb(axis, momenta, energies,
                    kx_max = 0.25, E_min = 300, E_max = 500, E_line = None,
                    color_plot = "black", linestyle_plot = "-", marker_plot = None,
                    color_line = "red", linestyle_line = "--",label_plot = None):

    """
    Função para plot da estrutura de bandas
        * energies  : matrix com valores das energias
        * momenta   : array com os valores de kx
        * kx_max    : valor max. para plot de "|k_x|"
        * E_min     : valor min. para plot de energia
        * E_max     : valor max. para plot de energia
        * E_line    : Nível de energia de interesse (traçar linha)

    Output: Gráfico da estrutura de bandas dependente do
    momentum na direção-x. A unidade de energia depende
    das unidades dos parâmetros adotados no Hamiltoniano enquanto
    que a unidade de momentum é o inverso do parâmetro de rede "a",
    esse, por sua vez, tem unidade também dependente dos parâmetros
    do Hamiltoniano.
    """

    # font_labels = 22
    # font_titles = 22
    font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
    matplotlib.rc('font', **font)


    momenta_trans = momenta * (shapes.A_STD*shapes.A0*10**(-1))**(-1) # conversion to nm^{-1}

    # fig = plt.figure()
    # axis = fig.add_subplot(111)
    energies = np.array(energies)

    # The first first sub-band is labeled and it will represent the bandstructure
    axis.plot(momenta_trans, energies[:, 0],
                linewidth = 1.0, color = color_plot, linestyle = linestyle_plot,
                marker = marker_plot, markevery=20, label = label_plot)
    axis.plot(momenta_trans, energies[:, 1:],
                linewidth = 1.0, color = color_plot, linestyle = linestyle_plot,
                marker = marker_plot, markevery=20)

    axis.hlines(E_line, -kx_max, kx_max,
                linewidth = 2.0, color = color_line, linestyle = linestyle_line)
    axis.grid()
    axis.set_xlim(-kx_max, kx_max)
    axis.set_ylim(E_min, E_max)
    axis.set_xlabel(r'k$_x$ [$nm^{-1}$]', fontsize=FONT_TITLES)
    axis.set_ylabel(r'$\varepsilon$ [meV]', fontsize=FONT_TITLES)
    # fig = plt.gcf()
    # fig.set_size_inches(6, 6)
    # plt.tight_layout()
    # plt.show()
    return 0

def bands_cont2D_and_discr(axis, free_elec_energies, confined_elec_energies, momenta, kx_max=0.25, E_min=300, E_max = 500, E_line = None):

    momenta_trans = momenta * (shapes.A_STD * shapes.A0*10**(-1))**(-1) # conversion to nm^{-1}
    axis.plot(momenta_trans, free_elec_energies, linewidth=2.5,color="blue")
    axis.plot(momenta_trans, confined_elec_energies, linewidth = 0.8, color="black", linestyle="-")
    axis.hlines(E_line, -kx_max, kx_max, linewidth = 2.0, linestyle = "--", color = "red")
    axis.grid()
    axis.set_title("(a)", fontsize=FONT_TITLES)
    axis.set_xlim(-kx_max, kx_max)
    axis.set_ylim(E_min, E_max)
    axis.set_xlabel(r'k$_x$ [$nm^{-1}$]', fontsize=FONT_TITLES)
    axis.set_ylabel(r'$\varepsilon$ [meV]', fontsize=FONT_TITLES)
    # plt.tight_layout()
    # plt.show()
    return 0

def continuous_levels_eF(kx_value, ky_value, hamiltonian, params, eF_array):
    params_witout_eF = params.copy()
    params_witout_eF.pop("eF")
    H_k = hamiltonian.subs(params_witout_eF)
    H_k = kwant.continuum.lambdify(H_k)

    if len(params) > 25 : # len(params_97) = 22, len(params_103) = 33, len(params_110) = 31
        levels_cont = np.array([scipy.linalg.eigvalsh(H_k(k_x = kx_value, k_y = ky_value,
            AlphaC = αCGeral, AlphaV = αVGeral,
            DeltaGamma = ΔγGeral, DeltaR = ΔE_ΔR,
            GammaC = γCGeral, GammaV = γVGeral,
            Px = PxGeral, eF=eField_value)) for eField_value in eF_array])
    else:
        levels_cont = np.array([scipy.linalg.eigvalsh(H_k(k_x = kx_value, k_y = ky_value,
            AlphaC = αCGeral, AlphaV = αVGeral,
            DeltaGamma = ΔγGeral, DeltaE = ΔE_ΔR,
            Px = PxGeral, eF=eField_value)) for eField_value in eF_array])
    return levels_cont

def kx_zero_levels(syst_lead_like, params):
    '''
    This function GETS a
        * translational symmetric finalized system: syst_lead_like; and a
        * dictionary of parameters: params;
    and RETURNS all the energy levels for k_x = 0.
    '''
    bands = kwant.physics.Bands(syst_lead_like, params = params)
    energies = bands(0) # This object is callable: given a momentum
                        # it returns a array containing the eigenenergies
                        # ao all modes at this momentum.
    return energies

def collected_energies(syst_lead_like, params, eF_array):
    '''
    This function generates a matrix containing in each column the eigenenergies
    for each value for electric field present in the array eF_array.

    We must pass to this function
        * the translational symmetric finalized system: syst_lead_like;
        * the dictionary of parameters: params; and
        * the array containing the desired values of electric field
    '''
    energies_collected = []

    for eField in eF_array:
        params["eF"] = eField
        energies_collected.append(kx_zero_levels(syst_lead_like, params))

    return energies_collected
