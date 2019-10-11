'''
    GaSb-Geometry (GaSb/InAs/GaSb) Hamiltonian definitions.

    Originaly it was named as "systemsGaSb.py" on "GaSb_InAs_project".

    Module to store functions to build systems adopted by kwant for transport
    calculations in GaSb-InAs-GaSb Quantum Wells.

    Hamiltonians and coeficients updated in July 24th of 2019;
    Test version of the Hamiltonian for the system with 103 Å typed in August
    the 6th of 2019.

'''
import kwant, kwant.continuum
import numpy as np

from system_geometry.shapes import *

"""
    # Default values for building the systems:
    #
    # Remembering that the Hamiltonian's parameters
    # adopt units which lengths are mesured in units
    # of Bohr's radius.

"""
# A0 = 0.529167 # raio de Bohr em Å = 10e-10 m
# L_STD  = 500    # nm = 10e-9
# W_STD  = 300    # nm = 10e-9
# A_STD  = 60     # units of Bohr's radius
# L_STD *= 10/A0  # convertion into units of Bohr's radius
# W_STD *= 10/A0  # convertion into units of Bohr's radius

# Funtions of electric field magnitude:
def ΔE_ΔR(x, A):
    # Pode ser usada para 'Delta E': sistemas 97 A
    # ou para 'Delta R': sistemas com 103 ou 110 A
    return A * x

def ΔγGeral(x, *args):
    if len(args) == 3:
        # Sistema com 97 A
        A, B, C = args
        test = False
    else:
        # Sistemas com 103 ou 110 A
        A, B, C, D, E, F = args
        test = True
    if (abs(x) >= 6) and test:
        # Só será executado para 103 A ou 110 A
        return D + E*x + F*x**2
    else:
        return A + B*np.exp(C*x)

def PxGeral(x, *args):
    if len(args) == 3:
        # 103 A (x é o campo)
        A, B, C = args
        return A + B*x + C*x**2
    else:
        # Para 97 e 110 A (x é um parâmetro)
        return x

def αCGeral(x, *args):
    if len(args) == 3:
        # 103 A ou 110 A
        A, B, C = args
        return A + B*x + C*x**2
    else:
        # 97 A
        A, B = args
        return A + B*x

def αVGeral(x, *args):
    # Essa função é igual para todos os sistemas
    A, B, C = args
    return A+B*np.exp(C*x)

def γCGeral(x, *args):
    if len(args) == 2:
        A, B = args
        return A + B*x**2
    else:
        return x

def γVGeral(x, *args):
    if len(args) == 3:
        A, B, C = args
        return A + B*x + C*x**2
    else:
        return x

# Hamiltonians for the systems :

def hamiltonian_97_k_plus():
    '''
    Hamiltonian given by Guilherme:

    This function produces the hamiltonian for the system with
    width of 97 Å.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x**2 + k_y**2) * GammaC + AlphaC(eF, A3, B3) * (k_x + k_y)',
    'H_12' : '+1j * (k_x + 1j * k_y) * Px(P)',
    'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_21' : '-1j * (k_x - 1j * k_y) * Px(P)',
    'H_22' : 'EV + DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_33' : 'EV - DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC - AlphaC(eF,A3,B3) * (k_x + k_y)',
    'H_45' : '+1j * (k_x - 1j * k_y) * Px(P)',
    'H_54' : '-1j * (k_x + 1j * k_y) * Px(P)',
    'H_55' : 'EV + DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_66' : 'EV - DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0,    0,    0, H_16],
       [H_21, H_22,    0,    0,    0,    0],
       [   0,    0, H_33, H_34,    0,    0],
       [   0,    0, H_43, H_44, H_45,    0],
       [   0,    0,    0, H_54, H_55,    0],
       [ H_61,   0,    0,    0,    0, H_66]
       ]
       """, locals = subs)

    return hamiltonian

def hamiltonian_97_k_minus():
    '''
    Hamiltonian given by Guilherme:

    This function produces the hamiltonian for the system with
    width of 97 Å.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x**2 + k_y**2) * GammaC + AlphaC(eF, A3, B3) * (k_x + k_y)',
    'H_12' : '+1j * (k_x - 1j * k_y) * Px(P)',
    'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_21' : '-1j * (k_x + 1j * k_y) * Px(P)',
    'H_22' : 'EV + DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_33' : 'EV - DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC - AlphaC(eF,A3,B3) * (k_x + k_y)',
    'H_45' : '+1j * (k_x + 1j * k_y) * Px(P)',
    'H_54' : '-1j * (k_x - 1j * k_y) * Px(P)',
    'H_55' : 'EV + DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_66' : 'EV - DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0,    0,    0, H_16],
       [H_21, H_22,    0,    0,    0,    0],
       [   0,    0, H_33, H_34,    0,    0],
       [   0,    0, H_43, H_44, H_45,    0],
       [   0,    0,    0, H_54, H_55,    0],
       [ H_61,   0,    0,    0,    0, H_66]
       ]
       """, locals = subs)

    return hamiltonian


def hamiltonian_97_up():
    '''
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 97 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x**2 + k_y**2) * GammaC + AlphaC(eF, A3, B3) * (k_x + k_y)',
    'H_12' : '+1j * (k_x + 1j * k_y) * Px(P)',
    'H_21' : '-1j * (k_x - 1j * k_y) * Px(P)',
    'H_22' : 'EV + DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_33' : 'EV - DeltaE(eF, A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) + AlphaV(eF,A4,B4,C4) * (k_x + k_y)'
    }

    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0],
       [H_21, H_22,    0],
       [   0,    0, H_33]
       ]
       """, locals = subs)

    return hamiltonian

def hamiltonian_97_down():

    """
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 97 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian."""

    sympify = kwant.continuum.sympify

    subs = {
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC - AlphaC(eF,A3,B3) * (k_x + k_y)',
    'H_45' : '+1j * (k_x - 1j * k_y) * Px(P)',
    'H_54' : '-1j * (k_x + 1j * k_y) * Px(P)',
    'H_55' : 'EV + DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV + DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)',
    'H_66' : 'EV - DeltaE(eF,A2) + (k_x**2 + k_y**2) * (GammaV - DeltaGamma(eF,A1,B1,C1)) - AlphaV(eF,A4,B4,C4) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_44, H_45,    0],
       [H_54, H_55,    0],
       [   0,    0, H_66]
       ]
       """, locals = subs)

    return hamiltonian


def hamiltonian_103_k_plus():
    '''
    Hamiltonian given by Guilherme with modifications:

    H_12 : k_x  -->  k_x + 1j * k_y
    H_21 = conj(H_12)

    H_45(k_x,k_y) = conj(H_12(-k_x, -k_y))
    H_54 = conj(H_45)

    I sincerely believed that this must be the correct formulation,
    due to its agreement with Krishtopenko formulation.Besides that,
    this Hamiltonian produces edges states in differents egdes of
    the sample.

    This function produces the hamiltonian for the system with
    width of 103 Å.
    '''
    sympify = kwant.continuum.sympify

    subs = {
        'H_11' : 'EC + (k_x + k_y) * AlphaC(eF, A6, B6, C6) + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) ',
        'H_12' : '+1j * (k_x + 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
        'H_21' : '-1j * (k_x - 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_23' : '+1j * DeltaR(eF, A1)',
        'H_32' : '-1j * DeltaR(eF, A1)',
        'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
        'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
        'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A6, B6, C6) * (k_x + k_y)',
        'H_45' : '+1j * (k_x - 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_54' : '-1j * (k_x + 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_56' : '-1j * DeltaR(eF, A1)',
        'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
        'H_65' : '+1j * DeltaR(eF, A1)',
        'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
           13.6 * 1000 *[
           [H_11, H_12,    0,    0,    0, H_16],
           [H_21, H_22, H_23,    0,    0,    0],
           [   0, H_32, H_33, H_34,    0,    0],
           [   0,    0, H_43, H_44, H_45,    0],
           [   0,    0,    0, H_54, H_55, H_56],
           [H_61,    0,    0,    0, H_65, H_66]
           ]
    """, locals = subs)

    return hamiltonian

def hamiltonian_103_k_minus():
    '''
    Hamiltonian given by Guilherme with some modifications:

    H_12 : k_x  -->  k_x + 1j * k_y
    H_21 = conj(H_12)

    H_45 : k_x  -->  k_x + 1j * k_y
    H_54 = conj(H_45)

    This function produces the hamiltonian for the system with
    width of 103 Å.
    '''
    sympify = kwant.continuum.sympify
    subs = {
        'H_11' : 'EC + (k_x + k_y) * AlphaC(eF, A6, B6, C6) + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) ',
        'H_12' : '+1j * (k_x - 1j*k_y) * Px(eF, A5, B5, C5)',
        'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
        'H_21' : '-1j * (k_x + 1j*k_y) * Px(eF, A5, B5, C5)',
        'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_23' : '+1j * DeltaR(eF, A1)',
        'H_32' : '-1j * DeltaR(eF, A1)',
        'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
        'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
        'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A6, B6, C6) * (k_x + k_y)',
        'H_45' : '+1j * (k_x + 1j*k_y) * Px(eF, A5, B5, C5)',
        'H_54' : '-1j * (k_x - 1j*k_y) * Px(eF, A5, B5, C5)',
        'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_56' : '-1j * DeltaR(eF, A1)',
        'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
        'H_65' : '+1j * DeltaR(eF, A1)',
        'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)'
        }
    hamiltonian = sympify("""
           13.6 * 1000 *[
           [H_11, H_12,    0,    0,    0, H_16],
           [H_21, H_22, H_23,    0,    0,    0],
           [   0, H_32, H_33, H_34,    0,    0],
           [   0,    0, H_43, H_44, H_45,    0],
           [   0,    0,    0, H_54, H_55, H_56],
           [H_61,    0,    0,    0, H_65, H_66]
           ]
    """, locals = subs)

    return hamiltonian


def hamiltonian_103_up():
    '''
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 103 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian.
    '''
    sympify = kwant.continuum.sympify

    subs = {
        'H_11' : 'EC + (k_x + k_y) * AlphaC(eF, A6, B6, C6) + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) ',
        'H_12' : '+1j * (k_x + 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_21' : '-1j * (k_x - 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_23' : '+1j * DeltaR(eF, A1)',
        'H_32' : '-1j * DeltaR(eF, A1)',
        'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A7, B7, C7) * (k_x + k_y)'
        }

    hamiltonian = sympify("""
           13.6 * 1000 *[
           [H_11, H_12,    0],
           [H_21, H_22, H_23],
           [   0, H_32, H_33]
           ]
    """, locals = subs)

    return hamiltonian

def hamiltonian_103_down():

    """
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 103 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian."""

    sympify = kwant.continuum.sympify

    subs = {
        'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A6, B6, C6) * (k_x + k_y)',
        'H_45' : '+1j * (k_x - 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_54' : '-1j * (k_x + 1j * k_y) * Px(eF, A5, B5, C5)',
        'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)',
        'H_56' : '-1j * DeltaR(eF, A1)',
        'H_65' : '+1j * DeltaR(eF, A1)',
        'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A7, B7, C7) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
           13.6 * 1000 *[
           [H_44, H_45,    0],
           [H_54, H_55, H_56],
           [   0, H_65, H_66]
           ]
    """, locals = subs)

    return hamiltonian


def hamiltonian_110_k_plus():
    '''
    Hamiltonian given by Guilherme with modifications:

    H_12 : k_x  -->  k_x + 1j * k_y
    H_21 = conj(H_12)

    H_45(k_x,k_y) = conj(H_12(-k_x, -k_y))
    H_54 = conj(H_45)

    I sincerely believed that this must be the correct formulation,
    due to its agreement with Krishtopenko formulation. Besides that,
    this Hamiltonian produces edges states in differents egdes of
    the sample.

    This function produces the hamiltonian for the system with
    width of 110 Å.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) + AlphaC(eF, A5, B5, C5) * (k_x + k_y)',
    'H_12' : '+1j * (k_x + 1j * k_y) * Px(P)',
    'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_21' : '-1j * (k_x - 1j * k_y) * Px(P)',
    'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_23' : '+1j * DeltaR(eF, A1)',
    'H_32' : '-1j * DeltaR(eF, A1)',
    'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A5, B5, C5) * (k_x + k_y)',
    'H_45' : '+1j * (k_x - 1j * k_y) * Px(P)',
    'H_54' : '-1j * (k_x + 1j * k_y) * Px(P)',
    'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_56' : '-1j * DeltaR(eF, A1)',
    'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_65' : '+1j * DeltaR(eF, A1)',
    'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0,    0,    0, H_16],
       [H_21, H_22, H_23,    0,    0,    0],
       [   0, H_32, H_33, H_34,    0,    0],
       [   0,    0, H_43, H_44, H_45,    0],
       [   0,    0,    0, H_54, H_55, H_56],
       [ H_61,   0,    0,    0, H_65, H_66]
       ]
       """, locals = subs)

    return hamiltonian

def hamiltonian_110_k_minus():
    '''
    Hamiltonian given by Guilherme with modifications:

    H_12 : k_x  -->  k_x + 1j * k_y
    H_21 = conj(H_12)

    H_45 : k_x  -->  k_x + 1j * k_y
    H_54 = conj(H_45)

    This function produces the hamiltonian for the system with
    width of 110 Å.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x + k_y) * AlphaC(eF, A5, B5, C5) + (k_x**2 + k_y**2) * GammaC(eF, A2, B2)',
    'H_12' : '+1j * (k_x - 1j*k_y) * Px(P)',
    'H_16' : '-(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_21' : '-1j * (k_x + 1j*k_y) * Px(P)',
    'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_23' : '+1j * DeltaR(eF, A1)',
    'H_32' : '-1j * DeltaR(eF, A1)',
    'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_34' : '(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_43' : '(k_x**2 - k_y**2) * Eta2 - 1j * k_x * k_y * Eta3',
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A5, B5, C5) * (k_x + k_y)',
    'H_45' : '+1j * (k_x + 1j*k_y) * Px(P)',
    'H_54' : '-1j * (k_x - 1j*k_y) * Px(P)',
    'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_56' : '-1j * DeltaR(eF, A1)',
    'H_61' : '-(k_x**2 - k_y**2) * Eta2 + 1j * k_x * k_y * Eta3',
    'H_65' : '+1j * DeltaR(eF, A1)',
    'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0,    0,    0, H_16],
       [H_21, H_22, H_23,    0,    0,    0],
       [   0, H_32, H_33, H_34,    0,    0],
       [   0,    0, H_43, H_44, H_45,    0],
       [   0,    0,    0, H_54, H_55, H_56],
       [ H_61,   0,    0,    0, H_65, H_66]
       ]
       """, locals = subs)

    return hamiltonian


def hamiltonian_110_up():
    '''
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 110 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian.
    '''
    sympify = kwant.continuum.sympify

    subs = {
    'H_11' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) + AlphaC(eF, A5, B5, C5) * (k_x + k_y)',
    'H_12' : '+1j * (k_x + 1j * k_y) * Px(P)',
    'H_21' : '-1j * (k_x - 1j * k_y) * Px(P)',
    'H_22' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_23' : '+1j * DeltaR(eF, A1)',
    'H_32' : '-1j * DeltaR(eF, A1)',
    'H_33' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) + AlphaV(eF, A6, B6, C6) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_11, H_12,    0],
       [H_21, H_22, H_23],
       [   0, H_32, H_33]
       ]
       """, locals = subs)

    return hamiltonian

def hamiltonian_110_down():

    """
    Hamiltonian given by Guilherme:

    This function produces the first block of the hamiltonian for the system
    with width of 110 Å. Note that we've adopted Eta3 == Eta2 == 0 here in order to
    uncouple the superior and inferior blocks of the original Hamiltonian.

    We've also considered the version "plus" of the Hamiltonian."""

    sympify = kwant.continuum.sympify

    subs = {
    'H_44' : 'EC + (k_x**2 + k_y**2) * GammaC(eF, A2, B2) - AlphaC(eF, A5, B5, C5) * (k_x + k_y)',
    'H_45' : '+1j * (k_x - 1j * k_y) * Px(P)',
    'H_54' : '-1j * (k_x + 1j * k_y) * Px(P)',
    'H_55' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) + DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)',
    'H_56' : '-1j * DeltaR(eF, A1)',
    'H_65' : '+1j * DeltaR(eF, A1)',
    'H_66' : 'EV + (k_x**2 + k_y**2) * (GammaV(eF, A3, B3, C3) - DeltaGamma(eF,A4,B4,C4,D4,E4,F4)) - AlphaV(eF, A6, B6, C6) * (k_x + k_y)'
    }



    hamiltonian = sympify("""
       13.6 * 1000 *[
       [H_44, H_45,    0],
       [H_54, H_55, H_56],
       [   0, H_65, H_66]
       ]
       """, locals = subs)

    return hamiltonian


# Builder:
def system_builder(hamiltonian, centralShape, a = A_STD):
    '''
    Here we define the system's geometry as a rectangular slab attached with leads
    of same width.

    The hamiltonian argument is supposed to be an output of one of the following
    function:
        * hamiltonian_97
        * hamiltonian_103
        * hamiltonian_110

    The subsequent arguments are:
        a -> lattice parameter: separation between sites
        L -> length of the scattering region
        W -> width of the leads and for scattering region

    All of these last parameter have default values
        A_STD   = 60        # units of Bohr's radius
        L_nano  = 500       # nanometers
        W_nano  = 300       # nanometers

        L_STD   = L_nano * 10/ A0   # units of Bohr's radius
        W_STD   = W_nano * 10/ A0   # units of Bohr's radius

    '''
    template = kwant.continuum.discretize(hamiltonian, grid=a)

    # def shape(site):
    #     (x, y) = site.pos
    #     return (-W/2 < y < W/2 and -L/2 < x < L/2)

    def lead_shape(site):
        (x, y) = site.pos
        return (-centralShape.Wmax/2 < y < centralShape.Wmax/2)

    syst = kwant.Builder()
    syst.fill(template, centralShape, (0, 0))

    lead = kwant.Builder(kwant.TranslationalSymmetry([-a, 0]))
    lead.fill(template, lead_shape, (-centralShape.Lmax, 0))

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    syst = syst.finalized()
    return syst

def just_lead_builder(hamiltonian, W = W_STD, a = A_STD, symmetry=-1):

    template = kwant.continuum.discretize(hamiltonian, grid=a)

    def lead_shape(site):
        (x, y) = site.pos
        return (-W/2 < y < W/2)

    lead = kwant.Builder(kwant.TranslationalSymmetry([symmetry * a, 0]))
    lead.fill(template, lead_shape, (0, 0))

    return lead.finalized()


params_97 = {
    "EC" : 0.0320879,
    "EV" : 0.0314029,
    "GammaC" : 36.917,
    "GammaV" : -22.4782,
    "P" : -0.108812,
    "Eta2" : -0.280174,
    "Eta3" : -1.16087,
    "A1": -0.231075,
    "B1": 4.11407,
    "C1":-0.255277,
    "A2": 2.79538e-5,
    "A3": 3.24001e-5,
    "B3": 8.28089e-5,
    "A4": 0.0219519,
    "B4": -0.02344,
    "C4": -0.286293,
    "DeltaE": ΔE_ΔR,
    "DeltaGamma": ΔγGeral,
    "AlphaC": αCGeral,
    "AlphaV": αVGeral,
    "Px": PxGeral,
    "eF": 0
}

params_103 = {
    "EC": 0.0313913,
    "EV": 0.0314028,
    "Eta2": -1.37783,
    "Eta3": 4.51996,
    "A1": -0.0000272091,
    "A2": 42.6859, "B2": -0.00694795,
    "A3": -19.1031, "B3": -0.0203693, "C3": 0.00417352,
    "A4": 1.31607, "B4": -3.7044, "C4": -0.140137,
    "D4": -0.158825, "E4": -0.0213621, "F4": -0.0031922,
    "A5": 0.0971718, "B5": 1.7474e-7, "C5": -5.92765e-7,
    "A6": -0.0000541767, "B6": 0.000248155, "C6": -3.66946e-7,
    "A7": 0.0164299, "B7": -0.0166324, "C7": -0.43601,
    "GammaC": γCGeral,
    "GammaV": γVGeral,
    "DeltaR": ΔE_ΔR,
    "DeltaGamma": ΔγGeral,
    "AlphaC": αCGeral,
    "AlphaV": αVGeral,
    "Px": PxGeral,
    "eF": 0
}

params_110 = {
    "EC": 0.030775,
    "EV": 0.0314028,
    "Eta2": 0.0866876,
    "Eta3": 0.398958,
    "P": 0.125543,
    "A1": -0.0000277126,
    "A2": 39.3037, "B2": -0.00892854,
    "A3": -16.9182, "B3": 0.0434725, "C3": 0.00284609,
    "A4": -0.0609416, "B4": -0.560249, "C4": -0.536596,
    "D4": 0.494559, "E4": -0.156445, "F4": -0.00221618,
    "A5": 0.0000647954, "B5": 0.000407089, "C5": -4.18461e-6,
    "A6": 0.0206344, "B6": -0.0210344, "C6": -0.270751,
    "GammaC": γCGeral,
    "GammaV": γVGeral,
    "DeltaR": ΔE_ΔR,
    "DeltaGamma": ΔγGeral,
    "AlphaC": αCGeral,
    "AlphaV": αVGeral,
    "Px": PxGeral,
    "eF": 0
}
