import numpy as np
import matplotlib.pyplot as plt

A0 = 0.529167 # raio de Bohr em Å = 10e-10 m
L_STD  = 500    # nm = 10e-9
W_STD  = 300    # nm = 10e-9
R_STD  = 200    # nm = 10e-9
A_STD  = 60     # units of Bohr's radius
L_STD *= 10/A0  # convertion into units of Bohr's radius
W_STD *= 10/A0  # convertion into units of Bohr's radius
R_STD *= 10/A0  # convertion into units of Bohr's radius

# classes (and function) for constricted shaped system
class SmoothRect:

    def __init__(self, Vmax, beta, x0, L):
        self.Vmax = Vmax
        self.beta = beta
        self.x0   = x0
        self.L    = L


    def __call__(self, x):
        V    = self.Vmax
        beta = self.beta
        x0   = self.x0
        L    = self.L
        return V * (expFermi(x, x0,-beta) + expFermi(x, x0+L, beta))

class ConstrictSmoothRect:

    def __init__(self, Wmax = W_STD, Wmin = 0.5 * W_STD,
                beta = A_STD, Lmax = L_STD, Lconst = 0.5 * L_STD):
        self.Wmax   = Wmax
        self.Wmin   = Wmin
        self.beta   = beta
        self.Lmax   = Lmax
        self.Lconst = Lconst

    def __call__(self, Site):
        (x, y) = Site.pos
        V    = (self.Wmax - self.Wmin)/2
        beta = self.beta
        x0   = -(self.Lconst/2)
        l    = self.Lconst
        F = SmoothRect( V, beta, x0, l)
        G = SmoothRect(-V, beta, x0, l)
        f = F(x) + self.Wmin/2 # limite superior
        g = G(x) - self.Wmin/2 # limite inferior
        # f =   self.Wmax/2
        # g = - self.Wmax/2
        return (g < y < f) and (-self.Lmax/2 < x < self.Lmax/2)

class Rect:

    def __init__(self, Wmax = W_STD, Lmax = L_STD):
        self.Wmax = Wmax
        self.Lmax = Lmax

    def __call__(self, Site):
        (x, y) = Site.pos
        W = self.Wmax
        L = self.Lmax
        return (-L/2 < x < L/2) and (-W/2 < y < W/2)

class Circular:

    def __init__(self, Radius = R_STD):
        self.Radius = Radius
    
    def __call__(self, Site):
        (x, y) = Site.pos
        R = self.Radius
        return x**2 + y**2 < R**2 

def expFermi(x, x0, beta):
    return 1/(np.exp(beta * (x0-x)) + 1)
