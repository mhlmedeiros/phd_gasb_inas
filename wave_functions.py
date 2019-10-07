'''
Script for analysis of wavefunctions on GaSb/InAs/GaSb simmetric quantum wells.

This piece code is part of the project "phd_gasb_inas", which comprises the work
related to the Phd. Dissertation named: "Quantum transport of charge and spin in
topological insulators 2D".

Author: Marcos Medeiros
email: mhlmedeiros@gmail.com
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Kwant related stuff
import kwant

# local application imports
from hamiltonians import gasb_hmailtonian as gasb
from system_geometry import shapes
