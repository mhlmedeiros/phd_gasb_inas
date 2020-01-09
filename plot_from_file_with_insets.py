import matplotlib.pyplot as plt
import matplotlib
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
import numpy as np

FONT_LABELS = 18
FONT_TITLES = 20
font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
matplotlib.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# Fermi = np.load("data_435_440_meV_Fermi_2001.npy")
# trans_num1 = np.load("data_435_440_meV_Fermi_Transport_UP_2001.npy")
# trans_num2 = np.load("data_435_440_meV_Fermi_Transport_DN_2001.npy")
# trans_total = trans_num1 + trans_num2

eF_total = np.load("data_62_64_eF_501_clean.npy")
trans_eF = np.load("data_62_64_transport_501_clean.npy")

eF_extrema = np.load("data_62_64_eF_501_clean_extrema.npy")
trans_eF_extrema = np.load("data_62_64_transport_501_clean_extrema.npy")

fig, ax1 = plt.subplots(figsize=(10,5))
ax1.plot(eF_total, trans_eF)
ax1.plot(eF_extrema, trans_eF_extrema, linestyle = ' ', marker = 'x')
# ax1.plot(Fermi, trans_total)
ax1.grid()
ax1.set_xlabel('eF [meV]')
# ax1.set_xlim((440, 435))
# ax1.set_xlim((440, 439.5))
# ax1.set_ylim((-2, 13))
ax1.set_ylabel(r'$\sigma$ $[e^2/h]$')
# ax1.legend()

# ax2 = plt.axes([0, 0, 1, 1])
# ip = InsetPosition(ax1, [0.5,0.07,0.3,0.3])
# ax2.set_axes_locator(ip)
# mark_inset(ax1, ax2, loc1=2, loc2=3, fc="none", ec='0.8')
#
# ax2.plot(Fermi[Fermi >= 439.6], trans_total[Fermi>=439.6])
# ax2.set_xlim((440, 439.6))
# ax2.tick_params(axis='both', which='major', labelsize=10)
# ax2.set_ylim((1.5, 3.0))
# ax2.grid()
#
#
# ax3 =  fig.add_axes([0.15, 0.65, 0.25, 0.25])
# ax3.plot(eF_total, trans_eF, color="red")
# ax3.tick_params(axis='both', which='major', labelsize=10)
# ax3.set_xlim((62,64))
# ax3.set_ylim((1.5, 3.0))
# ax3.set_xlabel(r"eF [meV]", fontsize = 10)
# ax3.grid()
plt.tight_layout()
plt.show()
