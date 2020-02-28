import matplotlib.pyplot as plt
import matplotlib
# from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
import numpy as np
import read_plot_current as read_plot

FONT_LABELS = 24
FONT_TITLES = 26
font = {'family' : 'serif', 'weight' : 'bold', 'size': FONT_LABELS}
matplotlib.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# Fermi = np.load("data_435_440_meV_Fermi_2001.npy")
# trans_num1 = np.load("data_435_440_meV_Fermi_Transport_UP_2001.npy")
# trans_num2 = np.load("data_435_440_meV_Fermi_Transport_DN_2001.npy")
# trans_total = trans_num1 + trans_num2

Fermi_total = np.load("./data/transport/data_439.7_440_meV_Fermi_1001.npy")
trans_eF = np.load("./data/transport/data_439.7_440_meV_Fermi_Transport_Total_1001.npy")

Fermi_extrema = np.load("./data/transport/data_439.7_440_meV_Fermi_1001_extrema_Fermi_439.7_440.npy")
trans_eF_extrema = np.load("./data/transport/data_439.7_440_meV_Fermi_Transport_Total_1001_extrema_Fermi_439.7_440.npy")

fig, ax1 = plt.subplots(figsize=(10,10))
ax1.plot(Fermi_total, trans_eF)
etiquetas_max = ['(a) 439.96','(b) 439.91','(c) 439.85','(d) 439.79','(e) 439.74']

for x, y, m, etiqueta in zip(Fermi_extrema[0::2], trans_eF_extrema[0::2], ["o","s","^","d","*"], etiquetas_max):
    ax1.plot(x, y, linestyle=' ', marker=m, markersize=10, label=etiqueta)


ax1.grid()
ax1.set_xlabel('Fermi [meV]')
ax1.set_xlim((440, 439.7))
ax1.set_ylim((0.5, 4))
ax1.set_ylabel(r'$G_{01}$ $[e^2/h]$')
ax1.legend(loc=4)

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
# ax1.arrow(eF_extrema[0], trans_eF_extrema[0]+0.02, -0.007, 0.21, head_width=0.005, head_length=0.02, fc='k', ec='k')

# ax3 =  fig.add_axes([0.11, 0.65, 0.25, 0.25])
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.96meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax3)
# # ax3.plot(eF_total, trans_eF, color="red")
# # ax3.tick_params(axis='both', which='major', labelsize=10)
# # ax3.set_xlim((62,64))
# # ax3.set_ylim((1.5, 3.0))
# # ax3.set_xlabel(r"eF [meV]", fontsize = 10)
# # ax3.grid()

title_font = 26
ax4 =  fig.add_axes([0.125, 0.70, 0.25, 0.25]) # for Min
read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.961meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax4) # Max
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.92meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax4) # Min
ax4.set_title("(a)",fontsize = title_font)


ax5 =  fig.add_axes([0.40, 0.70, 0.25, 0.25]) # for Min
read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.9073meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax5) # Max
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.875meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax5) # Min
ax5.set_title("(b)",fontsize = title_font)



ax6 =  fig.add_axes([0.675, 0.70, 0.25, 0.25]) # for Min
read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.8542meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax6) # Max
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.815meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax6) # Min
ax6.set_title("(c)",fontsize = title_font)


ax7 =  fig.add_axes([0.125, 0.10, 0.25, 0.25]) # for Min
read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.7948meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax7) # Max
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.765meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax7) # Min
ax7.set_title("(d)",fontsize = title_font)



ax8 =  fig.add_axes([0.40, 0.10, 0.25, 0.25]) # for Min
read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.7411meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax8) # Max
# read_plot.plot_map("./data/local_currents/97_currents_eF_62.0meV_Fermi_439.72meV_VShift_100_lead_0_gammaLead_36.917.npz", axis=ax8) # Min
ax8.set_title("(e)",fontsize = title_font)




plt.tight_layout()
plt.show()
