# Credits: A. Velberg and M. Calvo Carrera

import numpy as np
import matplotlib.pyplot as plt
# BoP implementation
# plt.rcParams.update({
#     "text.usetex": True
# })

#Pfus = 250e6 # W
def calcQe(Pfus):
    n_factor = 1.11  # mUltiplication factor
    Pthermal = Pfus * n_factor

    e_hx = 0.95  # Efficiency of heat exchangers between loops
    Pheat = Pthermal * e_hx ** 2

    e_rankine = .333  # ideal rankine efficiency
    # Pe_total = Pheat*e_rankine

    Paux = 10e6  # W
    e_aux = 0.5  # Aux power efficiency
    Pcryo = 2 * 200e3  # kW
    e_cryo = 0.025
    Ppump_Flibe = 2e5  # Vdot*rho*g*h, Vdot = 1
    e_pump_Flibe = 0.75

    h1 = 3540e3  # J/kgK
    h4 = 670.38e3  # J/kgK
    h3_liq = 292e3  # J/kgK
    h3_vap = 2626e3  # J/kgK
    x3 = 0.13
    x2 = 0.95
    h3 = (1 - x3) * h3_liq + x3 * h3_vap
    h2 = (1 - x2) * h3_liq + x2 * h3_vap
    mdot = Pheat / (h1 - h4)  # kg/s
    Vdot = mdot / 1000
    Ppump_water = mdot * (h4 - h3)
    e_pump_water = 0.75

    e_rankine_real = (h1 - h2) / (h1 - h4)
    Pe_total = Pheat * e_rankine_real

    Pe_systems = Paux / e_aux + Pcryo / e_cryo + 2 * Ppump_Flibe / e_pump_Flibe + Ppump_water / e_pump_water
    Pe_avail = Pe_total - Pe_systems

    Qe_rankine = Pe_total / (Pe_systems)


    # Brayton Analysis

    # h3_b = 1595e3 # J/kgK
    # h1_b = 4569e3 # J/kgK
    # h2_b = 3655e3 # J/kgK
    # h4_b = 2099e3 # J/kgK
    h1_b = 4560  # turbine inlet
    h2_b = 3015  # 2355 # turbine exit
    h3_b = 1571  # compressor inlet
    h4_b = 2355  # compressor exit

    mdot_b = Pheat / (h1_b - h4_b)
    Pcomp = mdot_b * (h4_b - h3_b)
    e_comp = 0.75
    Pturb = mdot_b * (h1_b - h2_b)
    e_brayton_actual = 1 - (448 - 300) / (873 - 578)
    e_brayton_actual_2 = (h1_b - h2_b - ((h4_b - h3_b) / e_comp)) / (h1_b - h4_b)
    e_brayton_actual_3 = (h1_b - h2_b) / (h1_b - h4_b)
    e_brayton = 0.4  # ideal brayton eff

    Pheat_available = Pheat  # - Pcomp/e_comp
    Pe_total_b = Pheat_available * e_brayton_actual_2

    Pe_systems_brayton = Paux / e_aux + Pcryo / e_cryo + 2 * Ppump_Flibe / e_pump_Flibe  # + Pcomp/e_comp
    Pe_avail_b = Pe_total_b - Pe_systems_brayton
    Qe_brayton = Pe_total_b / Pe_systems_brayton
    return Qe_rankine,Qe_brayton, Pe_avail,Pe_avail_b

Pfus = np.arange(300e6,500e6,10e6)

Qe_r_array = []
Qe_b_array = []
Pe_r_array = []
Pe_b_array = []
for Pfus_i in Pfus:
    Qe_r,Qe_b,Pe_r,Pe_b = calcQe(Pfus_i)
    Qe_r_array = np.append(Qe_r_array,Qe_r)
    Qe_b_array = np.append(Qe_b_array, Qe_b)
    Pe_r_array = np.append(Pe_r_array, Pe_r)
    Pe_b_array = np.append(Pe_b_array, Pe_b)

# plt.figure()
# plt.title('Pfus Scan')
# plt.plot(Pfus,Qe_r_array,label='Qe Rankine')
# plt.plot(Pfus,Qe_b_array,label='Qe Brayton')
# plt.legend()
#
# plt.figure()
# plt.title('Pfus Scan')
# plt.plot(Pfus,Pe_r_array/1e6,label='Pe Rankine')
# plt.plot(Pfus,Pe_b_array/1e6,label='Pe Brayton')
# plt.legend()



fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel(r'$P_{fus}$ [W]',fontsize=15)
ax1.set_ylabel(r'$Q_e$', color=color,fontsize=15)
ax1.plot(Pfus, Qe_r_array, label=r'$Q_e$ Rankine', color=color)
ax1.plot(Pfus, Qe_b_array, label=r'$Q_e$ Brayton', linestyle='dashdot', color=color)
ax1.tick_params(axis='y', labelcolor=color,labelsize=15)
ax1.tick_params(axis='x',labelsize=15)
ax1.legend(loc='upper left',fontsize=15)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel(r'$P_e$ [MWe]', color=color,fontsize=15)  # we already handled the x-label with ax1
ax2.plot(Pfus, Pe_r_array/1e6, label=r'$P_e$ Rankine', color=color)
ax2.plot(Pfus, Pe_b_array/1e6, label=r'$P_e$ Brayton', linestyle='dashed' ,color=color)
ax2.tick_params(axis='y', labelcolor=color,labelsize=15)
ax2.legend(loc='lower right',fontsize=15)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('rankine_brayton_comp.png',dpi=250)
plt.show()


