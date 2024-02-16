# Script to optimize a rankine cycle with regeneration for arbitrary starting parameters
# A. Velberg - MIT PSFC

from pyXSteam.XSteam import XSteam
import numpy as np
import matplotlib.pyplot as plt


# Temperatures in K, Pressures in bar, Pin in Watts
def calcRankine(T1,P1,P3,eff_turb,eff_pump,Pin):
    st = XSteam(XSteam.UNIT_SYSTEM_MKS) # Steam table with units # m/kg/sec/°C/bar/W
    h1 = st.h_pt(P1,T1)
    s1 = st.s_pt(P1,T1)
    h3_ideal = st.h_ps(P3,s1)
    h3 = h1-eff_turb*(h1-h3_ideal)
    h4 = st.hL_p(P3) # Just liquid, quality = 0
    s4 = st.sL_p(P3)

    # Must run optimization loop to determine properties at 5,6,7. 5 is calculated using entropy and efficiency of 4 at
    # pressure of 2 = 6. then 7 is determined at pressure of 1 using the entropy of 6.

    # Optimize the regenerative heat, Qregen
    pressures_2_6 = np.linspace(P3,P1,100)
    s1_array = np.ones(100)*s1
    f = np.frompyfunc(st.h_ps,2,1)
    h2_array_ideal = f(pressures_2_6,s1)
    h2_array = h1-eff_turb*(h1-h2_array_ideal)
    print('here1')
    f2 = np.frompyfunc(st.hL_p,1,1)
    h6_array = f2(pressures_2_6)
    print('here2')
    y_array = (h6_array-h4)/(h2_array-h4)# fraction of fluid tapped into regenerator
    QR_array = (h1-h6_array)/(1-y_array) # regenerative heats
    QR_max_idx = np.argmax(QR_array) # This is the optimal Qr wrt h6 variation. Ie this chooses the correct idx to use in other arrays to satisfy dQ/dh6 = 0
    plt.figure()
    plt.plot(QR_array)
    plt.show()
    h6_opt = h6_array[QR_max_idx]
    h2_opt = h2_array[QR_max_idx]
    P2_6_opt = pressures_2_6[QR_max_idx]
    y_opt = y_array[QR_max_idx]
    print('h2, h6, y optimal', h2_opt,h6_opt, y_opt)
    print('p26',P2_6_opt,'opt idx', QR_max_idx)
    # Now calculate the remaining enthalpies

    # h5 is at pressure of 2,6, use s4
    h5_ideal = st.h_ps(P2_6_opt,s4)
    h5 = h4 + (h5_ideal-h4)/eff_pump

    # h7 is at pressure of 1, use s6
    s6 = st.sL_p(P2_6_opt)
    h7_ideal = st.h_ps(P1,s6)
    h7 = h6_opt + (h7_ideal-h6_opt)/eff_pump
    print("T7",st.t_ph(P1,h7))
    # Can now calculate efficiencies, power consumption for individual components, etc.
    eff_rankine = 1 - (1-y_opt)*(h3-h4)/(h1-h7) # Taken from my ME 251 notes :)\
    eff_rankine2 = (h1-h2_opt + (1-y_opt)*(h2_opt-h3))/(h1-h7)
    mdot = Pin/((h1-h7)*1e3) # total mdot
    P_pump1 = mdot*(1-y_opt)*(h5-h4)/1e3
    P_pump2 = mdot*(h7-h6_opt)/1e3
    P_turb = (mdot*(h1-h2_opt) + mdot*(1-y_opt)*(h2_opt-h3))/1e3

    # Mechanical turbine power (efficiency of mech-electrical conversion is assumed 1) - aux - flibe pumps - cryo - rankine pumps
    P_systems = 57e3 + 0.883e3 + 1.64e3 + P_pump2 + P_pump1
    P_electric = P_turb-P_systems/1e3
    # Electrical Gain
    Qe = P_turb/P_systems*1e3

    print('Efficiency, Pump_1 power [MW], Pump_2 power [MW], P_turb (both turbines combined) [MW], mdot [kg/s]',eff_rankine,P_pump1,P_pump2,P_turb,mdot,eff_rankine2)
    print('P elec, Qe, Psystems', P_electric, Qe, P_systems/1e3)

    return eff_rankine, P_pump1, P_pump2,P_turb,mdot,P_electric, Qe


#eff_rankine, P_pump1, P_pump2,P_turb,mdot,P_electric, Qe = calcRankine(550,300,0.6,0.9,0.75,389.524e6)
eff_rankine, P_pump1, P_pump2,P_turb,mdot,P_electric, Qe = calcRankine(550,150,0.6,0.9,0.75,426.71e6)

# use p1=150 for subcrit, p1 = 300 for supercrit

###eff_rankine, P_pump1, P_pump2,P_turb,mdot,P_electric, Qe = calcRankine(550,300,0.6,0.9,0.75,395e6*1.2*.95**2)
