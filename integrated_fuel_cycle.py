# Script to optimize a rankine cycle with regeneration for arbitrary starting parameters
# A. Velberg - MIT PSFC

from pyXSteam.XSteam import XSteam
import numpy as np
import matplotlib.pyplot as plt


'''
Calc Rankine:
    Inputs from loop: T1 (max temp)(C), salt pump powers (W), P_in (input power to rankine cycle) (W)
    Inputs chosen by user: P1, P3 (high and low pressures)(bar), eta_t,eta_p_rankine (rankine component isentropic eff), 
                            P_rf, P_cryo, P_mag (subsystem power requirements)(W)

'''
def calcRankine(T1,P_in,p1,p3,eta_t,eta_p_rankine,P_rf,P_cryo,P_mag,P_pump_salts):
    st = XSteam(XSteam.UNIT_SYSTEM_MKS) # Steam table with units # m/kg/sec/°C/bar/W
    h1 = st.h_pt(p1,T1)
    s1 = st.s_pt(p1,T1)
    h3_ideal = st.h_ps(p3,s1)
    h3 = h1-eta_t*(h1-h3_ideal)
    h4 = st.hL_p(p3) # Just liquid, quality = 0
    s4 = st.sL_p(p3)

    # Must run optimization loop to determine properties at 5,6,7. 5 is calculated using entropy and efficiency of 4 at
    # pressure of 2 = 6. then 7 is determined at pressure of 1 using the entropy of 6.

    # Optimize the regenerative heat, Qregen
    pressures_2_6 = np.linspace(p3,p1,100)
    s1_array = np.ones(100)*s1
    f = np.frompyfunc(st.h_ps,2,1)
    h2_array_ideal = f(pressures_2_6,s1)
    h2_array = h1-eta_t*(h1-h2_array_ideal)
    f2 = np.frompyfunc(st.hL_p,1,1)
    h6_array = f2(pressures_2_6)
    y_array = (h6_array-h4)/(h2_array-h4)# fraction of fluid tapped into regenerator
    QR_array = (h1-h6_array)/(1-y_array) # regenerative heats
    QR_max_idx = np.argmax(QR_array) # This is the optimal Qr wrt h6 variation. Ie this chooses the correct idx to use in other arrays to satisfy dQ/dh6 = 0
    # plt.figure()
    # plt.plot(QR_array)
    # plt.show()
    h6_opt = h6_array[QR_max_idx]
    h2_opt = h2_array[QR_max_idx]
    P2_6_opt = pressures_2_6[QR_max_idx]
    y_opt = y_array[QR_max_idx]
    #print('h2, h6, y optimal', h2_opt,h6_opt, y_opt)
    #print('p26',P2_6_opt,'opt idx', QR_max_idx)
    # Now calculate the remaining enthalpies

    # h5 is at pressure of 2,6, use s4
    h5_ideal = st.h_ps(P2_6_opt,s4)
    h5 = h4 + (h5_ideal-h4)/eta_p_rankine

    # h7 is at pressure of 1, use s6
    s6 = st.sL_p(P2_6_opt)
    h7_ideal = st.h_ps(p1,s6)
    h7 = h6_opt + (h7_ideal-h6_opt)/eta_p_rankine
    #print("T7",st.t_ph(P1,h7))
    # Can now calculate efficiencies, power consumption for individual components, etc.
    eff_rankine = 1 - (1-y_opt)*(h3-h4)/(h1-h7) # Taken from my ME 251 notes :)\
    eff_rankine2 = (h1-h2_opt + (1-y_opt)*(h2_opt-h3))/(h1-h7)
    mdot = P_in/((h1-h7)*1e3) # total mdot
    P_pump1 = mdot*(1-y_opt)*(h5-h4)*1e3
    P_pump2 = mdot*(h7-h6_opt)*1e3
    P_turb = (mdot*(h1-h2_opt) + mdot*(1-y_opt)*(h2_opt-h3))*1e3

    # Mechanical turbine power (efficiency of mech-electrical conversion is assumed 1) - aux - flibe pumps - cryo - rankine pumps
    P_systems = P_rf + P_pump_salts + P_cryo + P_mag + P_pump2/0.75 + P_pump1/0.75
    P_electric = (P_turb-P_systems)
    # Electrical Gain
    Qe = P_turb/P_systems

    # print('Efficiency, Pump_1 power [MW], Pump_2 power [MW], P_turb (both turbines combined) [MW], mdot [kg/s]',eff_rankine,P_pump1/1e6,P_pump2/1e6,P_turb/1e6,mdot,eff_rankine2)
    # print('P elec, Qe, Psystems', P_electric/1e6, Qe, P_systems/1e6)

    return P_electric, Qe, eff_rankine,P_turb,P_pump1,P_pump2,mdot

#P_electric, Qe, eff_rankine,P_turb,P_pump1,P_pump2,mdot = calcRankine(550,676.8285042031814e6,300,0.6,0.9,0.75,20e6,1.64e6,0,.5e6)

'''
Calc Secondary Salt Loop
    Inputs from loop: 
    Inputs chosen by user:

'''

def calcstorage(Pfus,tau_ft,tau_dt,Prf,eff_rankine,P_pump1,P_pump2):

    #tau_ft should be 20 min, tau_dt should be 2 min
    eff_hx=0.95                                                 #Effectiveness of HXs to be checked with Ansys (for the moment assumed to be equal)
    eff_1=eff_hx
    eff_2=eff_hx
    eff_3=eff_hx
    eff_4=eff_hx

    eff_t=eff_rankine
    eff_h2o_pumps = 0.75
    eff_ms_pumps = 0.75

    mdot_ms=321.7551935 #kg/s
    Cp_ms = 1.84e3 # J/kg-K
    deltaT_ms = 20 #K
    P_h=mdot_ms*Cp_ms*deltaT_ms ### aprox 12 MW
    eff_h=0.9  #efficiency of the heating ms process

    Plosses = 0.0 # = energy in the core/tau_losses

    Ppumps_on =  (P_pump1 + P_pump2)/eff_h2o_pumps + 0.5e6/eff_ms_pumps  #this could include cryo and magnets power if different with plasma on and off
    Ppumps_off = (P_pump1 + P_pump2)/eff_h2o_pumps + 0.5e6/eff_ms_pumps

    Psto_off=((Pfus*eff_1)-((Ppumps_on-Ppumps_off)/(eff_2*eff_t))-(Prf/(eff_2*eff_t))+(P_h*((1/(eff_h*eff_t*eff_2))-1))-(Plosses/eff_1))/(1+((tau_dt/tau_ft)*(eff_3/eff_4)))

    Psto_on = Psto_off*(tau_dt/tau_ft)*(eff_3/eff_4)

    Pin_on = Pfus*eff_1 - Psto_on

    return Psto_off,Psto_on,Pin_on, Ppumps_off,Ppumps_on, P_h



'''
Converging the fuel cycle


'''

# User chosen parameters
p1 = 300 # bar. Turbine entry
p3 = 0.6 # bar. Condenser entry
eta_t = 0.99 # turbine mech-elec conversion efficiency
eta_p_rankine = 0.75 #  elec-mech conversion efficiency for pumps in Rankine cycle
P_fus = 451e6 # W. Fusion power
P_rf = 40e6 # W. RF power
P_cryo = 1.64e6 # W. Cryo power
P_mag = 0.5e6 # W. Magnet power
tau_dt = 120 # s . Time of interpulse downtime
tau_ft = 20*60 # s . Pulse length


# Initial Guesses
P_in_0 = 424e6 # Guess initial power into rankine. This will be adjusted to account for power draw to thermal storage system during iteration
T1_0 = 550 # C . Guess for rankine cycle hot leg temperature
P_pump_salts_0 = 0.5e6/.75 # Pumping power required for molten salts (ie from pressure drops). Almost certainly an underestimate.


tolerance = 0.1 # How close we want to get in convergence
P_in_diff = [1] # Just arbitrary starting value
idx = 0
while  abs(P_in_diff[idx]) > tolerance and idx < 10:
    if idx == 0:
        Pin = P_in_0
        T1 = T1_0
        P_pump_salts = P_pump_salts_0 # Not strictly necessary as input to rankine until last iteration, but need some value so code will run

    # First, calculate rankine cyle
    P_electric, Qe, eff_rankine, P_turb, P_pump1, P_pump2, mdot = calcRankine(T1,Pin,p1,p3,eta_t,eta_p_rankine,P_rf,P_cryo,P_mag,P_pump_salts)

    # Next, calculate the salt loop parameters using the rankine information
    P_sto_off, P_sto_on, Pin_new,P_pumps_off,P_pumps_on,P_h = calcstorage(P_fus,tau_ft,tau_dt,P_rf,eff_rankine,P_pump1,P_pump2)
    #print('Pin',Pin/1e6,' Pin new', Pin_new/1e6,' Ppump 1', P_pump1/1e6, 'Ppump 2', P_pump2/1e6)

    # Compare the calculated powers

    P_in_diff_i = abs(Pin - Pin_new)
    P_in_diff = np.append(P_in_diff,P_in_diff_i)

    # Update Pin for next rankine iteration
    Pin = Pin_new
    P_pump_salts = 0.5e6/.75
    idx += 1


'''
Print Results
'''
print('P_electric [MW] =', P_electric/1e6 , '\n Qe = ',Qe, '\n Rankine Pump Power [MW]= ',((P_pump1+P_pump2)/.75 + P_pump_salts)/1e6,
      '\n Rankine Efficiency = ', eff_rankine, 'Ph', P_h/1e6)
plt.figure()
plt.plot(P_in_diff,'-o',color='magenta')
plt.yscale('log')
plt.show()