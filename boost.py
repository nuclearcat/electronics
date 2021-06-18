#!/usr/bin/env python3
import numpy as np
#import matplotlib.pyplot as plt

#Equations
#These equations assume that the coils dissipates all current on each cycle.
#Fmin= 1/(2*Tonmax)
#Ipk= V*Ton/L
#Tonmax=Isat*L/V
#Irms= Ipk/sqrt(3),
#Vrms= Vpk*sqrt(Duty_Cycle)
#P= I*V*sqrt(Duty_Cycle)/sqrt(3)
#in the case of a 50% duty cycle
#P= I*V/sqrt(6)
#E= L*I^2/2

#The energy stored in the magnetic field of an inductor can be calculated as
#W = 1/2*L*I^2

def inductor_max_duty(voltage, inductance, isat):
	t_on_max = isat * inductance / voltage
	f_on_max = 1.0/(2.0*t_on_max)
	print("Max on time "+str(t_on_max*1000*1000) + " uS")
	print("Min freq(0.5 duty cycle) "+str(f_on_max/1000) + " khz")
	return t_on_max



# Khan academy bullshit
#def inductor_current_slope_vfixed(voltage, inductance):
#	current_slope = voltage/inductance
#	print("Current slope "+str(current_slope/(1000*1000))+"A/usec")
#	return current_slope

# time constant(tau) for RL circuit
# Current declines to 0.368 of its initial value per each time constant
def inductor_params(inductance, resistance):
	time_constant = inductance / resistance
	print("time constant " + str(time_constant*1000) + "ms")
	inductance_params = {'inductance': inductance, 'resistance': resistance}
	return inductance_params

def current_inductor_time(inductance_params, mosfet, voltage, time):
	i_max = voltage / (inductance_params['resistance'] + mosfet['rds_on'])
	tau = inductance_params['inductance'] / (inductance_params['resistance'] + mosfet['rds_on'])
	expected_current = i_max*(1-np.exp(-1*time/tau))
	#i_max*np.exp((-inductance_params['resistance']/inductance_params['inductance'])*(i-a))
	return expected_current

# Stored energy at specific "current point"
def inductor_energy(inductance, current):
	energy = 1/2 * inductance * (current**2)
	print("Stored energy "+str(energy)+"J")
	return energy

def boost_dutycycle(vin, vout, efficiency):
	duty_cycle = 1-(vin*efficiency/vout)
	return duty_cycle

def energy_stored(inductance, slope, time):
	current = slope * time
	energy = 1/2 * inductance * (current**2)
	return energy

# https://www.powerelectronicsnews.com/the-dc-dc-boost-converter-power-supply-design-tutorial-section-5-1/


Vin_min = 90
Vin_max = 110
Iin_max = 20
Vout = 200
worst_efficiency = 0.8 # Typical efficiency is 0.85 - 0.95
max_duty_cycle = 1-(Vin_min*worst_efficiency/Vout)
min_duty_cycle = 1-(Vin_max*worst_efficiency/Vout)
available_inductance = 33.0 / (1000*1000) # uH
Fsw = 90*1000 # kHz
print("Max-min duty cycle: "+str(min_duty_cycle)+"-"+str(max_duty_cycle))

# Author note: most likely often "recommended inductance" it is relevant if you want to stay in CCM mode

#However, the suggestion of the 20%~40% current ripple ratio does not take in account the package size of
#inductor. At the small output current condition, following the suggestion may result in large inductor that is
#not applicable in a real circuit. Actually, the suggestion is only the start-point or reference for an inductor
#selection. It is not the only factor, or even not an important factor to determine the inductance in the low
#power application of a boost converter

#Keep in mind also that DCM is usually when the output power is low so there isn’t much energy in the DCM ringing. 
# DCM is, in fact, a case where the inductor ripple current is 200% of the average since it falls all the way to zero. 
#Really low power boosts are sometimes designed to run in DCM on purpose since you can really lower the inductance 
#and that makes for tiny and cheap inductors. 
#But the general trend is to stay in CCM and run at high frequency when a small size is necessary.

# Author note: It might be good idea to have RCD snubber if you run in DCM mode

ripple_1 = 0.2 * Iin_max
ripple_2 = 0.4 * Iin_max
inductance_1 = ((max_duty_cycle / Fsw) * Vin_min) / ripple_1
inductance_2 = ((max_duty_cycle / Fsw) * Vin_min) / ripple_2
print("Recommended inductance (read note!) for CCM: " + str(int(inductance_1 * 1000 * 1000)) + " - " + str(int(inductance_2 * 1000 * 1000)) + " uH")

# TODO: min duty cycle!
I_dcm_ccm_thr = min_duty_cycle * (Vin_max / (2*available_inductance*Fsw))
print("DCM/CCM threshold " + str(I_dcm_ccm_thr) + " A")

inductor_ripple = (Vin_min*max_duty_cycle)/(Fsw * available_inductance)
print("Inductor ripple current: "+str(inductor_ripple)+" A")

p2p_ripple1 = min_duty_cycle * (Vin_max / (available_inductance*Fsw))
p2p_ripple2 = max_duty_cycle * (Vin_min / (available_inductance*Fsw))
print("Peak to peak inductance ripple current: "+str(p2p_ripple2)+" to "+str(p2p_ripple1)+" A")

I_out = (Vin_max * Iin_max) / Vout

I_diode_avg = I_out / (1-max_duty_cycle)
I_diode_peak = I_diode_avg + (inductor_ripple / 2)
print("Diode avg current: "+str(I_diode_avg)+"A peak current: "+str(I_diode_peak)+" A")

exit(0);

# We have inductance 33uH with Isat 30A, input voltage 56V
#inductance = 33.0 / (1000*1000) # uH
#inductance_resistance = 1.55 / 1000.0 # mOhm
ind_params = inductor_params(33.0 / (1000*1000), 1.55 / 1000.0) #uH mOhm
mosfet_params = {'rds_on': (56.2/1000)}
isat = 30
voltage_in = 50
voltage_out = 130
freq = 100*1000 # Khz
efficiency = 0.8


ct = current_inductor_time(ind_params, mosfet_params, 56, 1/20000/2)
print(str(ct) + "A")

# test
#inductance = 7.50/1000
#inductance_resistance = 3
#inductor_params(inductance, inductance_resistance)

#current_slope = inductor_current_slope_vfixed(voltage_in,inductance)
#inductor_max_time = inductor_max_duty(voltage_in,inductance,isat)
#max_duty_cycle = boost_dutycycle(voltage_in, voltage_out, efficiency)
#print(str(max_duty_cycle)+" required duty cycle")
#on_time = 1/freq*max_duty_cycle
#energy_each_cycle = energy_stored(inductance, current_slope, on_time)
#print(str(energy_each_cycle*1000) + "mJ each cycle")


