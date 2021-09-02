#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

skip = 0
inp_file = np.loadtxt("energy.txt", skiprows=1)
time = inp_file[skip:,0]
total_energy = inp_file[skip:,1]
potential = inp_file[skip:,2]
kinetic = inp_file[skip:,3]

plt.plot(time, total_energy, label="Total Energy")
#plt.plot(time, potential*96.45, label="Potential")
#plt.plot(time, kinetic*96.45, label="Kinetic")
plt.title("Water")
plt.xlabel("Time [ps]")
plt.ylabel("Energy [kJ/mol]")
plt.legend()
plt.tight_layout()
plt.savefig("eng_pol_acnt.png")
