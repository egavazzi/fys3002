import h5py
import matplotlib.pyplot as plt
import numpy as n
import scipy.interpolate as sinterp
import scipy.integrate as sinteg

from pymsis import msis


# Precomputed ionization profile using Sergienko and Ivanov (1993)
# 
h=h5py.File("q_10keV.h5","r")
# ionization rate as a function of altitude of a
# precipitation flux of one electron per squaremeter per second
# primary auroral electron with energy of 10 keV
q=h["q"][()]
# height in meters
hgt=h["hgt"][()]*1e3
# energy of primary precipitating electron
E=h["energy"][()]
h.close()



# PyMSIS
# See documentation: https://pypi.org/project/pymsis/
msis_dates = n.array([n.datetime64("2024-01-01T00:00")])

lat = 69.65
lon = 18.96

data = msis.run(msis_dates, lon, lat , hgt/1e3, geomagnetic_activity=-1)
# N_2
n_N2=data[0,0,0,:,1]
# O_2
n_O2=data[0,0,0,:,2]
# O
n_O=data[0,0,0,:,3]
plt.semilogx(n_N2, hgt/1e3,label="$N_2$")
plt.semilogx(n_O2, hgt/1e3,label="$O_2$")
plt.semilogx(n_O, hgt/1e3,label="$O$")
plt.legend()
plt.xlabel("Neutral density (m$^{-3}$)")
plt.ylabel("Height (km)")
plt.tight_layout()
plt.show()
print(data.shape)

# E-region recombination coefficient
# Geldhill, The effective recombination coefficient of electrons
# in the ionosphere between 50 and 150 km (1986)
recombination_coeff=3e-13

# plot ionization rate
plt.semilogx(q,hgt/1e3)
plt.title("Ionization-rate $q(h)$\n(for an incident flux of 1 m$^{-2}$s$^{-1}$ precipitating electrons E=%1.1f keV)"%(E/1e3))
plt.xlabel("Ionization rate (m$^{-3}$ s$^{-1}$)")
plt.ylabel("Height (km)")
plt.tight_layout()
plt.show()

# altitude step in meters
dh=n.diff(hgt)[0]

# total number of ionizations per second per incident primary electron
# (should be about 300 for a 10 keV electron, source: Björn)
Q_tot=sinteg.trapezoid(q,hgt,dx=dh)

# steady state electron density with a flux of
# 1e13 precipitating electrons per squaremeter per second
# 
n_e = n.sqrt(1e13*q/recombination_coeff)

plt.semilogx(n_e,hgt/1e3)
plt.title("Steady-state electron density $n_e(h) = \\sqrt{q(h)/\\alpha}$")
plt.xlabel("Electron density ($m^{-3}$)")
plt.ylabel("Height (km)")
plt.show()

# collision frequency  (Brekke's book)
# assume constant temperature.
# Brekke's n_n is in 1/cm^3
Te=300.0
# use effective collision frequency when nu_en << \omega
# Hargreaves (1969)
nu_en = (5/2)*5.4e-10*(n_N2/1e6)*n.sqrt(Te)

plt.semilogx(nu_en,hgt/1e3)
plt.title("Electron-neutral collision frequency")
plt.xlabel("$\\frac{5}{2}\\nu_{en}$ (s$^{-1}$)")
plt.ylabel("Height (km)")
plt.tight_layout()
plt.show()



