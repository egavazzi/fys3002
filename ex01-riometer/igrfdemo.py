import numpy as n
import scipy.constants as c
import ppigrf
import matplotlib.pyplot as plt
from datetime import datetime

# Look at documentation here:
# https://github.com/klaundal/ppigrf/

lat = 69.65
lon = 18.96
# 0 to 500 km
h = n.linspace(0,500,500)
date = datetime(2024,1,17)

Be, Bn, Bu = ppigrf.igrf(lon,lat,h,date)
print(Be.shape)
plt.subplot(131)
plt.plot(n.sqrt(Be[0,:]**2.0 + Bn[0,:]**2.0 + Bu[0,:]**2.0),h)
plt.xlabel("$|\\vec{B}|$ (nT)")
plt.ylabel("Height (km)")
plt.subplot(132)
# Absolute value B in nanotesla
B_abs=n.sqrt(Be[0,:]**2.0 + Bn[0,:]**2.0 + Bu[0,:]**2.0)/1e9
plt.plot(c.e*B_abs/c.electron_mass,h)
plt.xlabel("$\\omega_c$ (rad/s)")
plt.ylabel("Height (km)")

plt.subplot(133)
# Absolute value B in nanotesla
plt.plot(c.e*B_abs/c.electron_mass/2/n.pi,h)
plt.xlabel("$f_c$ (Hz)")
plt.ylabel("Height (km)")
plt.tight_layout()
plt.savefig("gyrofreq.png")
plt.tight_layout()
plt.show()
