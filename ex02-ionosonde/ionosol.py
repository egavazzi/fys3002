#!/usr/bin/env python3

import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
# https://github.com/space-physics/iri2016
# pip install iri2016
# documentation:
# http://www.physics.mcgill.ca/mist/memos/MIST_memo_46.pdf
import iri2016.profile as iri

from datetime import datetime, timedelta
from matplotlib.pyplot import figure, show

time_date = datetime(2016,7,7,12,0,0)
alt_km_range = (0,1000,1)
glat=69
glon=19

sim = iri.IRI(time_date, [0,1000,1], glat, glon)

ne=sim["ne"]
alt_km=sim["alt_km"]

plt.semilogx(ne,alt_km)
plt.xlabel("$n_e$ (m$^{-3}$)")
plt.ylabel("Altitude (km)")
plt.show()

freqs=n.linspace(0.5e6,16e6,num=1000)
V = n.zeros([len(alt_km),1000])
for fi in range(len(freqs)):
    omp2=(c.e**2.0 * ne / (c.epsilon_0 * c.electron_mass))
    V[:,fi] = n.sqrt(1 - omp2/((2.0*n.pi*freqs[fi])**2.0))

fof2=9.0*n.sqrt(n.max(ne))/1e6
print(fof2)

plt.pcolormesh(freqs/1e6,alt_km,V)
plt.axvline(fof2,color="green",label="O-mode cutoff")
cb=plt.colorbar()
cb.set_label("$v_g/c$")
plt.xlabel("$n_e$ (m$^{-3}$)")
plt.ylabel("Altitude (km)")
plt.xlabel("Frequency (MHz)")
#plt.show()
# delta height in meters 
dx=n.diff(alt_km)[0]*1e3

reflection_height_idx=n.zeros(len(freqs),dtype=int)
group_delay=n.zeros(len(freqs))
for fi in range(len(freqs)):
    ridx=n.where(n.isnan(V[:,fi]))[0]
    if len(ridx)>0:
        reflection_height_idx[fi]=ridx[0]
        group_delay[fi]=2.0*n.sum(dx/(c.c*V[0:(ridx[0]-2),fi]))
    else:
        reflection_height_idx[fi]=-1
        group_delay[fi]=-1
plt.plot(freqs/1e6,alt_km[reflection_height_idx],color="red",label="Reflection height (O-mode)")
plt.plot(freqs/1e6,group_delay*c.c/2.0/1e3,color="blue",label="Virtual height (O-mode)")
plt.ylim([0,1000])
plt.legend()

plt.show()





