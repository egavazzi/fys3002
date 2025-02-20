#!/usr/bin/env python3

import numpy as n
import matplotlib.pyplot as plt
# https://github.com/space-physics/iri2016
# pip install iri2016
# documentation:
# http://www.physics.mcgill.ca/mist/memos/MIST_memo_46.pdf
import iri2016.profile as iri

from datetime import datetime, timedelta
from matplotlib.pyplot import figure, show

time_date = datetime(2016,7,7,12,0,0)
alt_km_range = (50,1000,1)
glat=69
glon=19

sim = iri.IRI(time_date, [50,1000,1], glat, glon)

ne=sim["ne"]
alt_km=sim["alt_km"]

plt.semilogx(ne,alt_km)
plt.xlabel("$n_e$ (m$^{-3}$)")
plt.ylabel("Altitude (km)")
plt.show()
