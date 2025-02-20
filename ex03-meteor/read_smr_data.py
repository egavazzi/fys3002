import h5py
import numpy as n
import matplotlib.pyplot as plt

# open file
h=h5py.File("doppler_data.h5","r")
# print data fields
print(h.keys())
# doppler shift estimates
doppler = h["doppler_hertz"][()]
# the three dimensional
# radar wave vector in rad/meter
k_radar_east = h["k_radar_east"][()]
k_radar_north = h["k_radar_north"][()]
k_radar_up = h["k_radar_up"][()]
height = h["height"][()]
# measurement times
t_unix = h["t_unix_sec"][()]

# this is how many seconds of time the measurements span
print(n.max(t_unix)-n.min(t_unix))

# this is how many doppler measurements there are
print(len(doppler))

# calculate the frequencies based on the Bragg wave vectors
# (propagation factors)
k_abs = n.sqrt(k_radar_east**2.0 + k_radar_north**2.0 + k_radar_up**2.0)
# 4.0*n.pi/lambda = |k_bragg|
# lambda = c/f
# c*|k_bragg|/4/pi = f
#print("Effective radar frequencies")
#print(k_abs*3e8/4/n.pi/1e6)
# print heights  of detections in km
#print(height/1e3)
# print doppler shift measurements (in units of hertz)

n_meas = len(doppler)
A = n.zeros([n_meas,2])
# TODO: estimate zonal and meridional wind
h.close()
