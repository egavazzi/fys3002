import h5py
import numpy as n
import matplotlib.pyplot as plt


h=h5py.File("ex03-meteor/doppler_data.h5","r")
print(h.keys())
doppler = h["doppler_hertz"][()]
# the three dimensional
# radar wave vector in rad/meter
k_radar_east = h["k_radar_east"][()]
k_radar_north = h["k_radar_north"][()]
k_radar_up = h["k_radar_up"][()]
height = h["height"][()]

# calculate the frequencies based on the Bragg wave vectors
# (propagation factors)
k_abs = n.sqrt(k_radar_east**2.0 + k_radar_north**2.0 + k_radar_up**2.0)
# 4.0*n.pi/lambda = |k_bragg|
# lambda = c/f
# c*|k_bragg|/4/pi = f
print("Effective radar frequencies")
print(k_abs*3e8/4/n.pi/1e6)
# print heights  of detections in km
print(height/1e3)
# print doppler shift measurements (in units of hertz)

n_meas = len(doppler)
A = n.zeros([n_meas,2])
A[:,0]=k_radar_east
A[:,1]=k_radar_north
# hertz to rad/s
m = -2.0*n.pi*doppler
xhat=n.linalg.lstsq(A,m)[0]
print(xhat)
# print(h["t_unix_sec"][()])

h.close()
