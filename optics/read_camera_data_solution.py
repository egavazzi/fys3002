#!/usr/bin/env python3
#
# Determine the location of an artificial aurora blob using
# triangulation with two cameras at two geographically separated
# locations: Silkkimuotka and Nikkaluotka.
#
import numpy as n
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.interpolate as sint

# Juha's coordinate system calculations
# WGS84 lat, lon, height to cartesian and back.
# Also observer centric coordinate systems.
import jcoord

# Data from ALIS during artificial aurora campaign
# this is a matlab file
d=sio.loadmat("fys3002/optics/Data_4_triangulation.mat")

# these are star positions for brightest stars in the image field of view,
# only for the Silkkimuotka station
# angles in radians
star_azimuth=d["CheckStars"][:,0]
star_zenith_angle=d["CheckStars"][:,1]

# load the image for station N (Nikkaluokta)
# both are 128x128 pixel images
N_img=d["dN"]
# load the image for station S (Silkkimuotka)
S_img=d["dS"]

# load the camera calibration station N (Nikkaluokta)
# convert from radians to degrees. mapping of (x,y) -> (az,za)
# (za = zenith angle)
# this has been determined by fitting a lens model to stars in the field
# of view of the camera
azimuth_N=180.0*d["azN"]/n.pi
zenith_angle_N=180.0*d["zeN"]/n.pi

# load the camera calibration station N (Silkkimuotka)
# convert from radians to degrees
# mapping of (x,y) -> (az,za)
azimuth_S=180.0*d["azS"]/n.pi
zenith_angle_S=180.0*d["zeS"]/n.pi
# also in radians, for easier comparison with reference stars
azimuth_S_rad=d["azS"]
zenith_angle_S_rad=d["zeS"]

# meshgrid for pixels
# this is needed for intepolation of x,y image pixels to az,za 
x,y=n.meshgrid(n.arange(128),n.arange(128))

# horizontal projection of az,za avoids a 2*pi phase wrap in azimuth angle near 360 deg
# line of sight vectors
hor_S = n.zeros([128*128,2])
hor_S[:,0]=n.sin(azimuth_S_rad.flatten())*n.sin(zenith_angle_S_rad.flatten())
hor_S[:,1]=n.cos(azimuth_S_rad.flatten())*n.sin(zenith_angle_S_rad.flatten())

# check star positions for Silkkimuotka
# interpolate the pixel values for stars based on their predicted azimuth and zenith angle
star_x=sint.griddata(hor_S,x.flatten(),(n.sin(star_azimuth)*n.sin(star_zenith_angle), n.cos(star_azimuth)*n.sin(star_zenith_angle)))
star_y=sint.griddata(hor_S,y.flatten(),(n.sin(star_azimuth)*n.sin(star_zenith_angle), n.cos(star_azimuth)*n.sin(star_zenith_angle)))

# observer location for Nikkaluokta
# lat, lon, alt in deg N, deg E, and height above sea level in meters.
N_lat = d["latlongN"][0,0]
N_lon = d["latlongN"][0,1]
N_alt = d["altN"][0,0]

# observer location for Silkkimuotka
S_lat = d["latlongS"][0,0]
S_lon = d["latlongS"][0,1]
S_alt = d["altS"][0,0]

print("Silkkimuotka lat %1.2f (deg) N lon %1.2f hgt (deg) E %1.2f meters"%(S_lat,S_lon,S_alt))
print("Nikkaluokta lat %1.2f (deg) N lon %1.2f (deg) E hgt %1.2f meters"%(N_lat,N_lon,N_alt))

# Reference: Björn Gustavsson, personal communications
observation_time="1999-02-16T17:40:45Z"

# cartesian position vecotrs for locations of observers
N_location = jcoord.geodetic2ecef(N_lat, N_lon, N_alt)
S_location = jcoord.geodetic2ecef(S_lat, S_lon, S_alt)

# eye-balled centroid location for artificial aurora blob in Nikkaluokta and Silkkimuotka
N_x = 74.7
N_y = 49.8
S_x = 83.00
S_y = 59.9

# plot the images of the artificial aurora blob seen from two stations
plt.subplot(121)
plt.pcolormesh(x,y,N_img,cmap="plasma")
plt.plot([N_x],[N_y],"x")
plt.colorbar()
plt.title("Nikkaluokta")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.subplot(122)
plt.pcolormesh(x,y,S_img,cmap="plasma")
plt.plot([S_x],[S_y],"x")
plt.plot(star_x,star_y,lw=0,marker="o",fillstyle="none",color="black")
plt.title("Silkkimuotka")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.colorbar()
plt.tight_layout()
plt.show()

# plot the calibration (x,y) -> (az,alt)
plt.subplot(221)
plt.pcolormesh(x,y,azimuth_N,cmap="turbo")
plt.colorbar()
plt.title("Nikkaluokta azimuth (deg)")
plt.xlabel("x pixel")
plt.ylabel("y pixel")

plt.subplot(222)
plt.pcolormesh(x,y,zenith_angle_N,cmap="turbo")
plt.title("Nikkaluokta zenith angle (deg)")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.colorbar()

plt.subplot(223)
plt.pcolormesh(x,y,azimuth_S,cmap="turbo")
plt.colorbar()
plt.title("Silkkimuotka azimuth (deg)")
plt.xlabel("x pixel")
plt.ylabel("y pixel")

plt.subplot(224)
plt.pcolormesh(x,y,zenith_angle_S,cmap="turbo")
plt.title("Silkkimuotka zenith angle (deg)")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.colorbar()
plt.tight_layout()
plt.show()

# get the unit vectors of look directions towards the centroid of the blob from
# nikkaluokta and silkkimuotka, so that a triangulation can be made.

# meshgrid for pixels
xy=n.zeros([128*128,2])

# interpolation from x,y to az,ze
xy[:,0]=x.flatten()
xy[:,1]=y.flatten()

# measured azimuth and elevation from Nikkaluotka
N_az=sint.griddata(xy,azimuth_N.flatten(),(N_x,N_y))
N_za=sint.griddata(xy,zenith_angle_N.flatten(),(N_x,N_y))

# measured azimuth and elevation from Nikkaluotka
S_az=sint.griddata(xy,azimuth_S.flatten(),(S_x,S_y))
S_za=sint.griddata(xy,zenith_angle_S.flatten(),(S_x,S_y))

# Unit vectors in ECEF (ITRS) from Nikkaluokta and Silkkimuotka towards artifical aurora blob center
N_unit_vec=jcoord.azel_ecef(N_lat, N_lon, N_alt, N_az, 90.0-N_za)
S_unit_vec=jcoord.azel_ecef(S_lat, S_lon, S_alt, S_az, 90.0-S_za)

# theory matrix
A = n.zeros([2,2])
m = n.zeros(2)

# TODO: populate A and m to estimate the distance to the center of the blob from Nikkaluokta and Silkkimuotka
# along the look direction (unit vector)
# shouldn't be more than about 6 lines of code
m[0] = n.dot(S_location, N_unit_vec) - n.dot(N_location, N_unit_vec)
m[1] = n.dot(S_location, S_unit_vec) - n.dot(N_location, S_unit_vec)
A[0, 0] = 1
A[1, 1] = -1
A[0, 1] = - n.dot(S_unit_vec, N_unit_vec)
A[1, 0] = n.dot(N_unit_vec, S_unit_vec)

# Estimate distance from Nikkaluokta and Silkkimuotka to the center of the blob
# using the least-squares method
xhat=n.linalg.lstsq(A,m)[0]

# evaluate the heating blob position from N and S
# they should be the same (approximately)
blob_N=N_location + xhat[0]*N_unit_vec
blob_S=S_location + xhat[1]*S_unit_vec

# convert to geodetic for N and S positions
r_blob_N = jcoord.ecef2geodetic(blob_N[0],blob_N[1],blob_N[2])
r_blob_S = jcoord.ecef2geodetic(blob_S[0],blob_S[1],blob_S[2])

print("Nikkaluokta blob position %1.2f N %1.2f E %1.2f km over sea level"%(r_blob_N[0],r_blob_N[1],r_blob_N[2]/1e3))
print("Silkkimuokta blob position %1.2f N %1.2f E %1.2f km over sea level"%(r_blob_S[0],r_blob_S[1],r_blob_S[2]/1e3))

# The result should be approximately this, if you've done the programming right:
# 
# Nikkaluokta blob position 69.27 N 19.18 E 236.45 km over sea level
# Silkkimuokta blob position 69.27 N 19.18 E 236.46 km over sea level
