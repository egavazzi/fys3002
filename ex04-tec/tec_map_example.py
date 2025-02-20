#!/usr/bin/env python3
import numpy as n
import matplotlib.pyplot as plt
import h5py
import cartopy.crs as ccrs
from datetime import datetime
from cartopy.feature.nightshade import Nightshade

# Get data from:
# I recommend the binned data, as it is a smaller file
# gps150317g.004.hdf5
# http://millstonehill.haystack.mit.edu/
fname="gps150317g.004.hdf5"

h = h5py.File(fname,"r")
lat=h["Data/Array Layout/gdlat"][()]
lon=h["Data/Array Layout/glon"][()]
timestamps=h["Data/Array Layout/timestamps"][()]
tec=h["Data/Array Layout/2D Parameters/tec"][()]
dtec=h["Data/Array Layout/2D Parameters/dtec"][()]


lon_grid,lat_grid=n.meshgrid(lon,lat)
#print(lat_grid.shape)
#print(lon_grid.shape)

#plt.pcolormesh(lon_grid,lat_grid,tec[:,:,0])
#plt.show()


lat_points=lat_grid.flatten()
lon_points=lon_grid.flatten()


# go through all timesteps
for i in range(len(timestamps)):

    # convert array to points to allow plotting data as scatter plot to allow larger data points 
    # at high latitudes where data is sparse
    tec_points=tec[:,:,i].flatten()
    # only plot thoses that are not nan
    good_idx=n.where(n.isnan(tec_points) != True)[0]
    tec_points2=tec_points[good_idx]
    lat_points2=lat_points[good_idx]
    lon_points2=lon_points[good_idx]

    # print progress
    print(i)
    fig = plt.figure(figsize=[10, 5])
    # ax1 for Northern Hemisphere
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    # ax2 for Southern Hemisphere
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(180, -90))

    data_projection = ccrs.PlateCarree()
    this_frame_date=datetime.utcfromtimestamp(timestamps[i])
    for ax in [ax1, ax2]:
        ax.coastlines(zorder=3)
        ax.stock_img()
        ax.gridlines()
        
        ax.add_feature(Nightshade(this_frame_date))

        # make a scatterplot
        mp=ax.scatter(lon_points2,
                      lat_points2,
                      c=tec_points2,
                      s=n.sqrt(n.abs(lat_points2)+1.0)*0.5, # scale scatter plot point size by latitude
                      vmin=0,vmax=50,
                      transform=data_projection,
                      zorder=2,
                      cmap="turbo")
        cb=plt.colorbar(mp,ax=ax)
        cb.set_label("TEC Units")
    # convert unix seconds to UTC string
    time_str=this_frame_date.strftime('%Y-%m-%d %H:%M:%S')
    ax1.set_title(time_str)
    # ensure small borders in plot
    plt.tight_layout()
    # save figures to compose an animation
    plt.savefig("fig-%06d.png"%(i))
    plt.close()
    plt.clf()

