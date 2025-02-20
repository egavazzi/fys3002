import Pkg; Pkg.activate(@__DIR__)

# import Pkg; Pkg.activate(@__DIR__)
using GLMakie
using GeoMakie
using Downloads
using HDF5
using Dates
using AstroLib
GLMakie.activate!()




## Extract the data
filename = joinpath(@__DIR__, "gps150317g.004.hdf5")
fid = h5open(filename, "r")
data = read(fid)
close(fid)

lat = data["Data"]["Array Layout"]["gdlat"]
lon = data["Data"]["Array Layout"]["glon"]
timestamps = data["Data"]["Array Layout"]["timestamps"]
tec = data["Data"]["Array Layout"]["2D Parameters"]["tec"]
dtec = data["Data"]["Array Layout"]["2D Parameters"]["dtec"]






## Plot 3D Earth
n = 1024 ÷ 1 # 2048
θ = LinRange(0, π, n)
φ = LinRange(0, 2π, 2 * n)
x = [cos(φ) * sin(θ) for θ in θ, φ in φ]
y = [sin(φ) * sin(θ) for θ in θ, φ in φ]
z = [cos(θ) for θ in θ, φ in φ]

fig = Figure(size = (800, 600), backgroundcolor = :grey80)
ax = LScene(fig[1, 1], show_axis = false)
surface!(ax, x, y, z;
    color = GeoMakie.earth(),
    shading = NoShading,
    )
rotate_cam!(ax.scene, (deg2rad(-40), deg2rad(150), 0))  # point on Svalbard

# Add a title
title_text = Observable(string(unix2datetime(timestamps[1])))
Label(fig[0, 1], title_text; tellwidth = false)


# And plot the tec
# Initialize Observables
tec_points = tec[1, :, :]
good_idx = findall(!isnan, tec_points)
good_lon = Observable([lon[i[1]] for i in good_idx])
good_lat = Observable([lat[i[2]] for i in good_idx])
good_tec = Observable(tec_points[good_idx])
# Switch to cartesian coordinates
function toCartesian(lon, lat, alt; cxyz = (0, 0, 0))
    x = cxyz[1] + (alt/6500 + 1) * cosd(lat) * cosd(lon+180)
    y = cxyz[2] + (alt/6500 + 1) * cosd(lat) * sind(lon+180)
    z = cxyz[3] + (alt/6500 + 1) * sind(lat)
    return (x, y, z)
end
positions = Observable(toCartesian.(good_lon[], good_lat[], 100))


# Plot the tec points
sc1 = scatter!(ax, positions; color = good_tec, colormap = :plasma,
               colorrange = (0, 40), depthsorting = true,)
Colorbar(fig[1, 2], sc1; label = "TEC units", tellheight = false)

# Add nightshade
function night_shade(unixtime)
    # This is a simplified model of the nightshade on Earth surface.
    # It is assuming that Earth rotation axis is not tilted.
    zen_lon, _ = zenpos(unix2datetime(unixtime), 0, 0, 0)
    zen_lon = rad2deg(zen_lon)
    midnight_lon = 180 - zen_lon
    nightshade_lon = (midnight_lon - 90):(midnight_lon + 90)

    return midnight_lon, nightshade_lon
end
nightshade_lon = Observable(night_shade(timestamps[1])[2])
θ_shade = LinRange(0, π, 360)
φ_shade = Observable(LinRange(deg2rad(nightshade_lon[][1]), deg2rad(nightshade_lon[][end]), 360))
x_shade = Observable([cos(φ) * sin(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500))
y_shade = Observable([sin(φ) * sin(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500))
z_shade = Observable([cos(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500))
surface!(ax, x_shade, y_shade, z_shade; color = fill((:black, 0.4), (180, 180)),
         shading = NoShading)






# Make a function to update the data points
function update_plot!(i_t, ax, nightshade_lon, tec, good_tec, timestamps, φ_shade, x_shade,
                      y_shade, z_shade)
    # update the tec points
    local tec_points = tec[i_t, :, :]
    local good_idx = findall(!isnan, tec_points)
    good_lon.val = [lon[i[1]] for i in good_idx]
    good_lat.val = [lat[i[2]] for i in good_idx]
    positions.val = toCartesian.(good_lon[], good_lat[], 100)
    good_tec[] = tec_points[good_idx]
    notify(positions)

    # update the nightshade
    _, nightshade_lon[] = night_shade(timestamps[i_t])
    φ_shade[] = LinRange(deg2rad(nightshade_lon[][1]), deg2rad(nightshade_lon[][end]), 360)
    x_shade[] = [cos(φ) * sin(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500)
    y_shade[] = [sin(φ) * sin(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500)
    z_shade[] = [cos(θ) for θ in θ_shade, φ in φ_shade[]] .* (1 + 200/6500)

    # update the title
    title_text[] = string(unix2datetime(timestamps[i_t]))
    return nothing
end
#
display(GLMakie.Screen(), fig)




##
for i_t in 1:length(timestamps)
    update_plot!(i_t, ax, nightshade_lon, tec, good_tec, timestamps, φ_shade, x_shade,
                 y_shade, z_shade)
    sleep(0.05)
end
# update_plot!(4, ax, nightshade_lon, tec, good_tec, timestamps, φ_shade, x_shade,
#                       y_shade, z_shade)
