#=
    License ID: SEAL_B
=#
using MAT
using RadarData
tsc_ext = Base.get_extension(RadarData, :TSC)

file_name = joinpath(ENV["DATA_DIR"], "TSC/mat/data_1_01.mat")

gps_week = 2260
@time tsc_vph = tsc_ext.load_tsc(file_name, gps_week, sub_inds = 1:10000)

# Make range-Doppler map
using SARImaging
using SARPlots
rd = vph_to_rd(tsc_vph)
save_sar("tsc_rd.png", rd)

# Apply autofocus
mapdrift!(
    tsc_vph,
    500,
    500;
    max_iters = 10,
    drift_func = fit_drift_spline(10),
    verbose = true
)
rd = vph_to_rd(tsc_vph)
save_sar("tsc_rd_mapdrift.png", rd)

# Try polar format algorithm
tsc_pfa = polar_format(tsc_vph; inner_box = false, distortion = false)
save_sar("tsc_pfa.png", tsc_pfa)

# Looks weird, try backprojecting
image_grid = default_grid(tsc_vph)
image_grid = downsample(image_grid, [4, 2])
tsc_bp = zeros(eltype(tsc_vph), size(image_grid))
@time backproject!(tsc_bp, image_grid, tsc_vph)

save_image("tsc_bp.png", hist_adaptive_scale(tsc_bp))
