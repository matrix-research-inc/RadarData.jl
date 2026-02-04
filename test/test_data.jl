#=
    License ID: SEAL_B
=#

using RadarData
using StaticArrays
using KernelAbstractions

center_freq = 8.0e9
bandwidth = 0.2e9
num_freqs = 601
freq_list_hz = RadarData.get_freq_list(center_freq, bandwidth, num_freqs)

# Scattering Centers
scene_radius_x = 35.0
scene_radius_y = 35.0
dx_scene = 8.0
dy_scene = 8.0
scat_centers = [
    SVector{3}(idx1, idx2, idx3) for
    idx1 in (-scene_radius_x):dx_scene:(scene_radius_x-5),
    idx2 in (-scene_radius_y):dy_scene:(scene_radius_y-5), idx3 in 0.0
]
scat_amps = ones(length(scat_centers))
scs = RadarData.ScatterCenters(vec(scat_centers), scat_amps)

# Sensor path
half_angle = 4 * pi / 180
center_angle = 0 * pi / 180
standoff = 150
num_pulses = 1001
angle_list = range(-half_angle, half_angle, length = num_pulses) .+ center_angle
el_angle = 45 * pi / 180
tx_pos_enu = vcat(
    standoff .* cos.(angle_list)' .* cos(el_angle),
    standoff .* sin.(angle_list)' .* cos(el_angle),
    fill(standoff * sin(el_angle), length(angle_list))'
)
tx_pos_enu[1, :] += range(0, 20, length = num_pulses)
rx_pos_enu = tx_pos_enu
spatial_freqs = 2 .* pi ./ RadarData.c0 .* freq_list_hz

# Generate data (CPU)
@time vph = RadarData.generate_vph(scs, spatial_freqs, tx_pos_enu)

# Float32 generation
scs = convert(ScatterCenters{Float32}, scs)
spatial_freqs = Vector{Float32}(spatial_freqs)
tx_pos_enu = Float32.(tx_pos_enu)
rx_pos_enu = tx_pos_enu

@time vph = RadarData.generate_vph(scs, spatial_freqs, tx_pos_enu)

# Generate on Apple Silicon, limited to Float32
using Metal
spatial_freqs_mtl = MtlArray(spatial_freqs)
tx_pos_enu_mtl = MtlArray(tx_pos_enu)

# Timing includes shuttling data back from GPU
@time vph = RadarData.generate_vph(scs, spatial_freqs_mtl, tx_pos_enu_mtl)

# Generate on CUDA 
using CUDA
spatial_freqs_cu = CuArray(spatial_freqs)
tx_pos_enu_cu = CuArray(tx_pos_enu)

@time vph = RadarData.generate_vph(scs, spatial_freqs_cu, tx_pos_enu_cu)

using InteractiveViz, WGLMakie
using FFTW
heatmap(log10.(abs2.(fftshift(fft(vph.data)))))