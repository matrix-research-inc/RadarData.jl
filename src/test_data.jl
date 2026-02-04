#=
	Generate simple radar data for testing

    License ID: SEAL_B
=#
export ScatterCenters, generate_vph

"""
    ScatterCenters{T<:AbstractFloat}

Struct containing 3D locations of point scatterers as static vectors, and a 
vector of real-valued amplitudes associated with each scatterer.
"""
struct ScatterCenters{T<:AbstractFloat}
    locations::Vector{SVector{3,T}}
    amplitudes::Vector{T}
end

Base.length(scs::ScatterCenters) = length(scs.locations)
Base.eltype(scs::ScatterCenters) = eltype(scs.amplitudes)
Base.convert(
    ::Type{ScatterCenters{T2}},
    scs::ScatterCenters{T1}
) where {T1,T2<:AbstractFloat} =
    ScatterCenters(map(x -> T2.(x), scs.locations), T2.(scs.amplitudes))

"""
	generate_vph(scs, spatial_freqs, tx_pos_enu, rx_pos_enu = tx_pos_enu, kernel_sign = -1)

Generate Video Phase History (VPH) data.
If a GPU array support package is included, computations will be performed on a GPU device.
VPH is always returned on CPU.
"""
function generate_vph(
    scs::ScatterCenters{T},
    spatial_freqs::AbstractVector{T},
    tx_pos_enu::AbstractMatrix{T},
    rx_pos_enu::AbstractMatrix{T} = tx_pos_enu;
    kernel_sign::T = T(-1.0)
) where {T<:AbstractFloat}
    # Allocate and initialize output memory depending on backend
    data = gpu_zeros(ComplexF32, length(spatial_freqs), size(tx_pos_enu, 2))
    ref_range_m = move_to_gpu(vec(sqrt.(sum(abs2, tx_pos_enu, dims = 1)) .+
                    sqrt.(sum(abs2, rx_pos_enu, dims = 1))))
    sp_freqs = move_to_gpu(spatial_freqs)
    tx_gpu = move_to_gpu(tx_pos_enu)
    rx_gpu = move_to_gpu(rx_pos_enu)
    kernel! = add_scatterer_data!(gpu_backend())
    for idx in 1:length(scs)
        kernel!(
            data,
            sp_freqs,
            tx_gpu,
            rx_gpu,
            ref_range_m,
            gpu_float_type().(kernel_sign),
            gpu_float_type().(scs.locations[idx]),
            gpu_float_type().(scs.amplitudes[idx]),
            ndrange = size(data)
        )
        KernelAbstractions.synchronize(gpu_backend())
    end
    slow_time_s = range(0, size(tx_pos_enu, 2) - 1, length = size(tx_pos_enu, 2))
    return VPH(
        Array(data),
        Array(spatial_freqs) .* c0 ./ 2 ./ pi,
        Array(ref_range_m),
        slow_time_s,
        Array(tx_pos_enu);
        rx_pos_enu = Array(rx_pos_enu),
        kernel_sign = kernel_sign
    )
end

@kernel function add_scatterer_data!(
    data::AbstractMatrix{ComplexF32},
    sp_freqs::AbstractVector{T},
    tx_pos_enu::AbstractMatrix{T},
    rx_pos_enu::AbstractMatrix{T},
    ref_range_m::AbstractVector{T},
    kernel_sign::T,
    sc_loc::AbstractVector{T},
    sc_amp::T
) where {T<:AbstractFloat}
    ridx, cidx = @index(Global, NTuple)
    @inbounds range_m =
        sqrt(abs2(tx_pos_enu[1, cidx] - sc_loc[1]) +
             abs2(tx_pos_enu[2, cidx] - sc_loc[2]) +
             abs2(tx_pos_enu[3, cidx] - sc_loc[3])) +
        sqrt(abs2(rx_pos_enu[1, cidx] - sc_loc[1]) +
             abs2(rx_pos_enu[2, cidx] - sc_loc[2]) +
             abs2(rx_pos_enu[3, cidx] - sc_loc[3])) - ref_range_m[cidx]
    @inbounds data[ridx, cidx] += sc_amp * cis(kernel_sign * sp_freqs[ridx] * range_m)
end