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
struct ScatterCenters{T <: AbstractFloat}
    locations::Vector{SVector{3, T}}
    amplitudes::Vector{T}
end

Base.length(scs::ScatterCenters) = length(scs.locations)
Base.eltype(scs::ScatterCenters) = eltype(scs.amplitudes)
function Base.convert(
        ::Type{ScatterCenters{T2}},
        scs::ScatterCenters{T1}
) where {T1, T2 <: AbstractFloat}
    ScatterCenters(map(x -> T2.(x), scs.locations), T2.(scs.amplitudes))
end

"""
	generate_vph(scs, spatial_freqs, tx_pos_enu;
                 rx_pos_enu = tx_pos_enu,
                 slow_time_s = 0:(size(tx_pos_enu, 2) - 1,
                 kernel_sign = -1,
                 target_motion = (identity))

Generate Video Phase History (VPH) data.
If a GPU array support package is included, computations will be performed on a GPU device.
VPH is always returned on CPU.
Optionally provide `target_motion`, which is a function of the form `(location, time) -> new_location`.
This function must not allocate on the heap (tip: use StaticArrays).
The default is no motion.
"""
function generate_vph(
        scs::ScatterCenters{Float64},
        spatial_freqs::AbstractVector{Float64},
        tx_pos_enu::AbstractMatrix{Float64};
        rx_pos_enu::AbstractMatrix{Float64} = tx_pos_enu,
        slow_time_s::AbstractVector{Float64} =
        range(0, size(tx_pos_enu, 2) - 1, size(tx_pos_enu, 2)),
        kernel_sign::Float64 = -1.0,
        target_motion::Function = (x, y) -> x
)
    # Allocate and initialize output memory depending on backend
    data = gpu_zeros(ComplexF32, length(spatial_freqs), size(tx_pos_enu, 2))
    ref_range_m = move_to_gpu(vec(sqrt.(sum(abs2, tx_pos_enu, dims = 1)) .+
                                  sqrt.(sum(abs2, rx_pos_enu, dims = 1))))
    sp_freqs = move_to_gpu(spatial_freqs)
    tx_gpu = move_to_gpu(tx_pos_enu)
    rx_gpu = move_to_gpu(rx_pos_enu)
    time_gpu = move_to_gpu(slow_time_s)
    kernel! = add_scatterer_data!(gpu_backend())
    for idx in eachindex(scs.locations)
        kernel!(
            data,
            sp_freqs,
            tx_gpu,
            rx_gpu,
            time_gpu,
            ref_range_m,
            gpu_float_type().(kernel_sign),
            gpu_float_type().(scs.locations[idx]),
            gpu_float_type().(scs.amplitudes[idx]),
            target_motion,
            ndrange = size(data)
        )
        KernelAbstractions.synchronize(gpu_backend())
    end
    return VPH(
        Array(data),
        Array(spatial_freqs) .* c0 ./ 2 ./ pi,
        Array(ref_range_m),
        Array(slow_time_s),
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
        slow_time_s::AbstractVector{T},
        ref_range_m::AbstractVector{T},
        kernel_sign::T,
        sc_loc::AbstractVector{T},
        sc_amp::T,
        target_motion::Function
) where {T <: AbstractFloat}
    ridx, cidx = @index(Global, NTuple)
    # Update scatterer location per pulse, `target_motion` must not allocate
    @inbounds up_loc = target_motion(sc_loc, slow_time_s[cidx])
    @inbounds range_m = sqrt(abs2(tx_pos_enu[1, cidx] - up_loc[1]) +
                             abs2(tx_pos_enu[2, cidx] - up_loc[2]) +
                             abs2(tx_pos_enu[3, cidx] - up_loc[3])) +
                        sqrt(abs2(rx_pos_enu[1, cidx] - up_loc[1]) +
                             abs2(rx_pos_enu[2, cidx] - up_loc[2]) +
                             abs2(rx_pos_enu[3, cidx] - up_loc[3])) - ref_range_m[cidx]
    @inbounds data[ridx, cidx] += sc_amp * cis(kernel_sign * sp_freqs[ridx] * range_m)
end
