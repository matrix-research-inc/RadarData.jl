#=
	Shortcuts to common FFT-based methods

	License ID: SEAL_B
=#
export fft_spacing, sfft, sifft, vph_to_rd, rd_to_vph

"""
	fft_spacing(fft_len)

Relative indices for window of length `fft_len` centered on the DC-bin of an equivalently sized FFT.
"""
fft_spacing(fft_len) = (-fld(fft_len, 2)):(cld(fft_len, 2)-1)

"""
	sfft(in_array, [dims])

Symmetric FFT, equivalent to `fftshift(fft(ifftshift(input)))`
"""
sfft(in_array, dims = 1:ndims(in_array)) =
    fftshift(fft(ifftshift(in_array, dims), dims), dims)

"""
	sifft(in_array, [dims])

Symmetric IFFT, equivalent to `fftshift(ifft(ifftshift(input)))`
"""
sifft(in_array, dims = 1:ndims(in_array)) =
    fftshift(ifft(ifftshift(in_array, dims), dims), dims)

"""
	vph_to_rd(in_array, kernel_sign = -1)

Create a range-Doppler image from VPH.data using FFTs / IFFTs according to VPH.kernel_sign. 
"""
vph_to_rd(in_array::Matrix, kernel_sign::Float64 = -1.0) =
    kernel_sign < 0 ? sifft(in_array, 1:2) : sfft(in_array, 1:2)

"""
	rd_to_vph(in_array, kernel_sign = -1)

Convert complex range-Doppler image back to spatial frequency domain (VPH.data domain) 
using FFTs / IFFTs based on VPH.kernel_sign. 
"""
rd_to_vph(in_array::Matrix, kernel_sign::Float64 = -1.0) =
    kernel_sign < 0 ? sfft(in_array, 1:2) : sifft(in_array, 1:2)
