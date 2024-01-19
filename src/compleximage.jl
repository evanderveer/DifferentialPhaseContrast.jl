norm(x::Complex) = real(x)^2 + imag(x)^2
real_image(pixels::Matrix{<:Complex}) = real.(pixels)
imaginary_image(pixels::Matrix{<:Complex}) = imag.(pixels)


function fft_cooley_inverse(
    pixels::Matrix{<:Complex}
)

end

function fft_cooley_inverse_recursion(
    pixels::Matrix{<:Complex},
    N::Int,
    delta::Int,
    shifts::Tuple{Int, Int}
)

end

function omega_inverse(
    N::Int,
    pow::Int
)
    num = 2*π*pow/N
    return Complex(cos(num), sin(num))
end