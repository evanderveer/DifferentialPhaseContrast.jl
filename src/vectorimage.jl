function make_vector_image(
    image_x::Matrix{<:Real},
    image_y::Matrix{<:Real}
)

function fft_cooley(
    vectors::Matrix{<:Vector{<:Complex}}
)

end

function fft_cooley_recursion(
    vectors::Matrix{<:Vector{<:Complex}},
    N::Int,
    delta::Int,
    shifts::Tuple{Int, Int}
)

end

function omega_forward(
    N::Int,
    pow::Int
)
    num = -2*π*pow/N
    return Complex(cos(num), sin(num))
end