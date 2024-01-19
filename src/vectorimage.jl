function make_vector_image(
    image_y::Matrix{T},
    image_x::Matrix{T}
) where T<:Real
    if size(image_x) != size(image_y)
        throw(ArgumentError("images must be the same size"))
    end

    vectors = Matrix{Vector{Complex{T}}}(undef, size(image_x)...)

    for index in CartesianIndices(vectors)
        vectors[index] = Complex.([image_y[index], image_x[index]])
    end

    vectors
end

function fft_cooley(
    vectors::Matrix{<:Vector{<:Complex}}
)
    check_image(vectors)
    height, width = size(vectors)

    intermediate_vectors = fft_cooley_recursion(vectors, height, 1, (0, 0))
    output = similar(vectors)

    num = height / 2

    for index in CartesianIndices(output)
        index1, index2 = Tuple(index)
        output[index] = intermediate_vectors[(index1 - num + height) % height, 
                                             (index2 - num + height) % height]
    end

    output
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