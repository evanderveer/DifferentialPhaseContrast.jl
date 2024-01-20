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

    intermediate_vectors = fft_cooley_recursion(vectors, height, 1, (1, 1))
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
    if N == 1; return vectors[shifts...]; end;
    
    N1 = Int(N / 2)

    fft_intermediates = []
    for extra_shift in [(0, 0), (0, delta), (delta, 0), (delta, delta)]
        push!(fft_intermediates, fft_cooley_recursion(vectors, N1, 2*delta, shifts .+ extra_shift))
    end


    output = similar(vectors)

    for pow1 in 1:N1
    for pow2 in 1:N1
        current_intermediates = getindex.(fft_intermediates, pow1, pow2)
        current_intermediates[2] *= omega_forward(N, pow2)
        current_intermediates[3] *= omega_forward(N, pow1)
        current_intermediates[4] *= omega_forward(N, pow1 + pow2)

        output[pow1, pow2] = sum(current_intermediates)
        output[pow1, pow2 + N1] = sum(current_intermediates .* [1,-1,1,-1])
        output[pow1 + N1, pow2] = sum(current_intermediates .* [1,1,-1,-1])
        output[pow1 + N1, pow2 + N1] = sum(current_intermediates .* [1,-1,-1,1])
    end
    end
    
    output
end

function omega_forward(
    N::Int,
    pow::Int
)
    num = -2*π*pow/N
    return Complex(cos(num), sin(num))
end