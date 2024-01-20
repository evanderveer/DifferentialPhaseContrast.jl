norm(x::Complex) = real(x)^2 + imag(x)^2
real_image(pixels::Matrix{<:Complex}) = real.(pixels)
imaginary_image(pixels::Matrix{<:Complex}) = imag.(pixels)


function fft_cooley_inverse(
    pixels::Matrix{<:Complex}
)
    check_image(pixels)
    height, width = size(pixels)
    num = height / 2

    pixels1 = similar(pixels)

    for index in CartesianIndices(pixels1)
        pixels1[index] = fftImage[(Tuple(index)[1] - num + height) % height + 1,
                                  (Tuple(index)[2] - num + width) % width + 1]
    end
    
    pixels2 = fftCooleyRecursion(height, 1, 0, 0, pixels1, true)
    pixels2 ./ (height^2)
end

function fft_cooley_inverse_recursion(
    pixels::Matrix{<:Complex},
    N::Int,
    delta::Int,
    shifts::Tuple{Int, Int}
)
    if N == 1; return pixels[shifts...]; end;
    
    N1 = N / 2

    fft_intermediates = []
    for extra_shift in [(0, 0), (0, delta), (delta, 0), (delta, delta)]
        push!(fft_intermediates, fft_cooley_inverse_recursion(pixels, N1, 2*delta, shifts .+ extra_shift))
    end

    output = similar(pixels)

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

function omega_inverse(
    N::Int,
    pow::Int
)
    num = 2*π*pow/N
    return Complex(cos(num), sin(num))
end

function check_image(
    matrix::Matrix
)
    height, width = size(matrix)
    if width != height
        throw(ArgumentError("image width and height are not the same"))
    end
    if !ispow2(width)
        throw(ArgumentError("image size must be a power of 2")) 
    end
end