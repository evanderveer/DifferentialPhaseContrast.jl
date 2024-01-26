const IM2PI = Complex(0, 2*π)

real_image(pixels::Matrix{<:Complex}) = real.(pixels)
imaginary_image(pixels::Matrix{<:Complex}) = imag.(pixels)

function dpc(
    images::Vararg{AbstractMatrix{<:Real}, 4};
    order::Vector = [1, 2, 3, 4]
    )
    check_images(images...)

    image_x = images[order[1]] .- images[order[3]]
    image_y = images[order[2]] .- images[order[4]]

    fft_image_x = image_x |> fft |> fftshift
    fft_image_y = image_y |> fft |> fftshift

    vector_image = make_vector_image(fft_image_y, fft_image_x)

    idpc = vector_image |> 
           integrate_vectors |> 
           ifftshift |> 
           ifft |> 
           real_image |> 
           normalize_image

    ddpc = vector_image |> 
           differentiate_vectors |> 
           ifftshift |> 
           ifft |> 
           real_image |> 
           normalize_image

    (idpc, ddpc)
end

function normalize_image(
    image::AbstractMatrix{<:Real}
    )
    if any(size(image) .< 2)
        throw(ArgumentError("image size must be > 1"))
    end
    image_shifted = image .- minimum(image)
    image_shifted ./ maximum(image_shifted)
end

function check_images(
    images::Vararg{T}
    ) where T<:Matrix

    if any(x -> size(x) != size(images[1]), images)
        throw(ArgumentError("all images must be the same size"))
    end
end

function make_vector_image(
    image_y::Matrix{T},
    image_x::Matrix{T}
    ) where T<:Complex
    check_images(image_y, image_x)

    vectors = Matrix{Vector{T}}(undef, size(image_x)...)

    for index in CartesianIndices(vectors)
        vectors[index] = [image_y[index], image_x[index]]
    end

    vectors
end

function integrate_vectors(
    vectors::Matrix{<:Vector{<:Complex}}
    )
    height, width = size(vectors)
    output = Matrix{ComplexF64}(undef, height, width)

    for index in CartesianIndices(output)
        f1 = height/2 - Tuple(index)[1]
        f2 = Tuple(index)[2] - width/2
        if f1 == f2 == 0
            output[index] = zero(Complex)
        else
            output[index] = scalar_multiply(vectors[index], [f1, f2]) / (IM2PI*(f1^2 + f2^2))
        end
    end
    output
end

function differentiate_vectors(
    vectors::Matrix{<:Vector{<:Complex}}
    )
    height, width = size(vectors)
    output = Matrix{ComplexF64}(undef, height, width)

    for index in CartesianIndices(output)
        f1 = height/2 - Tuple(index)[1]
        f2 = Tuple(index)[2] - width/2
        if f1 == f2 == 0
            output[index] = zero(Complex)
        else
            output[index] = scalar_multiply(vectors[index], [f1, f2]) * IM2PI
        end
    end
    output
end

function scalar_multiply(
    vector1::Vector{<:Union{<:Real, <:Complex}},
    vector2::Vector{<:Union{<:Real, <:Complex}}
    )
    vector1 = Complex.(vector1)
    vector2 = Complex.(vector2)
    sum(vector1 .* conj.(vector2))
end

