function normalize_image(
    image::Matrix{<:Real}
)
    image_shifted = image .- minimum(image)
    image_shifted ./ maximum(image_shifted)
end