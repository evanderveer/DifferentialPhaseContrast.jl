# DifferentialPhaseContrast.jl

[![Build Status](https://github.com/evanderveer/CTFFT.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/evanderveer/CTFFT.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/evanderveer/DifferentialPhaseContrast.jl/graph/badge.svg?token=TH59ZWLKTI)](https://codecov.io/gh/evanderveer/DifferentialPhaseContrast.jl)

DifferentialPhaseContrast.jl provides functions for calculating integrated and differentiated differential phase contrast scanning transmission electron microscopy (iDPC-STEM) images from data from a four-segment annular dark field detector. 

### Example

```
> using DifferentialPhaseContrast, DelimitedFiles

# Load the data
> A = readdlm("./A.dat") .|> Float64
> B = readdlm("./B.dat") .|> Float64
> C = readdlm("./C.dat") .|> Float64
> D = readdlm("./D.dat") .|> Float64

# Do the calculation. Data must be supplied as <:Matrix{<:Real}. The optional *order* parameter
# determines in which order the images are used, [1, 2, 3, 4] by default.
> idpc, ddpc = dpc(A, B, C, D, order=[1, 2, 3, 4])

# Save the data
> writedlm("./iDPC.dat", idpc)
> writedlm("./dDPC.dat", ddpc)
```