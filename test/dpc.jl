@testset "dpc_functions.jl" begin
    @testset "real_image" begin
        for size in [(1,1), (1,2), (2,1), (2,2), (1,10), (10,1), (10,10)]
            @test DifferentialPhaseContrast.real_image(zeros(Complex, size...)) == zeros(Int, size...)
            @test DifferentialPhaseContrast.real_image(ones(Complex, size...)) == ones(Int, size...)

            Random.seed!(1)
            x = DifferentialPhaseContrast.real_image(Complex.(rand(Float64, size...), rand(Float64, size...)))
            Random.seed!(1)
            y = rand(Float64, size...)
            @test x == y
        end

    end

    @testset "imaginary_image" begin
        for size in [(1,1), (1,2), (2,1), (2,2), (1,10), (10,1), (10,10)]
            @test DifferentialPhaseContrast.imaginary_image(zeros(Complex, size...)) == zeros(Int, size...)
            @test DifferentialPhaseContrast.imaginary_image(ones(Complex, size...)) == zeros(Int, size...)

            Random.seed!(1)
            a = rand(Float64, size...)
            x = DifferentialPhaseContrast.imaginary_image(Complex.(rand(Float64, size...), a))
            Random.seed!(1)
            y = rand(Float64, size...)
            @test x == y
        end

    end
end

@testset "integration" begin
    
    A = readdlm("./images/A.dat") .|> Float64
    B = readdlm("./images/B.dat") .|> Float64
    C = readdlm("./images/C.dat") .|> Float64
    D = readdlm("./images/D.dat") .|> Float64

    idpc, ddpc = DifferentialPhaseContrast.dpc(A,B,C,D)

    idpc_test = readdlm("./images/idpc.dat") .|> Float64
    ddpc_test = readdlm("./images/ddpc.dat") .|> Float64

    @test all(idpc .≈ idpc_test)
    @test all(ddpc .≈ ddpc_test)
end