@testset "idpc_functions.jl" begin
    @testset "real_image" begin
        for size in [(1,1), (1,2), (2,1), (2,2), (1,10), (10,1), (10,10)]
            @test iDPC.real_image(zeros(Complex, size...)) == zeros(Int, size...)
            @test iDPC.real_image(ones(Complex, size...)) == ones(Int, size...)

            Random.seed!(1)
            x = iDPC.real_image(Complex.(rand(Float64, size...), rand(Float64, size...)))
            Random.seed!(1)
            y = rand(Float64, size...)
            @test x == y
        end

    end

    @testset "imaginary_image" begin
        for size in [(1,1), (1,2), (2,1), (2,2), (1,10), (10,1), (10,10)]
            @test iDPC.imaginary_image(zeros(Complex, size...)) == zeros(Int, size...)
            @test iDPC.imaginary_image(ones(Complex, size...)) == zeros(Int, size...)

            Random.seed!(1)
            a = rand(Float64, size...)
            x = iDPC.imaginary_image(Complex.(rand(Float64, size...), a))
            Random.seed!(1)
            y = rand(Float64, size...)
            @test x == y
        end

    end
end

@testset "integration" begin
    
    A = load("./images/seg0.png") .|> Float64
    B = load("./images/seg1.png") .|> Float64
    C = load("./images/seg2.png") .|> Float64
    D = load("./images/seg3.png") .|> Float64

    idpc = iDPC.idpc(A,B,C,D)

    idpc_test = readdlm("./images/idpc.dat") .|> Float64

    @test all(idpc .≈ idpc_test)
end