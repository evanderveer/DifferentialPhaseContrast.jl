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

    @testset "normalize_image" begin
        # Test case 1: Normalizing a positive image
        img1 = [1.0 2.0; 3.0 4.0]
        @test isapprox(DifferentialPhaseContrast.normalize_image(img1), [0.0 1/3; 2/3 1.0]) 

        # Test case 2: Normalizing a negative image
        img2 = [-1.0 -2.0; -3.0 -4.0]
        @test isapprox(DifferentialPhaseContrast.normalize_image(img2), [1.0 2/3; 1/3 0.0]) 

        # Test case 3: Normalizing an image with zero values
        img3 = [0.0 1.0; 2.0 3.0]
        @test isapprox(DifferentialPhaseContrast.normalize_image(img3), [0.0 1/3; 2/3 1.0])

        # Test case 4: Normalizing an empty image
        img4 = zeros(0, 0)
        @test_throws ArgumentError DifferentialPhaseContrast.normalize_image(img4)
    end

    @testset "check_images" begin
        # Test case 1: Images with the same size
        img1 = rand(3, 3)
        img2 = rand(3, 3)
        img3 = rand(3, 3)
        DifferentialPhaseContrast.check_images(img1, img2, img3)  # Should not throw an error
        DifferentialPhaseContrast.check_images(img1)  # Should not throw an error
        DifferentialPhaseContrast.check_images(img1, img2)  # Should not throw an error

        # Test case 2: Images with different sizes
        img4 = rand(4, 4)
        img5 = rand(3, 3)
        @test_throws ArgumentError DifferentialPhaseContrast.check_images(img1, img4)
        @test_throws ArgumentError DifferentialPhaseContrast.check_images(img1, img2, img4)
    end

    @testset "make_vector_image" begin
        # Test case 1: Valid complex images
        img_y1 = [1.0 + 2.0im 3.0 + 4.0im; 5.0 + 6.0im 7.0 + 8.0im]
        img_x1 = [9.0 + 10.0im 11.0 + 12.0im; 13.0 + 14.0im 15.0 + 16.0im]
        result1 = DifferentialPhaseContrast.make_vector_image(img_y1, img_x1)
        @test size(result1) == (2, 2)
        @test typeof(result1) == Matrix{Vector{ComplexF64}}
        @test result1[CartesianIndex((1, 1))] == [1.0 + 2.0im, 9.0 + 10.0im]

        # Test case 2: Invalid real images
        img_y2 = [1.0 2.0; 3.0 4.0]
        img_x2 = [5.0 6.0; 7.0 8.0]
        @test_throws MethodError DifferentialPhaseContrast.make_vector_image(img_y2, img_x2)
       
        # Test case 3: Images with different sizes
        img_y3 = [1.0 + 2.0im 3.0 + 4.0im; 5.0 + 6.0im 7.0 + 8.0im]
        img_x3 = [9.0 + 10.0im 11.0 + 12.0im 1.0 + 2.0im; 13.0 + 14.0im 15.0 + 16.0im 3.0 + 1.0im]
        @test_throws ArgumentError DifferentialPhaseContrast.make_vector_image(img_y3, img_x3)

        # Test case 4: Empty images
        img_y4 = []
        img_x4 = []
        @test_throws MethodError DifferentialPhaseContrast.make_vector_image(img_y4, img_x4)
    end

    @testset "integrate_vectors" begin
        A = zeros(ComplexF64, 2, 2)
        B = zeros(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test DifferentialPhaseContrast.integrate_vectors(C) == zeros(ComplexF64, 2, 2)

        A = ones(ComplexF64, 2, 2)
        B = ones(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test isapprox(DifferentialPhaseContrast.integrate_vectors(C), [0.0+0.0im 0.0-0.159155im; 0.0+0.159155im  0.0-0.0im]) atol=0.00001

        A = ones(ComplexF64, 1, 1)
        B = ones(ComplexF64, 1, 1)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test DifferentialPhaseContrast.integrate_vectors(C) == zeros(ComplexF64, 1, 1)

        A = ones(ComplexF64, 2, 2)
        B = zeros(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test isapprox(DifferentialPhaseContrast.integrate_vectors(C), [0.0+0.0im 0.0-0.0im; 0.0+0.159155im 0.0+0.0795775im]) atol=0.00001
    end

    @testset "differentiate_vectors" begin
        A = zeros(ComplexF64, 2, 2)
        B = zeros(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test DifferentialPhaseContrast.differentiate_vectors(C) == zeros(ComplexF64, 2, 2)

        A = ones(ComplexF64, 2, 2)
        B = ones(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test isapprox(DifferentialPhaseContrast.differentiate_vectors(C), [0.0+0.0im 0.0+6.28319im; 0.0-6.28319im  0.0-0.0im]) atol=0.00001

        A = ones(ComplexF64, 1, 1)
        B = ones(ComplexF64, 1, 1)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test DifferentialPhaseContrast.differentiate_vectors(C) == zeros(ComplexF64, 1, 1)

        A = ones(ComplexF64, 2, 2)
        B = zeros(ComplexF64, 2, 2)
        C = DifferentialPhaseContrast.make_vector_image(A, B)

        @test isapprox(DifferentialPhaseContrast.differentiate_vectors(C), [0.0+0.0im 0.0-0.0im; 0.0-6.28319im 0.0-6.28319im]) atol=0.00001
    end

    @testset "scalar_multiply" begin
        # Test case 1: Real vectors
        vector1 = [1.0, 2.0]
        vector2 = [3.0, 4.0]
        @test DifferentialPhaseContrast.scalar_multiply(vector1, vector2) == 11.0

        # Test case 2: Complex vectors
        vector3 = [1.0 + 2.0im, 3.0 - 1.0im]
        vector4 = [2.0 - 1.0im, 4.0 + 3.0im]
        @test DifferentialPhaseContrast.scalar_multiply(vector3, vector4) == 9.0 - 8.0im

        # Test case 3: Mixed real and complex vectors
        vector5 = [1.0, 2.0]
        vector6 = [3.0 + 1.0im, 4.0 - 2.0im]
        @test DifferentialPhaseContrast.scalar_multiply(vector5, vector6) == 11.0 + 3.0im

        # Test case 4: Empty vectors
        vector7 = []
        vector8 = []
        @test_throws MethodError DifferentialPhaseContrast.scalar_multiply(vector7, vector8)
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