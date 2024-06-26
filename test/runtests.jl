using RGRDMT
using Test

@testset "Loading MPS from MATLAB Code" begin
    include("loader.jl")
end

@testset "Constructing Isometries from MPS" begin
    include("isometry.jl")
end

