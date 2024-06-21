using RGRDMT
using Test

@testset "Core Features" begin
    include("Core.jl")
end

@testset "Loading MPS from MATLAB Code" begin
    include("loader.jl")
end
