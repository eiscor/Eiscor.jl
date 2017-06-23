module Eiscor

# list of support real types
SupportedTypes = [Float16,Float32,Float64,BigFloat]

# submodules
include("Rotation.jl")
include("RootFreeUnitary.jl")

end
