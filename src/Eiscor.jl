module Eiscor

# list of support real types
SupportedTypes = [Float16,Float32,Float64,BigFloat]

# submodules
include("RootFreeUnitary.jl")
include("Rotation.jl")

end
