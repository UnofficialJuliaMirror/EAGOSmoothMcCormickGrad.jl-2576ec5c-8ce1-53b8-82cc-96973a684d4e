using EAGOSmoothMcCormickGrad
using Base.Test

# write your own tests here
println("Testing Utility Functions...")
t = @elapsed include("Utilities.jl")
println("done (took $t seconds).")

#=
println("Testing McCormick Operators...")
t = @elapsed include("Operators.jl")
println("done (took $t seconds).")
=#

#=
println("Implicit Bounding Utilities...")
t = @elapsed include("D1_Interval_Test.jl")
println("done (took $t seconds).")

println("Implicit Function Routines...")
t = @elapsed include("D1_Interval_Test.jl")
println("done (took $t seconds).")
=#
