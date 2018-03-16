#workspace()

using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays
using DataFrames
using CSV

# Plots the differentiable McCormick relaxations and
# original function x*y over [0,200] by [0,400], respectively
EAGOSmoothMcCormickGrad.set_diff_relax(1)

Vals = zeros(Float64,40*80,9)

count = 1
for i=1:20:200
    for j=1:20:400
        temp1 = SMCg{2,Float64}(Float64(j-1),Float64(j-1),seed_g(Float64,1,2),
                                 seed_g(Float64,1,2),Interval(0.0,400.0),
                                 false, [Interval(0.0,400.0),Interval(0.0,200.0)],[200.0,100.0])
        temp2 = SMCg{2,Float64}(Float64(i-1),Float64(i-1),seed_g(Float64,2,2),
                                 seed_g(Float64,2,2),Interval(0.0,400.0),
                                 false, [Interval(0.0,400.0),Interval(0.0,200.0)],[200.0,100.0])
        Vals[count,1] = j-1 # x
        Vals[count,2] = i-1 # y
        Vals[count,3] = (j-1)*(i-1) #x*y
        temp_MC = temp1*temp2 # relax x*y
        Vals[count,4] = temp_MC.cc
        Vals[count,5] = temp_MC.cv
        Vals[count,6] = temp_MC.cc_grad[1]
        Vals[count,7] = temp_MC.cc_grad[2]
        Vals[count,8] = temp_MC.cv_grad[1]
        Vals[count,9] = temp_MC.cv_grad[2]
        count = count + 1
        println("dcvdx: $(temp_MC.cv_grad[1])")
        println("dcvdy: $(temp_MC.cv_grad[2])")
        println("i: $i")
        println("j: $j")
    end
end
DF = convert(DataFrame,Vals)
CSV.write("/home/mewilhel/Desktop/MultiplicationPlot3.csv", DF)
