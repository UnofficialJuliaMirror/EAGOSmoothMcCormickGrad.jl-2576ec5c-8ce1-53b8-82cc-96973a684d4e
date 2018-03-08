workspace()

#module Operators

#using Compat
#using Compat.Test
using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays

SMCglist_U = []
SMCglist_B = []
SMCglist_BC = []

# create seed gradient
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)

# create SmoothMcCormick seed object for x1 = 2.0 on [1.0,3.0] for relaxing
# a function f(x1,x2) on the interval box xIbox using mBox as a reference point
x = 0.6
xIntv1 = Interval(0.5,0.75)
xIBox = [xIntv1;xIntv1]
mBox = mid.(xIBox)
SMCg1 = SMCg{2,Float64}(x,x,a,a,xIntv1,false,xIBox,mBox)


# Test binary operators
binary_list = [+,-,/,*,min,max]
for i = 1:length(binary_list)
    push!(SMCglist_B,(binary_list[i])(SMCg1,SMCg1))
end

# Test unary operators
unary_list1 = [-, sin, cos, tan, asin,
             acos, atan, sinh, cosh, tanh,
             asinh, atanh, exp, exp2,
             exp10, log, log2, log10, sqrt,
             step,sign]
for i = 1:length(unary_list1)
    push!(SMCglist_U,(unary_list1[i])(SMCg1))
end

x = 2.6
xIntv1 = Interval(2.5,2.75)
xIBox = [xIntv1;xIntv1]
mBox = mid.(xIBox)
SMCg2 = SMCg{2,Float64}(x,x,a,a,xIntv1,false,xIBox,mBox)

unary_list2 = [acosh, exp, exp2,
             exp10, log, log2, log10, sqrt]
for i = 1:length(unary_list2)
    push!(SMCglist_U,(unary_list2[i])(SMCg2))
end

expo_list = Any[-3,-2,-1,0,1,2,3,4,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0]
for i = 1:length(expo_list)
    push!(SMCglist_BC,SMCg2^expo_list[i])
end

#end
