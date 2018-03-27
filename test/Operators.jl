module Operators

using Compat
using Compat.Test
using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays


function about(calc,val,tol)
    return (val - tol <= calc <= val + tol)
end

################################################################################
######################## Tests Standard McCormick Relaxations ##################
################################################################################
EAGOSmoothMcCormickGrad.set_diff_relax(0)

######## tests division of same object ######
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(1.0,5.0);Interval(-1.0,2.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(3.0,3.0,a,a,xIBox[1],false,xIBox,mBox)
out = X/X
@test out == 1.0

######## tests nonsmooth division ######
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(-5.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = X/Y
@test about(out.cc,0.6,1E-6)
@test about(out.cv,0.41666666,1E-6)
@test about(out.cc_grad[1],-0.2,1E-6)
@test about(out.cc_grad[2],0.2,1E-6)
@test about(out.cv_grad[1],-0.333333,1E-6)
@test about(out.cv_grad[2],0.1875,1E-6)
@test about(out.Intv.lo,-1.33333333,1E-6)
@test about(out.Intv.hi,1.0,1E-6)

######## tests exponent on product ######
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(-5.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = exp(X*Y)
@test about(out.cc,2.708614394334035e6,1E-1)
@test about(out.cv,1096.6331584284585,1E-1)
@test about(out.cc_grad[1],-2.80201e5,1E1)
@test about(out.cc_grad[2],-2.80201e5,1E1)
@test about(out.cv_grad[1],-5483.17,1E-1)
@test about(out.cv_grad[2],-3289.9,1E-1)
@test about(out.Intv.lo,2.06115e-09,1E-1)
@test about(out.Intv.hi,3.26902e+06,1E1)

######## tests log on product ######
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(3.0,7.0);Interval(3.0,9.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(7.0,7.0,b,b,xIBox[2],false,xIBox,mBox)
out = log(X*Y)
@test about(out.cc,3.4011973816621555,1E-8)
@test about(out.cv,2.7377551742960287,1E-8)
@test about(out.cc_grad[1],0.3,1E-1)
@test about(out.cc_grad[2],0.1,1E-1)
@test about(out.cv_grad[1],0.108106,1E-5)
@test about(out.cv_grad[2],0.108106,1E-5)
@test about(out.Intv.lo,2.19722,1E-4)
@test about(out.Intv.hi,4.14314,1E-4)

# tests powers (square)
a = seed_g(Float64,1,2)
xIBox = [Interval(3.0,7.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
out = X^2
@test about(out.cc,19,1E-8)
@test about(out.cv,16,1E-8)
@test about(out.cc_grad[1],10.0,1E-1)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],8.0,1E-5)
@test about(out.cv_grad[2],0.0,1E-5)
@test about(out.Intv.lo,9,1E-4)
@test about(out.Intv.hi,49,1E-4)

# tests powers (^-2 on positive domain)
a = seed_g(Float64,1,2)
xIBox = [Interval(3.0,7.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
out = X^(-2)
@test about(out.cc,0.08843537414965986,1E-8)
@test about(out.cv,0.0625,1E-8)
@test about(out.cc_grad[1],-0.0226757,1E-1)
@test about(out.cc_grad[2],0.0,1E-4)
@test about(out.cv_grad[1],-0.03125,1E-5)
@test about(out.cv_grad[2],0.0,1E-4)
@test about(out.Intv.lo,0.0204081,1E-4)
@test about(out.Intv.hi,0.111112,1E-4)

# tests powers (^-2 on negative domain)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(-2)
@test about(out.cc,0.08246527777777776,1E-8)
@test about(out.cv,0.04938271604938271,1E-8)
@test about(out.cc_grad[1],0.0190972,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],0.0219479,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,0.015625,1E-4)
@test about(out.Intv.hi,0.111112,1E-4)

# tests powers (^0)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(0)
@test out == 1.0

# tests powers (^1)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(1)
@test about(out.cc,-4.5,1E-8)
@test about(out.cv,-4.5,1E-8)
@test about(out.cc_grad[1],1.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],1.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,-8.0,1E-4)
@test about(out.Intv.hi,-3.0,1E-4)

# tests powers (^2)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(2)
@test about(out.cc,25.5,1E-8)
@test about(out.cv,20.25,1E-8)
@test about(out.cc_grad[1],-11.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],-9.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,9.0,1E-4)
@test about(out.Intv.hi,64.0,1E-4)

# tests powers (^3)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(3)
@test about(out.cc,-91.125,1E-8)
@test about(out.cv,-172.5,1E-8)
@test about(out.cc_grad[1],60.75,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],97.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,-512,1E-4)
@test about(out.Intv.hi,-27,1E-4)

# tests powers (^4)
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(4)
@test about(out.cc,1285.5,1E-8)
@test about(out.cv,410.0625,1E-8)
@test about(out.cc_grad[1],-803.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],-364.5,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,81,1E-4)
@test about(out.Intv.hi,4096,1E-4)

# tests sqrt
a = seed_g(Float64,1,2)
xIBox = [Interval(3.0,9.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = sqrt(X)
@test about(out.cc,2.1213203435596424,1E-8)
@test about(out.cv,2.049038105676658,1E-8)
@test about(out.cc_grad[1],0.235702,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],0.211325,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,1.73205,1E-4)
@test about(out.Intv.hi,3,1E-4)

# tests powers (^3 greater than zero ISSUE WITH CC)
a = seed_g(Float64,1,2)
xIBox = [Interval(3.0,8.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(3)
@test about(out.cc,172.5,1E-8)
@test about(out.cv,91.125,1E-8)
@test about(out.cc_grad[1],97.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],60.75,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,27,1E-4)
@test about(out.Intv.hi,512,1E-4)

# tests powers (^4 greater than zero)

a = seed_g(Float64,1,2)
xIBox = [Interval(3.0,8.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(4)
@test about(out.cc,1285.5,1E-1)
@test about(out.cv,410.0625,1E-1)
@test about(out.cc_grad[1],803.0,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],364.5,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,81,1E-4)
@test about(out.Intv.hi,4096,1E-4)


# tests powers (^4 zero in range)

a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out = X^(4)
@test about(out.cc,2818.5,1)
@test about(out.cv,410.0625,1)
@test about(out.cc_grad[1],-365.0,1)
@test about(out.cc_grad[2],0.0,1)
@test about(out.cv_grad[1],-364.5,1)
@test about(out.cv_grad[2],0.0,1)
@test about(out.Intv.lo,0,1)
@test about(out.Intv.hi,4096,1)

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(-5.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = exp2(X*Y)
@test about(out.cc,27150.62857159206,1E-1)
@test about(out.cv,128.0,1E-1)
@test about(out.cc_grad[1],-2808.69,1E1)
@test about(out.cc_grad[2],-2808.69,1E1)
@test about(out.cv_grad[1],-443.614,1E-1)
@test about(out.cv_grad[2],-266.169,1E-1)
@test about(out.Intv.lo,9.53674e-07,1E-1)
@test about(out.Intv.hi,32768,1E1)

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(-5.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = exp10(X*Y)
@test about(out.cc,8.285714285714286e14,1E-1)
@test about(out.cv,1.0e7,1E-1)
@test about(out.cc_grad[1],-8.57143e13,1E8)
@test about(out.cc_grad[2],-8.57143e13,1E8)
@test about(out.cv_grad[1],-1.15129e8,1E3)
@test about(out.cv_grad[2],-6.90776e7,1E3)
@test about(out.Intv.lo,9.99999e-21,1E-1)
@test about(out.Intv.hi,1e+15,1E1)

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(3.0,7.0);Interval(3.0,9.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(7.0,7.0,b,b,xIBox[2],false,xIBox,mBox)
out = log2(X*Y)
@test about(out.cc,4.906890595608519,1E-7)
@test about(out.cv,3.94974581312498,1E-7)
@test about(out.cc_grad[1],0.432809,1E-5)
@test about(out.cc_grad[2],0.14427,1E-5)
@test about(out.cv_grad[1],0.155964,1E-5)
@test about(out.cv_grad[2],0.155964,1E-5)
@test about(out.Intv.lo,3.16992,1E-5)
@test about(out.Intv.hi,5.97728,1E-5)

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(3.0,7.0);Interval(3.0,9.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.0,4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(7.0,7.0,b,b,xIBox[2],false,xIBox,mBox)
out = log10(X*Y)
@test about(out.cc,1.4771212547196624,1E-7)
@test about(out.cv,1.1889919649988407,1E-7)
@test about(out.cc_grad[1],0.130288,1E-5)
@test about(out.cc_grad[2],0.0434294,1E-5)
@test about(out.cv_grad[1],0.0469499,1E-5)
@test about(out.cv_grad[2],0.0469499,1E-5)
@test about(out.Intv.lo,0.954242,1E-5)
@test about(out.Intv.hi,1.79935,1E-5)

a = seed_g(Float64,1,2)
xIBox = [Interval(-3.0,8.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(4.5,4.5,a,a,xIBox[1],false,xIBox,mBox)
out = abs(X)
@test about(out.cc,6.409090909090908,1E-1)
@test about(out.cv,4.5,1E-1)
@test about(out.cc_grad[1],0.454545,1E-4)
@test about(out.cc_grad[2],0.0,1E-1)
@test about(out.cv_grad[1],1.0,1E-4)
@test about(out.cv_grad[2],0.0,1E-1)
@test about(out.Intv.lo,0,1E-4)
@test about(out.Intv.hi,8,1E-4)

end

# tests powers (^3 zero in range, ISSUE) TO DO
#=
a = seed_g(Float64,1,2)
xIBox = [Interval(-8.0,3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.5,-4.5,a,a,xIBox[1],false,xIBox,mBox)
out2 = X^(3)
ans1 = about(out.cc,-91.125,1E-8)
ans2 = about(out.cv,-172.5,1E-8)
ans3 = about(out.cc_grad[1],60.75,1E-4)
ans4 = about(out.cc_grad[2],0.0,1E-1)
ans5 = about(out.cv_grad[1],97.0,1E-4)
ans6 = about(out.cv_grad[2],0.0,1E-1)
ans7 = about(out.Intv.lo,-512,1E-4)
ans8 = about(out.Intv.hi,-27,1E-4)
=#
