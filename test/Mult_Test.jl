module Mult_Test

using Compat
using Compat.Test
using IntervalArithmetic
using EAGOSmoothMcCormickGrad

################################################################################
################################################################################
##############      Testing for Nonsmooth Standard Mult           ##############
################################################################################
################################################################################
EAGOSmoothMcCormickGrad.set_diff_relax(0)

################################################################################
################### Test Nonsmooth Zero in Both Case (Failing)   ###############
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-2.0,1.0);Interval(-1.0,2.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(0.0,0.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(1.0,1.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

println("out: $out")

@test out.cc == 2.0
@test out.cv == -1.0
@test out.cc_grad[1] == 2.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == 2.0
@test out.cv_grad[2] == 1.0
@test out.Intv.lo == -4.0
@test out.Intv.hi == 2.0


################################################################################
###################### Test Nonsmooth X1.l>0   (Passing)  ######################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(1.0,5.0);Interval(-1.0,2.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(3.0,3.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(1.0,1.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 5.0
@test out.cv == 1.0
@test out.cc_grad[1] == 2.0
@test out.cc_grad[2] == 1.0
@test out.cv_grad[1] == 2.0
@test out.cv_grad[2] == 5.0
@test out.Intv.lo == -5.0
@test out.Intv.hi == 10.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  X2.l>0 (Passing)  ######################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-6.0,-2.0);Interval(1.0,3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(2.0,2.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == -6.0
@test out.cv == -10.0
@test out.cc_grad[1] == 1.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == 3.0
@test out.cv_grad[2] == -2.0
@test out.Intv.lo == -18.0
@test out.Intv.hi == -2.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  X2.h<0 (Passing)  ######################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-6.0,-2.0);Interval(-7.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-5.0,-5.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 24.0
@test out.cv == 16.0
@test out.cc_grad[1] == -7.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == -3.0
@test out.cv_grad[2] == -2.0
@test out.Intv.lo == 6.0
@test out.Intv.hi == 42.0

################################################################################
############## Test Nonsmooth X1.h<0  &&  0 in X2 (Passing)  ###################
################################################################################

a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-6.0,-2.0);Interval(-7.0,4.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-4.0,-4.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-5.0,-5.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 24.0
@test out.cv == 16.0
@test out.cc_grad[1] == -7.0
@test out.cc_grad[2] == -2.0
@test out.cv_grad[1] == -7.0
@test out.cv_grad[2] == -6.0
@test out.Intv.lo == -24.0
@test out.Intv.hi == 42.0

################################################################################
############## Test Nonsmooth 0 in X1  &&  X2.l > 0 ()  #################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(1.0,4.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(3.0,3.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == -5.0
@test out.cv == -8.0
@test out.cc_grad[1] == 4.0
@test out.cc_grad[2] == -3.0
@test out.cv_grad[1] == 1.0
@test out.cv_grad[2] == -3.0
@test out.Intv.lo == -12.0
@test out.Intv.hi == 16.0

################################################################################
############## Test Nonsmooth 0 in X1  &&  X2.h < 0 ()         #################
################################################################################
a = seed_g(Float64,1,2)
b = seed_g(Float64,2,2)
xIBox = [Interval(-3.0,4.0);Interval(-5.0,-3.0)]
mBox = mid.(xIBox)
X = SMCg{2,Float64}(-2.0,-2.0,a,a,xIBox[1],false,xIBox,mBox)
Y = SMCg{2,Float64}(-4.0,-4.0,b,b,xIBox[2],false,xIBox,mBox)
out = X*Y

@test out.cc == 9.0
@test out.cv == 7.0
@test out.cc_grad[1] == -3.0
@test out.cc_grad[2] == -3.0
@test out.cv_grad[1] == -5.0
@test out.cv_grad[2] == -3.0
@test out.Intv.lo == -20.0
@test out.Intv.hi == 15.0
end
