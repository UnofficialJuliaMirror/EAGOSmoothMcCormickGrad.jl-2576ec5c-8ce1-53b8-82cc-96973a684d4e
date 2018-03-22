#module ImplicitBnd

#using Compat
#using Compat.Test
using IntervalArithmetic
using EAGOSmoothMcCormickGrad

nx = 3
np = 2
opt = Any[nx,np,0.5]
epsv = 1E-8

seed1a = 0.5*seed_g(Float64,1,2)
seed1b = 0.75*seed_g(Float64,1,2)
seed2a = 0.5*seed_g(Float64,2,2)
seed2b = 0.75*seed_g(Float64,2,2)

#=
x = [SMCg{np,Float64}()
     SMCg{np,Float64}()
     SMCg{np,Float64}()]

p = [SMCg{np,Float64}()
     SMCg{np,Float64}()]

p_ref = [SMCg{np,Float64}()
         SMCg{np,Float64}()]

xa = [SMCg{np,Float64}()
      SMCg{np,Float64}()
      SMCg{np,Float64}()]

xA = [SMCg{np,Float64}()
      SMCg{np,Float64}()
      SMCg{np,Float64}()]

z = [SMCg{np,Float64}()
     SMCg{np,Float64}()
     SMCg{np,Float64}()]

Affine_Exp!(x,p,p_ref,xa,xA,z,opt)
xa
xA
z

Correct_Exp!(z_mc,x_mc,X,nx,np,epsv)
x_mc
=#
#end
