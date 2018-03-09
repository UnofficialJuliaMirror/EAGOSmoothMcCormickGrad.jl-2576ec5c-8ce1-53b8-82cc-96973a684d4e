workspace()

using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays

# sets up initial options
opts = mc_opts()          # sets options for relaxtion
opts.kmax = 3             # sets number of iterations
opts.style = "KrawczykCW"   # sets style of contractor

#=
generates the expansion point parameters for the function using the opts
options
=#
f(x,p) = x[1]*p[1]+p[1]
g(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h(x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    return [t5]
end
hj(x,p) = [2*x[1]+p[1]]
P = [Interval(6.0,9.0)]
X = [Interval(-0.78,-0.4)]
p = [7.5]
pmid = mid.(P)
param = GenExpansionParams(h,hj,X,P,pmid,opts)

#=
relaxes the equality h(x,p)
=#
np = 1
szero = @SVector zeros(np)
sone = @SVector ones(np)
p_mc = [SMCg{np,Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
hbnds = MC_impRelax(h,hj,p_mc,pmid,X,P,opts,param)

#=
relaxation of f(x,p) & g(x,p) at (x(p),p)
=#
fbnds = impRelax_f(f,h,hj,X,P,p,pmid,opts,param)
fgbnds = impRelax_fg(f,g,h,hj,X,P,p,pmid,opts,param)
