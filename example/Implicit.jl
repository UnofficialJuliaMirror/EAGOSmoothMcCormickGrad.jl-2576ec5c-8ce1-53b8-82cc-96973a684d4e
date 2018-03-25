#workspace()

using EAGOSmoothMcCormickGrad
using IntervalArithmetic
using StaticArrays

# sets up initial options
opts = mc_opts(Float64)          # sets options for relaxtion
opts.kmax = 3             # sets number of iterations
#opts.style = "NewtonGS"   # sets style of contractor
opts.style = "KrawczykCW"   # sets style of contractor

#=
generates the expansion point parameters for the function using the opts
options using inverse preconditioner
=#
f(x,p) = x[1]*p[1]+p[1]
g(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h1(x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    return [t5]
end
hj1(x,p) = [2*x[1]+p[1]]
P = [Interval(6.0,9.0)]
X = [Interval(-0.78,-0.4)]
p = [7.5]
pmid = mid.(P)
param = GenExpansionParams(h1,hj1,X,P,pmid,opts)

#=
relaxes the equality h(x,p)

np = 1
szero = @SVector zeros(np)
sone = @SVector ones(np)
p_mc = [SMCg{np,Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
hbnds = MC_impRelax(h,hj,p_mc,pmid,X,P,opts,param)


relaxation of f(x,p) & g(x,p) at (x(p),p)

fbnds = impRelax_f(f,h,hj,X,P,p,pmid,opts,param)
fgbnds = impRelax_fg(f,g,h,hj,X,P,p,pmid,opts,param)




generates the expansion point parameters for the function using the opts
options using in place LDU full pivot conditioner with sparse calc

f1(x,p) = x[1]*p[1]+p[1]
g1(x,p) = [x[1]*p[1]+p[1];
          x[1]*p[1]+2*p[1]]
function h1!(hout,x,p)
    t1 = x[1]^2
    t2 = x[1]*p[1]
    t3 = 4.0
    t4 = t1 + t2
    t5 = t2 + t3
    hout[1] = t5
end
function hj1!(hout,x,p)
    hout[1] = 2*x[1]+p[1]
end
P1 = [Interval(6.0,9.0)]
X1 = [Interval(-0.78,-0.4)]
p1 = [7.5]
pmid1 = mid.(P1)
param1 = InGenExpansionParams(h1!,hj1!,X1,P1,pmid1,opts)



relaxes the equality h(x,p)

np1 = 1
szero1 = @SVector zeros(np1)
sone1 = @SVector ones(np1)
p_mc1 = [SMCg{np,Float64}(p1[i],p1[i],sone1,sone1,@interval(P1[i].lo,P1[i].hi),false,[∅],[1.0]) for i=1:np1]
hbnds = MC_NimpRelax(h1!,hj1!,p_mc1,pmid1,X1,P1,opts,param1)


relaxation of f(x,p) & g(x,p) at (x(p),p)
fbnds = NimpRelax_f(f1,h1!,hj1!,X1,P1,p,pmid1,opts,param1)
fgbnds = NimpRelax_fg(f1,g1,h1!,hj1!,X1,P1,p,pmid1,opts,param1)
=#
