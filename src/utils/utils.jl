# sets subgradient (length n) of x to be 1 at j
function grad(x::SMCg{N,T},j,n) where {N,T}
  cc_grad = [i == j ? 1.0 : 0.0 for i=1:n]
  cv_grad = [i == j ? 1.0 : 0.0 for i=1:n]
  return SMCg{N,T}(x.cc,x.cv,cc_grad,cv_grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end
# sets subgradient (length n) of x to zero
function zgrad(x::SMCg{N,T},n) where {N,T}
  grad = @SVector zeros(n)
  return SMCg{N,T}(x.cc,x.cv,grad,grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end
function SMCg{N,T}(cc,cv,gp,Intv,IntvBox,xref) where {N,T}
  temp = SMCg{N,T}(cc,cv,[],[],Intv,x.cnst,IntvBox,xref)
  return grad(temp,gp[1],gp[2],IntvBox)
end

function convert(::Type{SMCg{N,T}},x::S) where {S<:Integer,N,T}
          seed = @SVector zeros(Float64,N)
          SMCg{N,T}(x,x,seed,seed,Interval(x),false,[Interval(0.0,1.0)],[1.0])
end
function convert(::Type{SMCg{N,T}},x::S) where {S<:AbstractFloat,N,T}
          seed = @SVector zeros(Float64,N)
          SMCg{N,T}(x,x,seed,seed,Interval(x),false,[Interval(0.0,1.0)],[1.0])
end
function convert(::Type{SMCg{N,T}},x::S) where {S<:Interval,N,T}
          seed = @SVector zeros(Float64,N)
          SMCg{N,T}(x.hi,x.lo,seed,seed,x,false,[Interval(0.0,1.0)],[1.0])
end
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Integer,N,T} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:AbstractFloat,N,T} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Interval,N,T} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Real,N,T} = SMCg{N,T}

# Defines functions smooth McCormick relaxations depend on
########### defines middle operator
@inline function mid3(x,y,z)
  (((x>=y)&&(y>=z))||((z>=y)&&(y>=x))) && (return y,2)
  (((y>=x)&&(x>=z))||((z>=x)&&(x>=y))) && (return x,1)
  return z,3
end
function mid_grad(cc_grad, cv_grad, id)
  if (id == 1)
    return cc_grad
  elseif (id == 2)
    return cv_grad
  elseif (id == 3)
    return zero(cc_grad)
  else
    error("Invalid mid3 position")
  end
end

# function for computing line segment between points (x1,y1) & (x2,y2)
function line_seg(x,x1,y1,x2,y2)
   if (x2-x1) == 0.0
     return y1
   else
     return y1*((x2-x)/(x2-x1)) + y2*((x-x1)/(x2-x1))
   end
end
function dline_seg(x,x1,y1,x2,y2)
    return (y2-y1)/(x2-x1)
end
function grad_calc(cv,cc,int1,int2,dcv,dcc)
  cv_grad = dcv*( int1==1 ? cv :( int1==2 ? cv : zeros(cv)))
  cc_grad = dcc*( int2==1 ? cc :( int2==2 ? cc : zeros(cv)))
  return cv_grad, cc_grad
end

"""
--------------------------------------------------------------------------------
Function: tighten_subgrad
--------------------------------------------------------------------------------
Description:
Tightens the interval bounds using subgradients.
--------------------------------------------------------------------------------
Inputs:
cc         Float64 - concave bound
cv         Float64 - convex bound
cc_grad    SVector{N,Float64} - subgradient/gradient of concave bound
cv_grad    SVector{N,Float64} - subgradient/gradient of convex bound
Xintv      Interval - Interval domain of function
Xbox       IntervalBox - Original decision variable bounds
xref       Vector{Float64} - Reference point in Xbox
--------------------------------------------------------------------------------
Returns:
The updated interval.
--------------------------------------------------------------------------------
"""
function tighten_subgrad(cc,cv,cc_grad,cv_grad,Xintv,Xbox,xref)
  if (length(Xbox)>0)
    upper_refine = cc
    lower_refine = cv
    for i=1:length(Xbox)
      upper_refine = upper_refine + cc_grad[i]*(Xbox[i]-xref[i])
      lower_refine = lower_refine + cv_grad[i]*(Xbox[i]-xref[i])
    end
    return Interval(max(lower_refine.lo,Xintv.lo),min(upper_refine.hi,Xintv.hi))
  else
    return Xintv
  end
end

"""
--------------------------------------------------------------------------------
Function: cut_bnds
--------------------------------------------------------------------------------
Description:
Cuts the relaxations using the interval bounds
--------------------------------------------------------------------------------
Inputs:
cc         Float64 - concave bound
cv         Float64 - convex bound
cc_grad    SVector{N,Float64} - subgradient/gradient of concave bound
cv_grad    SVector{N,Float64} - subgradient/gradient of convex bound
Intv       Interval - Interval domain of function
--------------------------------------------------------------------------------
Returns:
The updated tuple (cc,cv,cc_grad,cv_grad).
--------------------------------------------------------------------------------
"""
function cut_bnds(cc,cv,cc_grad,cv_grad,Intv)
  if (cc > Intv.hi)
    cc = Intv.hi
    cc_grad = zero(cc_grad)
  end
  if (cv < Intv.lo)
    cv = Intv.lo
    cv_grad = zero(cc_grad)
  end
  return cc,cv,cc_grad,cv_grad
end

function outer_rnd(Intv)
  return Interval(Intv.lo-MC_param.outer_param,Intv.hi+MC_param.outer_param)
end

function isequal(x,y,atol,rtol)
  dist = abs(x-y)
  avg = 0.5*abs(x+y)
  return (dist < (atol + avg*rtol))
end

function seed_g(n,j)
    temp = zeros(n)
    temp[j] = 1.0
    return SVector{n}(temp)
end
