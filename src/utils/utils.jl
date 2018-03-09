"""
    seed_g(T::Type,j::Int64,N::Int64)

Creates a `x::SVector{N,T}` object that is one at `x[j]` and zero everywhere else.
"""
function seed_g(T::Type,j::Int64,N::Int64)
    return SVector{N,T}([i == j ? 1.0 : 0.0 for i=1:N])
end

"""
    grad(x::SMCg{N,T},j,n) where {N,T}

sets convex and concave (sub)gradients of length `n` of `x` to be `1` at index `j`
"""
function grad(x::SMCg{N,T},j::Int64) where {N,T<:AbstractFloat}
  sv_grad::SVector{N,T} = seed_g(T,j,N)
  return SMCg{N,T}(x.cc,x.cv,sv_grad,sv_grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end

"""
    zgrad(x::SMCg{N,T},n::Int64) where {N,T}

sets convex and concave (sub)gradients of length `n` to be zero
"""
function zgrad(x::SMCg{N,T}) where {N,T<:AbstractFloat}
  grad::SVector{N,T} = @SVector zeros(N)
  return SMCg{N,T}(x.cc,x.cv,grad,grad,x.Intv,x.cnst,x.IntvBox,x.xref)
end

function convert(::Type{SMCg{N,T}},x::S) where {S<:Integer,N,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,T}(convert(T,x),convert(T,x),seed,seed,Interval(convert(Interval{T},x)),
                    false,[emptyinterval(T)],[zero(T)])
end
function convert(::Type{SMCg{N,T}},x::S) where {S<:AbstractFloat,N,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,T}(convert(T,x),convert(T,x),seed,seed,Interval(convert(Interval{T},x)),
                    false,[emptyinterval(T)],[zero(T)])
end
function convert(::Type{SMCg{N,T}},x::S) where {S<:Interval,N,T<:AbstractFloat}
          seed::SVector{N,T} = @SVector zeros(T,N)
          SMCg{N,T}(convert(T,x.hi),convert(T,x.lo),seed,seed,convert(Interval{T},x),
                    false,[emptyinterval(T)],[zero(T)])
end

promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Integer,N,T<:AbstractFloat} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:AbstractFloat,N,T<:AbstractFloat} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Interval,N,T<:AbstractFloat} = SMCg{N,T}
promote_rule(::Type{SMCg{N,T}}, ::Type{S}) where {S<:Real,N,T<:AbstractFloat} = SMCg{N,T}

"""
    mid3(x::T,y::T,z::T)

Calculates the midpoint of three numbers returning the value and the index.
"""
function mid3(x::T,y::T,z::T) where {T<:AbstractFloat}
  (((x>=y)&&(y>=z))||((z>=y)&&(y>=x))) && (return y,2)
  (((y>=x)&&(x>=z))||((z>=x)&&(x>=y))) && (return x,1)
  return z,3
end

"""
    mid_grad(cc_grad::SVector{N,T}, cv_grad::SVector{N,T}, id::Int64)

Takes the concave relaxation gradient 'cc_grad', the convex relaxation gradient
'cv_grad', and the index of the midpoint returned 'id' and outputs the appropriate
gradient according to McCormick relaxation rules.
"""
function mid_grad(cc_grad::SVector{N,T}, cv_grad::SVector{N,T}, id::Int64) where {N,T<:AbstractFloat}
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

"""
    line_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the line segment between `(x1,y1)` and `(x2,y2)` at `x = x0`.
"""
function line_seg(x0::T,x1::T,y1::T,x2::T,y2::T) where {T<:AbstractFloat}
   if (x2-x1) == zero(T)
     return y1
   else
     return y1*((x2-x0)/(x2-x1)) + y2*((x0-x1)/(x2-x1))
   end
end

"""
    dline_seg(x0::T,x1::T,y1::T,x2::T,y2::T)

Calculates the value of the slope line segment between `(x1,y1)` and `(x2,y2)`
defaults to evaluating the derivative of the function if the interval is tight.
"""
function dline_seg(x0::T,x1::T,y1::T,x2::T,y2::T,d::T) where {T<:AbstractFloat}
    if (x2 == x1)
      return d
    else
      return (y2-y1)/(x2-x1)
    end
end

"""
    grad_calc(cv::T,cc::T,int1::Int64,int2::Int64,dcv::SVector{N,T},dcc::SVector{N,T}) where {N,T}

(Sub)gradient calculation function. Takes the convex gradient, 'cv', the
concave gradient, 'cc', the mid index values 'int1,int2', and the derivative of
the convex and concave envelope functions 'dcv,dcc'.
"""
function grad_calc(cv::SVector{N,T},cc::SVector{N,T},int1::Int64,int2::Int64,dcv::T,dcc::T) where {N,T<:AbstractFloat}
  cv_grad::SVector{N,T} = dcv*( int1==1 ? cv :( int1==2 ? cv : zeros(cv)))
  cc_grad::SVector{N,T} = dcc*( int2==1 ? cc :( int2==2 ? cc : zeros(cv)))
  return cv_grad, cc_grad
end

"""
    tighten_subgrad(cc,cv,cc_grad,cv_grad,Xintv,Xbox,xref)

Tightens the interval bounds using subgradients. Inputs:
* `cc::T`: concave bound
* `cv::T`: convex bound
* `cc_grad::SVector{N,T}`: subgradient/gradient of concave bound
* `cv_grad::SVector{N,T}`: subgradient/gradient of convex bound
* `Xintv::Interval{T}`: Interval domain of function
* `Xbox::Vector{Interval{T}}`: Original decision variable bounds
* `xref::Vector{T}`: Reference point in Xbox
"""
function tighten_subgrad(cc::T,cv::T,cc_grad::SVector{N,T},cv_grad::SVector{N,T},
                         Xintv::Interval{T},Xbox::Vector{Interval{T}},xref::Vector{T}) where {N,T<:AbstractFloat}
  if (length(Xbox)>0 && Xbox[1]!=∅)
    upper_refine::Interval{T} = Interval{T}(cc)
    lower_refine::Interval{T} = Interval{T}(cv)
    for i=1:length(Xbox)
      upper_refine = upper_refine + cc_grad[i]*(Xbox[i]-xref[i])
      lower_refine = lower_refine + cv_grad[i]*(Xbox[i]-xref[i])
    end
    return Interval{T}(max(lower_refine.lo,Xintv.lo),min(upper_refine.hi,Xintv.hi))
  else
    return Xintv
  end
end

"""
    outer_rnd!(Intv::Interval{T})

Outer rounds the interval `Intv` by `MC_param.outer_param`.
"""
function outer_rnd(Intv::Interval{T}) where {T<:AbstractFloat}
  return Interval{T}(Intv.lo-MC_param.outer_param, Intv.hi+MC_param.outer_param)
end

"""
    isequal(x::S,y::S,atol::T,rtol::T)

Checks that `x` and `y` are equal to with absolute tolerance `atol` and relative
tolerance `rtol`.
"""
function isequal(x::S,y::S,atol::T,rtol::T) where {S,T<:AbstractFloat}
  return (abs(x-y) < (atol + 0.5*abs(x+y)*rtol))
end

"""
    Intv(x::SMCg{N,T})
"""
function Intv(x::SMCg{N,T}) where {N,T<:AbstractFloat}
  return x.Intv
end
