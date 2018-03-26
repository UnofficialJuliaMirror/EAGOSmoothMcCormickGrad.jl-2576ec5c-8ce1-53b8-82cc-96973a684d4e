# convex/concave relaxation (Khan 3.1-3.2) of integer powers of 1/x for positive reals
function cv_negpowpos(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return x^n,n*x^(n-1)
end
function cc_negpowpos(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end
# convex/concave relaxation of integer powers of 1/x for negative reals
@inline function cv_negpowneg(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (isodd(n))
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    return x^n,n*x^(n-1)
  end
end
@inline function cc_negpowneg(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (isodd(n))
    return x^n,n*x^(n-1)
  else
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  end
end
# convex/concave relaxation of even powers greater than or equal to 4
@inline function cv_pow4(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return x^n,n*x^(n-1)
end
@inline function cc_pow4(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
end
# convex/concave relaxation of odd powers
@inline function cv_powodd(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
          if (xU<=zero(T))
             return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
          elseif (zero(T)<=xL)
            return x^n,n*x^(n-1)
          else
            val::T = (xL^n)*(xU-x)/(xU-xL)+(max(zero(T),x))^n
            dval::T = -(xL^n)/(xU-xL)+n*(max(zero(T),x))^(n-1)
            return val,dval
          end
end
@inline function cc_powodd(x::T,xL::T,xU::T,n::Integer) where {T<:AbstractFloat}
  if (xU<=zero(T))
    return x^n,n*x^(n-1)
  elseif (zero(T)<=xL)
    return line_seg(x,xL,xL^n,xU,xU^n),dline_seg(x,xL,xL^n,xU,xU^n,n*x^(n-1))
  else
    val::T = (xU^n)*(x-xL)/(xU-xL)+(min(zero(T),x))^n
    dval::T = (xU^n)/(xU-xL)+n*(min(zero(T),x))^(n-1)
    return val,dval
  end
end
# power of a generalized McCormick object raised to c
function pow(x::SMCg{N,T},c::Integer) where {N,T<:AbstractFloat}
  mid1::T = zero(T)
  mid2::T = zero(T)
  cc::T = zero(T)
  cv::T = zero(T)
  if (c==0)
    return one(T)
  elseif (c>0)
    if (c==2)
      return sqr(x)
    elseif (c==1)
      return x
    elseif (isodd(c))
      eps_max::T = x.Intv.hi
      eps_min::T = x.Intv.lo
      midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
      midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
      if (MC_param.mu >= 1)
        cc,dcc::T = cc_powodd(midcc,x.Intv.lo,x.Intv.hi,c)
        cv,dcv::T = cv_powodd(midcv,x.Intv.lo,x.Intv.hi,c)
        gcc1::T,gdcc1::T = cc_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcv1::T,gdcv1::T = cv_powodd(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcc2::T,gdcc2::T = cc_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  gcv2::T,gdcv2::T = cv_powodd(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    	  cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
      else
        cc,dcc = cc_powodd(midcc,x.Intv.lo,x.Intv.hi,c)
        cv,dcv = cv_powodd(midcv,x.Intv.lo,x.Intv.hi,c)
        cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
        cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      end
      intv::Interval{T} = pow(x.Intv,c)
      return SMCg{N,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
    else
      if (x.Intv.hi<zero(x.cc))
        eps_min = x.Intv.hi
        eps_max = x.Intv.lo
      elseif (x.Intv.lo>zero(x.cc))
        eps_min = x.Intv.lo
        eps_max = x.Intv.hi
      else
        eps_min = zero(x.cc)
        eps_max = (abs(x.Intv.lo)>=abs(x.Intv.hi))? x.Intv.lo : x.Intv.hi
      end
      midcc,cc_id = mid3(x.cc,x.cv,eps_max)
      midcv,cv_id = mid3(x.cc,x.cv,eps_min)
      cc,dcc = cc_pow4(midcc,x.Intv.lo,x.Intv.hi,c)
      cv,dcv = cv_pow4(midcv,x.Intv.lo,x.Intv.hi,c)
      if (MC_param.mu >= 1)
        gcc1,gdcc1 = cc_pow4(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcv1,gdcv1 = cv_pow4(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcc2,gdcc2 = cc_pow4(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  gcv2,gdcv2 = cv_pow4(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    	  cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
      else
        cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
        cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      end
      intv = x.Intv^c
      return SMCg{N,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
    end
  else
    if (x.Intv.hi<zero(x.cc))
      if (isodd(c))
        eps_min = x.Intv.hi
        eps_max = x.Intv.lo
      else
        eps_min = x.Intv.lo
        eps_max = x.Intv.hi
      end
      midcc,cc_id = mid3(x.cc,x.cv,eps_max)
      midcv,cv_id = mid3(x.cc,x.cv,eps_min)
      cc,dcc = cc_negpowneg(midcc,x.Intv.lo,x.Intv.hi,c)
      cv,dcv = cv_negpowneg(midcv,x.Intv.lo,x.Intv.hi,c)
      if (MC_param.mu >= 1)
        gcc1,gdcc1 = cc_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcv1,gdcv1 = cv_negpowneg(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcc2,gdcc2 = cc_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  gcv2,gdcv2 = cv_negpowneg(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    	  cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
      else
        cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
        cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      end
      intv = x.Intv^c
      return SMCg{N,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
    elseif (x.Intv.lo>zero(x.cc))
      eps_min = x.Intv.hi
      eps_max = x.Intv.lo
      midcc,cc_id = mid3(x.cc,x.cv,eps_max)
      midcv,cv_id = mid3(x.cc,x.cv,eps_min)
      cc,dcc = cc_negpowpos(midcc,x.Intv.lo,x.Intv.hi,c)
      cv,dcv = cv_negpowpos(midcv,x.Intv.lo,x.Intv.hi,c)
      if (MC_param.mu >= 1)
        gcc1,gdcc1 = cc_negpowpos(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcv1,gdcv1 = cv_negpowpos(x.cv,x.Intv.lo,x.Intv.hi,c)
    	  gcc2,gdcc2 = cc_negpowpos(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  gcv2,gdcv2 = cv_negpowpos(x.cc,x.Intv.lo,x.Intv.hi,c)
    	  cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
    	  cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
      else
        cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
        cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
      end
      intv = x.Intv^c
      return SMCg{N,T}(cc, cv, cc_grad, cv_grad, intv,x.cnst,x.IntvBox,x.xref)
    else
      error("Function unbounded on domain")
    end
  end
end
########### power of a generalized McCormick object raised to c
function  (^)(x::SMCg{N,T},c::Integer) where {N,T<:AbstractFloat}
  pow(x,c)
end
function  (^)(x::SMCg{N,T},c::Float64) where {N,T<:AbstractFloat}
  if (isinteger(c))
    return pow(x,Int64(c))
  else
    error("Relaxation of non-integer powers are unavailable.")
  end
end
########### Defines inverse
function inv(x::SMCg{N,T}) where {N,T<:AbstractFloat}
  #println("ran inverse")
  x^(-1)
end
