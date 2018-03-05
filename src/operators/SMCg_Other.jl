########### Defines differentiable step relaxations
@inline function cv_step(x,xL,xU)
  if (xU<=zero(xU))
    return zero(x),zero(x)
  elseif (xL>=zero(xL))
    return one(x),zero(x)
  elseif (x>=zero(x))
    return (x/xU)^2.0,2.0*x/xU^2.0
  else
    return zero(x),zero(x)
  end
end
@inline function cc_step(x,xL,xU)
  if (xU<=zero(xU))
     return zero(x),zero(x)
  elseif (xL>=zero(xL))
     return one(x),zero(x)
  elseif (x>=zero(x))
     return one(x),zero(x)
  else
    return one(x)-(x/xL)^2.0,-2.0*x/xL^2.0
  end
end
@inline function cv_step_NS(x,xL,xU)
  if (xU<=zero(xU))
    return zero(x),zero(x)
  elseif (xL>=zero(xL))
    return one(x),zero(x)
  elseif (x>=zero(x))
    return line_seg(x,0,0,xU,1),dline_seg(x,0,0,xU,1)
  else
    return zero(x),zero(x)
  end
end
@inline function cc_step_NS(x,xL,xU)
  if (xU<=zero(xU))
     return zero(x),zero(x)
  elseif (xL>=zero(xL))
     return one(x),zero(x)
  elseif (x>=zero(x))
     return one(x),zero(x)
  else
    return line_seg(x,xL,0,0,1),dline_seg(x,xL,0,0,1)
  end
end
@inline function step(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc,dcc = cc_step(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_step(midcv,x.Intv.lo,x.Intv.hi)
    gcc1,gdcc1 = cc_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = cv_step(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = cc_step(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = cv_step(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc,dcc = cc_step_NS(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = cv_step_NS(midcv,x.Intv.lo,x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, step(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline function step(x::Interval)
           isempty(x) && return emptyinterval(x)
           xmin = ((x.lo)<zero(x.lo)) ? zero(x.lo) : one(x.lo)
           xmax = ((x.hi)>=zero(x.lo)) ? one(x.lo) : zero(x.lo)
           return Interval(xmin,xmax)
end

########### Defines sign
@inline sign(x::SMCg{N,T}) where {N,T} = -step(x)+2*step(x)
