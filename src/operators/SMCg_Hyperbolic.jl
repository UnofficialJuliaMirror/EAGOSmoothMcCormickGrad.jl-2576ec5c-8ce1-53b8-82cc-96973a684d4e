@inline function sinh_env(x,y,z)
  return sinh(y)-sinh(x)-(y-x)*cosh(x)
end
@inline function sinh_envd(x,y,z)
  return (y-x)*sinh(x)
end
@inline function cv_sinh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return sinh(x),cosh(x)
  elseif (xU<=zero(x))
    return line_seg(x,xL,sinh(xL),xU,sinh(xU)),dline_seg(x,xL,sinh(xL),xU,sinh(xU))
  else
    try
      p = newton(xL,xL,zero(x),sinh_env,sinh_envd,xU,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),sinh_env,xU,zero(x))
      end
    end
    if (x>p)
      return sinh(x),cosh(x)
    else
      return line_seg(x,p,sinh(p),xU,sinh(xU)),dline_seg(x,p,sinh(p),xU,sinh(xU))
    end
  end
end
@inline function cc_sinh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,sinh(xL),xU,sinh(xU)),dline_seg(x,xL,sinh(xL),xU,sinh(xU))
  elseif (xU<=zero(x))
    return sinh(x),cosh(x)
  else
    try
      p = newton(xU,zero(x),xU,sinh_env,sinh_envd,xL,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,sinh_env,xL,zero(x))
      end
    end
    if (x>p)
      return line_seg(x,xL,sinh(xL),p,sinh(p))
    else
      return sinh(x),cosh(x)
    end
  end
end
@inline function sinh(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_sinh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_sinh(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_sinh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_sinh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_sinh(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_sinh(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, sinh(x.Intv),x.cnst, x.IntvBox, x.xref)
end

@inline function asinh_env(x,y,z)
  val = (asinh(z)-asinh(x))/(z-x)-1.0/sqrt(1.0+x^2.0)
  return val
end
@inline function asinh_envd(x,y,z)
  val = (asinh(z)-asinh(x))/(z-x)^2.0+x/(x^2.0+1.0)^(1.5)-1.0/((z-x)*sqrt(x^2.0+1.0))
  return val
end
@inline function cv_asinh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,asinh(xL),xU,asinh(xU)),dline_seg(x,xL,asinh(xL),xU,asinh(xU))
  elseif (xU<=zero(x))
    return asinh(x),1/sqrt(x^2+1)
  else
    try
      p = newton(xL/2.0,xL,zero(x),asinh_env,asinh_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),asinh_env,xL,xU)
      end
    end
    if (x<=p)
      return asinh(x),1/sqrt(x^2+1)
    else
      return line_seg(x,p,asinh(p),xU,asinh(xU)),dline_seg(x,p,asinh(p),xU,asinh(xU))
    end
  end
end
@inline function cc_asinh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return asinh(x),1/sqrt(x^2+1)
  elseif (xU<=zero(x))
    return line_seg(x,xL,asinh(xL),xU,asinh(xU)),dline_seg(x,xL,asinh(xL),xU,asinh(xU))
  else
    try
      p = newton(xU/2,zero(x),xU,asinh_env,asinh_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,asinh_env,xL,xU)
      end
    end
    if (x<=p)
      return line_seg(x,xL,asinh(xL),p,asinh(p)),dline_seg(x,xL,asinh(xL),p,asinh(p))
    else
      return asinh(x),1/sqrt(x^2+1)
    end
  end
end
@inline function asinh(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_asinh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_asinh(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_asinh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_asinh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_asinh(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_asinh(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, asinh(x.Intv),x.cnst, x.IntvBox, x.xref)
end

@inline function tanh_env(x,y,z)
  return (tanh(y)-tanh(x))/(1.0-tanh(x)^2)-y+x
end
@inline function tanh_envd(x,y,z)
  return -((1.0-tanh(x)^2)^(-2.0))*(tanh(y)-tanh(x))
end
@inline function cv_tanh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,tanh(xL),xU,tanh(xU)),dline_seg(x,xL,tanh(xL),xU,tanh(xU))
  elseif (xU<=zero(x))
    return tanh(x),sech(x)^2
  else
    try
      p = newton(xL,xL,zero(x),tanh_env,tanh_envd,xU,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),tanh_env,xU,zero(x))
      end
    end
    if (x<=p)
      return tanh(x),sech(x)^2
    else
      return line_seg(x,p,tanh(p),xU,tanh(xU)),dline_seg(x,p,tanh(p),xU,tanh(xU))
    end
  end
end
@inline function cc_tanh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return tanh(x),sech(x)^2
  elseif (xU<=zero(x))
    return line_seg(x,xL,tanh(xL),xU,tanh(xU)),dline_seg(x,xL,tanh(xL),xU,tanh(xU))
  else
    try
      p = newton(xU/2,zero(x),xU,tanh_env,tanh_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,tanh_env,xL,zero(x))
      end
    end
    if (x<=p)
      return line_seg(x,xL,tanh(xL),p,tanh(p)),dline_seg(x,xL,tanh(xL),p,tanh(p))
    else
      return tanh(x),sech(x)^2
    end
  end
end
@inline function tanh(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_tanh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_tanh(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_tanh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_tanh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_tanh(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_tanh(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, tanh(x.Intv),x.cnst, x.IntvBox, x.xref)
end

@inline function atanh_env(x,y,z)
  return (x-y)-(1.0-x^2.0)*(atan(x)-atan(y))
end
@inline function atanh_envd(x,y,z)
  return 1+2.0*x*(atan(x)-atan(y))
end
@inline function cv_atanh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return atanh(x),1/(1-x^2)
  elseif (xU<=zero(x))
    return line_seg(x,xL,atanh(xL),xU,atanh(xU)),dline_seg(x,xL,atanh(xL),xU,atanh(xU))
  else
    try
      p = newton(xU,zero(x),xU,atanh_env,atanh_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),atanh_env,xU,zero(x))
      end
    end
    if (x>p)
      return atanh(x),1/(1-x^2)
    else
      return line_seg(x,p,atanh(p),xU,atanh(xU)),dline_seg(x,p,atanh(p),xU,atanh(xU))
    end
  end
end
@inline function cc_atanh(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,atanh(xL),xU,atanh(xU)),dline_seg(x,xL,atanh(xL),xU,atanh(xU))
  elseif (xU<=zero(x))
    return atanh(x),1/(1-x^2)
  else
    try
      p = newton(xL,xL,zero(x),atanh_env,atanh_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,atanh_env,xL,zero(x))
      end
    end
    if (x>p)
      return line_seg(x,xU,atanh(xU),p,atanh(p)),dline_seg(x,xU,atanh(xU),p,atanh(p))
    else
      return atanh(x),1/(1-x^2)
    end
  end
end
@inline function atanh(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_atanh(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_atanh(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_atanh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_atanh(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_atanh(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_atanh(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, atanh(x.Intv),x.cnst, x.IntvBox, x.xref)
end
