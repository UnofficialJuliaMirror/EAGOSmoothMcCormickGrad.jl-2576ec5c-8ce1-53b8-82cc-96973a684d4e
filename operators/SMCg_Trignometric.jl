# convex relaxation (envelope) of cos function
function cv_cos(x,xL,xU)
  r = zero(x)
  kL = Base.ceil((-0.5-xL/(2.0*pi)))
  if (x<=(pi-2.0*pi*kL))
    xL1 = xL+2.0*pi*kL
    if (xL1 >= pi/2.0)
      return cos(x),-sin(x)
    end
    xU1 = min(xU+2.0*pi*kL,pi)
    if ((xL1>=(-pi/2))&&(xU1<=(pi/2)))
      if (abs(xL-xU)<MC_param.env_tol)
        r = zero(x)
      else
        r = (cos(xU)-cos(xL))/(xU-xL)
      end
      return cos(xL)+r*(x-xL),r
    end
    return cv_cosin(x+(2.0*pi)*kL,xL1,xU1)
  end
  kU = Base.floor((0.5-xU/(2.0*pi)))
  if (x>=(-pi-2.0*pi*kU))
    xU2 = xU+2.0*pi*kU
    if (xU2<=-pi/2.0)
      return cos(x),-sin(x)
    end
    return cv_cosin(x+2.0*pi*kU,max(xL+2.0*pi*kU,-pi),xU2)
  end
  return -1.0,0.0
end
# function for computing convex relaxation over nonconvex and nonconcave regions
function cv_cosin(x,xL,xU)
  xj = -Inf
  if (abs(xL)<=abs(xU))
    left::Bool = false
    x0 = xU
    xm = xL
  else
    left = true
    x0 = xL
    xm = xU
  end
  try
    xj = newton(x0,xL,xU,cv_cosenv,dcv_cosenv,xm,zero(x0))
  catch e
    if isa(e, ErrorException)
      xj = golden_section(xL,xU,cv_cosenv,xm,zero(x0))
    end
  end
  if ((left && x<=xj)||((~left) && x>=xj))
    return cos(x),-sin(x)
  else
    if abs(xm-xj)<MC_param.env_tol
      r = zero(x0)
    else
      r = (cos(xm)-cos(xj))/(xm-xj)
    end
    return cos(xm)+r*(x-xm),r
  end
end
# pivot point calculation function for convex relaxation of cosine
function cv_cosenv(x,y,z)
  return (x-y)*sin(x)+cos(x)-cos(y)
end
function dcv_cosenv(x,y,z)
  return (x-y)*cos(x)
end
# concave relaxation (envelope) of cos function
function cc_cos(x,xL,xU)
  temp = cv_cos(x-pi,xL-pi,xU-pi)
  return -temp[1],-temp[2]
end
function cos_arg(xL,xU)
  kL = Base.ceil(-0.5-xL/(2.0*pi))
  xL1 = xL+2.0*pi*kL
  xU1 = xU+2.0*pi*kL
  if ~((xL1>=-pi)&&(xL1<=pi))
    error("Cosine Argument Calculation: xL out of bounds.")
  end
  if (xL1<=zero(xL))
    if (xU1<=zero(xL))
      arg1 = xL
      arg2 = xU
    elseif (xU1>=pi)
      arg1 = pi-2.0*pi*kL
      arg2 = -2.0*pi*kL
    else
      arg1 = (cos(xL1)<=cos(xU1)) ? xL : xU
      arg2 = -2.0*pi*kL
    end
  end
  if (xU1<=pi)
    arg1 = xU
    arg2 = xL
  elseif (xU1>=(2.0*pi))
    arg1 = pi-2.0*pi*kL
    arg2 = 2.0*pi-2.0*pi*kL
  else
    arg1 = pi-2.0*pi*kL
    arg2 = (cos(xL1)>=cos(xU1)) ? xL : xU
  end
  return [arg1,arg2]
end

function cos(x::SMCg{N,T}) where {N,T}
  #println("x for cos: ", x)
  eps_max,eps_min = cos_arg(x.Intv.lo,x.Intv.hi)
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
	midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_cos(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_cos(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu>=1)
    gcc1,gdcc1 = cc_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = cv_cos(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = cc_cos(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = cv_cos(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, cos(x.Intv),x.cnst, x.IntvBox, x.xref)
end

function sin(x::SMCg{N,T}) where {N,T}
  cos(x-pi/2.0)
end

# pivot point calculation function for convex relaxation of complimentary error function
@inline function tan_env(x,y,z)
  return (x-y)-(tan(x)-tan(y))/(1.0+tan(x)^2.0)
end
# derivative of pivot point calculation function for convex relaxation of complimentary error function
@inline function tan_envd(x,y,z)
  return 2.0*tan(x)/(1.0+tan(x)^2.0)*(tan(x)-tan(y))
end
# convex relaxation (envelope) of tangent function
@inline function cv_tan(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return tan(x),sec(x)^2
  elseif (xU<=zero(x))
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU))
  else
    try
      p = secant(zero(x),xU,zero(x),xU,tan_env,xL,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,tan_env,xL,zero(x))
      end
    end
    if (x<=p)
      return line_seg(x,xL,tan(xL),p,tan(p)),dline_seg(x,xL,tan(xL),p,tan(p))
    else
      return tan(x),sec(x)^2
    end
  end
end
# concave relaxation (envelope) of tangent function
@inline function cc_tan(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,tan(xL),xU,tan(xU)),dline_seg(x,xL,tan(xL),xU,tan(xU))
  elseif (xU<=zero(x))
    return tan(x),sec(x)^2
  else
    try
      p = secant(zero(x),xL,xL,zero(x),tan_env,xU,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),tan_env,xU,zero(x))
      end
    end
    if (x<=p)
       return tan(x),sec(x)^2
    else
       return line_seg(x,p,tan(p),xU,tan(xU)),dline_seg(x,p,tan(p),xU,tan(xU))
     end
  end
end
@inline function tan(x::SMCg{N,T}) where {N,T}
  if ((x.Intv.lo==-Inf)||(x.Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_tan(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_tan(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_tan(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_tan(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_tan(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_tan(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, tan(x.Intv),x.cnst, x.IntvBox, x.xref)
end

@inline function acos(x::SMCg{N,T}) where {N,T}
  return asin(-x)+pi/2.0
end

# pivot point calculation function for convex relaxation of arcsine
@inline function asin_env(x,y,z)
  return (asin(x)-asin(y))/(x-y)-1.0/sqrt(1.0-x^2.0)
end
# derivative of pivot point calculation function for convex relaxation of arcsine
@inline function asin_envd(x,y,z)
  return 1.0/((x-y)*sqrt(1.0-x^2.0))-x/(1.0-x^2.0)^(1.5)-(asin(x)-asin(y))/(x-y)^2.0
end
# convex relaxation (envelope) of arcsine function
@inline function cv_asin(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return asin(x),1/sqrt(1-x^2)
  elseif (xU<=zero(x))
    return line_seg(x,xL,asin(xL),xU,asin(xU)),dline_seg(x,xL,asin(xL),xU,asin(xU))
  else
    try
      p = newton(xU,xL,xU,asin_env,asin_envd,xL,xU)
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,asin_env,xL,zero(x))
      end
    end
    if (x<=p)
      return line_seg(x,xL,asin(xL),p,asin(p)),dline_seg(x,xL,asin(xL),p,asin(p))
    else
      return asin(x),1/sqrt(1-x^2)
    end
  end
end
# concave relaxation (envelope) of arcsine function
@inline function cc_asin(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,asin(xL),xU,asin(xU)),dline_seg(x,xL,asin(xL),xU,asin(xU))
  elseif (xU<=zero(x))
    return asin(x),1/sqrt(1-x^2)
  else
    try
      p = secant(zero(x),xL,xL,zero(x),asin_env,xU,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),asin_env,xU,zero(x))
      end
    end
    if (x<=p)
      return asin(x),1/sqrt(1-x^2)
    else
      return line_seg(x,p,asin(p),xU,asin(xU)),dline_seg(x,p,asin(p),xU,asin(xU))
    end
  end
end
@inline function asin(x::SMCg{N,T}) where {N,T}
  if ((x.Intv.lo==-Inf)||(x.Intv.hi==Inf))
    error("Function unbounded on domain")
  end
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_asin(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_asin(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_asin(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_asin(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_asin(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_asin(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, asin(x.Intv),x.cnst, x.IntvBox, x.xref)
end

# pivot point calculation function for convex relaxation of arctangent
@inline function atan_env(x,y,z)
  return (x-y)-(1.0+x^2.0)*(atan(x)-atan(y))
end
# derivative of pivot point calculation function for convex relaxation of arctangent
@inline function atan_envd(x,y,z)
  return -2.0*x*(atan(x)-atan(y))
end
# convex relaxation (envelope) of arctan function
@inline function cv_atan(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return line_seg(x,xL,atan(xL),xU,atan(xU)),dline_seg(x,xL,atan(xL),xU,atan(xU))
  elseif (xU<=zero(x))
    return atan(x),1/(x^2+1)
  else
    try
      p = newton(xL,xL,zero(x),atan_env,atan_envd,xU,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(xL,zero(x),atan_env,xU,zero(x))
      end
    end
    if (x<=p)
      return atan(x),1/(x^2+1)
    else
      return line_seg(x,p,atan(p),xU,atan(xU)),dline_seg(x,p,atan(p),xU,atan(xU))
    end
  end
end
# concave relaxation (envelope) of arctan function
@inline function cc_atan(x,xL,xU)
  p = zero(x)
  if (xL>=zero(x))
    return atan(x),1/(x^2+1)
  elseif (xU<=zero(x))
    return line_seg(x,xL,atan(xL),xU,atan(xU)),dline_seg(x,xL,atan(xL),xU,atan(xU))
  else
    try
      p = newton(xU,zero(x),xU,atan_env,atan_envd,xL,zero(x))
    catch e
      if isa(e, ErrorException)
        p = golden_section(zero(x),xU,atan_env,xL,zero(x))
      end
    end
    if (x<=p)
      return line_seg(x,xL,atan(xL),p,atan(p)),dline_seg(x,xL,atan(xL),p,atan(p))
    else
      return atan(x),1/(x^2+1)
    end
  end
end
@inline function atan(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cc_atan(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cv_atan(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cc_atan(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcv1,gdcv1 = cv_atan(x.cv,x.Intv.lo,x.Intv.hi,c)
    gcc2,gdcc2 = cc_atan(x.cc,x.Intv.lo,x.Intv.hi,c)
    gcv2,gdcv2 = cv_atan(x.cc,x.Intv.lo,x.Intv.hi,c)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, atan(x.Intv),x.cnst, x.IntvBox, x.xref)
end
