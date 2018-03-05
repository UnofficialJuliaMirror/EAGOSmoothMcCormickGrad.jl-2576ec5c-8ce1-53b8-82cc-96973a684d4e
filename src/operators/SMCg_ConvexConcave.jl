@inline exp_cc(x,xL,xU) = line_seg(x,xL,exp(xL),xU,exp(xU)),dline_seg(x,xL,exp(xL),xU,exp(xU))
@inline exp_cv(x,xL,xU) = exp(x),exp(x)
function exp(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = exp_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = exp_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = exp_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = exp_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = exp_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = exp_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, exp(x.Intv),x.cnst, x.IntvBox, x.xref)
end

@inline exp2_cc(x,xL,xU) = line_seg(x,xL,exp2(xL),xU,exp2(xU)),dline_seg(x,xL,exp2(xL),xU,exp2(xU))
@inline exp2_cv(x,xL,xU) = exp2(x),exp2(x)*log(2)
function exp2(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = exp2_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = exp2_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = exp2_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = exp2_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = exp2_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = exp2_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, exp2(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline exp10_cc(x,xL,xU) = line_seg(x,xL,exp10(xL),xU,exp10(xU)),dline_seg(x,xL,exp10(xL),xU,exp10(xU))
@inline exp10_cv(x,xL,xU) = exp10(x),exp10(x)*log(10)
function exp10(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = exp10_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = exp10_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = exp10_cc(x.cv,x.Intv.lo,x.Intv.hi)
    gcv1,gdcv1 = exp10_cv(x.cv,x.Intv.lo,x.Intv.hi)
    gcc2,gdcc2 = exp10_cc(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = exp10_cv(x.cc,x.Intv.lo,x.Intv.hi)
    cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, exp10(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline sqrt_cc(x,xL,xU) = sqrt(x),1/(2*sqrt(x))
@inline sqrt_cv(x,xL,xU) = line_seg(x,xL,sqrt(xL),xU,sqrt(xU)),dline_seg(x,xL,sqrt(xL),xU,sqrt(xU))
function sqrt(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = sqrt_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = sqrt_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = sqrt_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = sqrt_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = sqrt_cc(x.cc,x.Intv.lo,x.Intv.hi)
    gcv2,gdcv2 = sqrt_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, sqrt(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline log_cc(x,xL,xU) = log(x),1/x
@inline log_cv(x,xL,xU) = line_seg(x,xL,log(xL),xU,log(xU)),dline_seg(x,xL,log(xL),xU,log(xU))
function log(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = log_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = log_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = log_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = log_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = log_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = log_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
    cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, log(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline log2_cc(x,xL,xU) = log2(x),1/(log(2)*x)
@inline log2_cv(x,xL,xU) = line_seg(x,xL,log2(xL),xU,log2(xU)),dline_seg(x,xL,log2(xL),xU,log2(xU))
function log2(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = log2_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = log2_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = log2_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = log2_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = log2_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = log2_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, log2(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline log10_cc(x,xL,xU) = log10(x),1/(log(10)*x)
@inline log10_cv(x,xL,xU) = line_seg(x,xL,log10(xL),xU,log10(xU)),dline_seg(x,xL,log10(xL),xU,log10(xU))
function log10(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = log10_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = log10_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = log10_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = log10_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = log10_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = log10_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, log10(x.Intv),x.cnst,x.IntvBox, x.xref)
end

@inline acosh_cc(x,xL,xU) = acosh(x),1/sqrt(x^2 - 1)
@inline acosh_cv(x,xL,xU) = line_seg(x,xL,acosh(xL),xU,acosh(xU)),dline_seg(x,xL,acosh(xL),xU,acosh(xU))
function acosh(x::SMCg{N,T}) where {N,T}
  eps_max = x.Intv.hi
  eps_min = x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = acosh_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = acosh_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = acosh_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = acosh_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = acosh_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = acosh_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, acosh(x.Intv),x.cnst,x.IntvBox, x.xref)
end

function abs_cv(x,xL,xU)
  if (x>=zero(x))
    return xU*(x/xU)^(MC_param.mu+1),(MC_param.mu+1)*(x/xU)^MC_param.mu
  else
    return -xL*(x/xL)^(MC_param.mu+1), -(MC_param.mu+1)*(x/xL)^MC_param.mu
  end
end
@inline abs_cc(x,xL,xU) = line_seg(x,xL,abs(xL),xU,abs(xU)),dline_seg(x,xL,abs(xL),xU,abs(xU))
function abs_cc_NS(x,lo,hi) # DONE (- Gradient)
  return line_seg(x,lo,abs(lo),hi,abs(hi)),dline_seg(x,lo,abs(lo),hi,abs(hi))
end
function abs_cv_NS(x,lo,hi) # DONE (- Gradient)
  return abs(x),sign(x)
end
function abs(x::SMCg{N,T}) where {N,T}
  #println("abs x in: ", x)
  eps_min,blank = mid3(x.Intv.lo,x.Intv.hi,zero(x.Intv.lo))
  eps_max = (abs(x.Intv.hi)>=abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
    cc,dcc = abs_cc(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = abs_cv(midcv,x.Intv.lo,x.Intv.hi)
    gcc1,gdcc1 = abs_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = abs_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = abs_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = abs_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc,dcc = abs_cc_NS(midcc,x.Intv.lo,x.Intv.hi)
    cv,dcv = abs_cv_NS(midcv,x.Intv.lo,x.Intv.hi)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, abs(x.Intv),x.cnst ,x.IntvBox, x.xref)
end

@inline cosh_cv(x,xL,xU) = cosh(x),-sinh(x)
@inline cosh_cc(x,xL,xU) = line_seg(x,xL,cosh(xL),xU,cosh(xU)),dline_seg(x,xL,cosh(xL),xU,cosh(xU))
function cosh(x::SMCg{N,T}) where {N,T}
  eps_min,blank = mid3(x.Intv.lo,x.Intv.hi,zero(x.Intv.lo))
  eps_max = (abs(x.Intv.hi)>=abs(x.Intv.lo)) ? x.Intv.hi : x.Intv.lo
  midcc,cc_id = mid3(x.cc,x.cv,eps_max)
  midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  cc,dcc = cosh_cc(midcc,x.Intv.lo,x.Intv.hi)
  cv,dcv = cosh_cv(midcv,x.Intv.lo,x.Intv.hi)
  if (MC_param.mu >= 1)
    gcc1,gdcc1 = cosh_cc(x.cv,x.Intv.lo,x.Intv.hi)
	  gcv1,gdcv1 = cosh_cv(x.cv,x.Intv.lo,x.Intv.hi)
	  gcc2,gdcc2 = cosh_cc(x.cc,x.Intv.lo,x.Intv.hi)
	  gcv2,gdcv2 = cosh_cv(x.cc,x.Intv.lo,x.Intv.hi)
	  cv_grad = max(0.0,gdcv1)*x.cv_grad_grad + min(0.0,gdcv2)*x.cc_grad
	  cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, cosh(x.Intv), x.cnst ,x.IntvBox, x.xref)
end
