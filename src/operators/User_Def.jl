export FDI_invApprox

function FDI_invApprox_cc(y, yL, yU)
 Ap = 0.00005
 g = 9.8

 xL = -(0.5/g)*(yL/Ap)^2
 xU = (0.5/g)*(yU/Ap)^2

 val = (xU-xL)*sqrt((y-xL)/(yU-yL))+xL
 dval = 0.5*((xU-xL)/(yU-yL))/sqrt((y-yL)/(yU-yL))
 return val,dval
end

function FDI_invApprox_cv(y, yL, yU)
  Ap = 0.00005
  g = 9.8

  xL = -(0.5/g)*(yL/Ap)^2
  xU = (0.5/g)*(yU/Ap)^2

  val = xL + (xU - xL)*((y - yL)/(yU - yL))^2
  dval = 2.0*(xU - xL)*((y - yL)/(yU - yL))/(yU - yL)
  return val,dval
end

function FDI_invApprox(x::SMCg{N,T}) where {N,T}

  # parameters
  Ap = 0.00005
  g = 9.8

  # calculating relaxation on either side of interval
  upper_calc = ((x/Ap)^2)/(2.0*g)
  lower_calc = (-(x/Ap)^2)/(2.0*g)

  if (x.Intv.lo>0.0) # return upper relaxation if x above zero
    return upper_calc
  elseif (x.Intv.hi<0.0) # return lower relaxation if x below zero
    return lower_calc
  else
    eps_min = 0.0
    eps_max = ((x.Intv.lo/Ap)^2)/(2.0*g)>(-(x.Intv.hi/Ap)^2)/(2.0*g) ? x.Intv.lo : x.Intv.hi
    cc_mid3,cc_id = mid3(x.cc,x.cv,eps_min)
    cv_mid3,cv_id = mid3(x.cc,x.cv,eps_max)
    cc,dcc = FDI_invApprox_cc(cc_mid3, x.Intv.lo, x.Intv.hi)
    cv,dcv = FDI_invApprox_cv(cv_mid3, x.Intv.lo, x.Intv.hi)
    if (MC_param.mu >= 1)
      gcc1,gdcc1 = FDI_invApprox_cc(x.cv,x.Intv.lo,x.Intv.hi)
      gcv1,gdcv1 = FDI_invApprox_cv(x.cv,x.Intv.lo,x.Intv.hi)
      gcc2,gdcc2 = FDI_invApprox_cc(x.cc,x.Intv.lo,x.Intv.hi)
      gcv2,gdcv2 = FDI_invApprox_cv(x.cc,x.Intv.lo,x.Intv.hi)
      cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
      cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
    else
      cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
      cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
    end
    Intv = Interval(lower_calc.Intv.lo,upper_calc.Intv.hi)
    return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,x.cnst,x.IntvBox,x.xref)
  end
end

# inverse of functional term
function FDI_invApprox(x)
  Ap = 0.00005
  g = 9.8
  if x>=0.0
    return ((x/Ap)^2)/(2.0*g)
  else
    return (-(x/Ap)^2)/(2.0*g)
  end
end

# inverse of functional term
function FDI_invApprox(x::Interval)
  Ap = 0.00005
  g = 9.8
  upper = ((x/Ap)^2)/(2.0*g)
  lower = (-(x/Ap)^2)/(2.0*g)
  if x.lo>=0.0
    return upper
  elseif x.hi<=0
    return lower
  else
    return Interval(lower.lo,upper.hi)
  end
end
