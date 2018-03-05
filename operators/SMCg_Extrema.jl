@inline function cv_max(x,xL,xU,a)
        if (xU<=a)
          return a, 0
        elseif (a<=xL)
          return x, 1
        else
          val = a + (xU-a)*(max(zero(a),((x-a)/(xU-a))))^(MC_param.mu+1)
          dval = ((xU-a)/(xU-a))*(MC_param.mu+1)*(max(zero(a),((x-a)/(xU-a))))^(MC_param.mu)
          return val, dval
        end
end
@inline function cc_max(x,xL,xU,a)
        return line_seg(x,xL,max(xL,a),xU,max(xU,a)),dline_seg(x,xL,max(xL,a),xU,max(xU,a))
end
function cc_max_NS(x,lo,hi,c)
  return line_seg(x,lo,max(lo,c),hi,max(hi,c)),dline_seg(x,lo,max(lo,c),hi,max(hi,c))
end
function cv_max_NS(x,lo,hi,c)
  return max(x,c),max(1,0)
end

for i in union(int_list, float_list)
	eval( quote
	char = $i
@inline function max(x::SMCg{N,T},c::$i) where {N,T}
  eps_min = x.Intv.lo
  eps_max = x.Intv.hi
	midcc,cc_id = mid3(x.cc,x.cv,eps_max)
	midcv,cv_id = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
     #println("ran me max 1")
	   cc,dcc = cc_max(midcc,x.Intv.lo,x.Intv.hi,c)
     cv,dcv = cv_max(midcv,x.Intv.lo,x.Intv.hi,c)
     gcc1,gdcc1 = cc_max(x.cv,x.Intv.lo,x.Intv.hi)
 	   gcv1,gdcv1 = cv_max(x.cv,x.Intv.lo,x.Intv.hi)
 	   gcc2,gdcc2 = cc_max(x.cc,x.Intv.lo,x.Intv.hi)
 	   gcv2,gdcv2 = cv_max(x.cc,x.Intv.lo,x.Intv.hi)
 	   cv_grad = max(0.0,gdcv1)*x.cv_grad + min(0.0,gdcv2)*x.cc_grad
 	   cc_grad = min(0.0,gdcc1)*x.cv_grad + max(0.0,gdcc2)*x.cc_grad
  else
    #println("ran me max 2")
    cc,dcc = cc_max_NS(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_max_NS(midcv,x.Intv.lo,x.Intv.hi,c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end

  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, max(x.Intv),x.cnst,x.IntvBox,x.xref)
end
@inline max(c::$i,x::SMCg{N,T}) where {N,T} = max(x,c)
@inline min(c::$i,x::SMCg{N,T}) where {N,T} = -max(-x,-c)
@inline min(x::SMCg{N,T},c::$i) where {N,T} = -max(-x,-c)
         end )
 end
# defines functions on which bivariant maximum mapping from Khan 2016
@inline function psil_max(x,y,lambda::Interval,nu::Interval,f1,f2)
   if (nu.hi<=lambda.lo)
     val = x
   elseif (lambda.hi<=nu.lo)
     val = y
   elseif ((nu.lo<=lambda.lo)&&(lambda.lo<nu.hi))
     val = x+(nu.hi-lambda.lo)*max(0.0,((y-x)/(nu.hi-lambda.lo)))^(MC_param.mu+1)
   else
     val =  y + (lambda.hi-nu.lo)*max(0.0,(x-y)/(lambda.hi-nu.lo))^(MC_param.mu+1)
   end
   if (nu.hi <= lambda.lo)
     grad_val = f1.cv_grad
   elseif (lambda.hi <= nu.lo)
     grad_val = f1.cc_grad
   else
     grad_val = max(0,psil_max_dx(x,y,lambda,nu))*f1.cv_grad +
                min(0,psil_max_dx(x,y,lambda,nu))*f1.cc_grad +
                max(0,psil_max_dy(x,y,lambda,nu))*f2.cv_grad +
                min(0,psil_max_dy(x,y,lambda,nu))*f2.cc_grad
   end
   return val,grad_val
end
@inline function thetar(x,y,lambda::Interval,nu::Interval)
    return (max(lambda.lo,nu.lo) + max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi) +
    max(lambda.hi,nu.lo))*max(0.0,((lambda.hi-x)/(lambda.hi-lambda.lo)-(y-nu.lo)/(nu.hi-nu.lo)))^(MC_param.mu+1)
end
function psil_max_dx(x,y,lambda,nu)
  if (nu.lo <= lambda.lo < nu.hi)
    return 1-(MC_param.mu+1)*max(0,(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return (MC_param.mu+1)*max(0,(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end
function psil_max_dy(x,y,lambda,nu)
  if (nu.lo <= lambda.lo < nu.hi)
    return (MC_param.mu+1)*max(0,(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return 1-(MC_param.mu+1)*max(0,(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end

@inline function psir_max(x,y,xgrad,ygrad,lambda::Interval,nu::Interval)
    if (nu.hi<=lambda.lo)
      return x,xgrad
    elseif (lambda.hi<=nu.lo)
      return y,ygrad
    else
      val = max(lambda.hi,nu.hi)-(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))*
          ((lambda.hi-x)/(lambda.hi-lambda.lo))-
          (max(lambda.hi,nu.hi)-max(lambda.hi,nu.lo))*((nu.hi-y)/(nu.hi-nu.lo))
          +thetar(x,y,nu,lambda)
      coeff = [(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))/(lambda.hi-lambda.lo)
               (max(lambda.hi,nu.hi)-max(nu.lo,lambda.hi))/(nu.hi-nu.lo)
               (MC_param.mu+1)*(max(lambda.hi,nu.lo)+max(lambda.lo,nu.hi)-max(lambda.lo,nu.lo)-max(lambda.hi,nu.hi))
               max(0,((lambda.hi-x)/(lambda.hi-lambda.lo))-((y-nu.lo)/(nu.hi-nu.lo)))^MC_param.mu
               1/(lambda.hi-lambda.lo)
               1/(nu.hi-nu.lo)
               ]
      grad_val = coeff[1]*xgrad + coeff[2]*ygrad +
                 coeff[3]*coeff[4]*(coeff[5]*xgrad+coeff[6]*ygrad)
      return val,grad_val
    end
end

@inline function max(x::SMCg{N,T},y::SMCg{N,T}) where {N,T}
    if (MC_param.mu >= 1)
      cc = zero(x.cc)
      cv = zero(x.cc)
      temp_mid = zero(x.cc)
      if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
        cv,cv_grad = psil_max(x.cv,y.cv,x.Intv,y.Intv,x,y)
      elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
        temp_mid,blank = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
        cv,cv_grad = psil_max(temp_mid,y.cv,x.Intv,y.Intv,x,y)
      elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
        temp_mid,blank = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
        cv,cv_grad = psil_max(x.cv,temp_mid,x.Intv,y.Intv,x,y)
      end
      cc,cc_grad = psir_max(x.cc,y.cc,x.cc_grad,y.cv_grad,x.Intv,y.Intv)
      return SMCg{T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),x.IntvBox,x.xref)
    elseif (x.Intv.hi <= y.Intv.lo)
      cc = y.cc
      cc_grad = y.cnst ? zeros(y.cc_grad) : y.cc_grad
    elseif (x.Intv.lo >= y.Intv.hi)
      cc = x.cc
      cc_grad = x.cnst ? zeros(x.cc_grad) : x.cc_grad
    elseif (MC_param.multivar_refine)
      maxLL = max(x.Intv.lo,y.Intv.lo)
      maxLU = max(x.Intv.lo,y.Intv.hi)
      maxUL = max(x.Intv.hi,y.Intv.lo)
      maxUU = max(x.Intv.hi,y.Intv.hi)
      thin1 = (diam(x.Intv) == 0.0)
      thin2 = (diam(y.Intv) == 0.0)
      r11 = thin1 ? 0.0 : (maxUL-maxLL)/diam(x.Intv)
      r21 = thin1 ? 0.0 : (maxLU-maxUU)/diam(x.Intv)
      r12 = thin2 ? 0.0 : (maxLU-maxLL)/diam(y.Intv)
      r22 = thin2 ? 0.0 : (maxUL-maxUU)/diam(y.Intv)
      cc1 = maxLL + r11*(x.cc-x.Intv.lo) + r12*(y.cc-y.Intv.lo)
      cc2 = maxUU - r21*(x.cc-x.Intv.hi) - r22*(y.cc-y.Intv.hi)
      if (cc1 <= cc2)
        cc = cc1
        cc_grad =  (x.cnst ? zeros(y.cc_grad) : r11*x.cc_grad) + (y.cnst ? zeros(x.cc_grad) : r12*y.cc_grad)
      else
        cc = cc2
        cc_grad = -(x.cnst ? zeros(y.cc_grad) : r21*x.cc_grad) - (y.cnst ? zeros(x.cc_grad) : r22*y.cc_grad)
      end
    else
      ccMC = 0.5*(x+y+abs(x-y))
      cc = ccMC.cc
      cc_grad = ccMC.cc_grad
    end
    cv = max(x.cv,y.cv)
    cv_grad = (x.cv > y.cv) ? (x.cnst ? zeros(x.cv_grad): x.cv_grad) :
                              (y.cnst ? zeros(y.cv_grad): y.cv_grad)
    cnst = y.cnst ? x.cnst : (x.cnst ? y.cnst : (x.cnst || y.cnst) )

    return SMCg{N,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),cnst,x.IntvBox,x.xref)
end

@inline function maxcv(x::SMCg{N,T},y::SMCg{N,T}) where {N,T}
        cv = zero(x.cc)
        temp_mid = zero(x.cc)
        if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
            cv = psil_max(x.cv,y.cv,x.Intv,y.Intv)
        elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
          temp_mid = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
          cv = psil_max(temp_mid,y.cv,x.Intv,y.Intv)
        elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
          temp_mid = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1.0/MC_param.mu))
          cv = psil_max(x.cv,temp_mid,x.Intv,y.Intv)
        end
        return cv
end
@inline function maxcc(x,y::SMCg{N,T}) where {N,T}
        cc = psir_max(x.cc,y.cc,x.Intv,y.Intv)
end
@inline mincv(x,y::SMCg{N,T}) where {N,T} = - maxcc(-x,-y)
@inline min(x::SMCg{N,T},y::SMCg{N,T}) where {N,T} = -max(-x,-y)
@inline max(x::SMCg{N,T},y::Interval{T}) where {N,T} = max(x,SMCg{N,T}(y))
@inline max(y::Interval{T},x::SMCg{N,T}) where {N,T} = max(x,SMCg{N,T}(y))
@inline min(x::SMCg{N,T},y::Interval{T}) where {N,T} = min(x,SMCg{N,T}(y))
@inline min(y::Interval{T},x::SMCg{N,T}) where {N,T} = min(x,SMCg{N,T}(y))
