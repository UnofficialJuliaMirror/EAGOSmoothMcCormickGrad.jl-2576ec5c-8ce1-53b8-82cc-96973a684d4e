@inline function cv_max(x::T,xL::T,xU::T,a::S) where {S,T<:AbstractFloat}
        ca::T = convert(T,a)
        if (xU<=ca)
          return ca, zero(T)
        elseif (ca<=xL)
          return x, one(T)
        else
          val::T = ca + (xU-ca)*(max(0.0,((x-ca)/(xU-ca))))^(MC_param.mu+1)
          dval::T = ((xU-ca)/(xU-ca))*(MC_param.mu+1)*(max(0.0,((x-ca)/(xU-ca))))^(MC_param.mu)
          return val, dval
        end
end
@inline function cc_max(x::T,xL::T,xU::T,a::S) where {S,T<:AbstractFloat}
  ca::T = convert(T,a)
  if (xU<=ca)
    d = zero(T)
  elseif (ca<=xL)
    d = one(T)
  else
    d = ((xU-ca)/(xU-ca))*(MC_param.mu+1)*(max(0.0,((x-ca)/(xU-ca))))^(MC_param.mu)
  end
  return line_seg(x,xL,max(xL,a),xU,max(xU,a)),dline_seg(x,xL,max(xL,a),xU,max(xU,a),d)
end
function cc_max_NS(x::T,lo::T,hi::T,c::S) where {S,T<:AbstractFloat}
  return line_seg(x,lo,max(lo,c),hi,max(hi,c)),dline_seg(x,lo,max(lo,c),hi,max(hi,c),((x>c) ? 1.0 : 0.0))
end
function cv_max_NS(x::T,xL::T,xU::T,c::S) where {S,T<:AbstractFloat}
  return max(x,c),((x>c) ? 1.0 : 0.0)
end

for i in union(int_list, float_list)
	eval( quote
	char = $i
@inline function max(x::SMCg{N,T},c::$i) where {N,T<:AbstractFloat}
  eps_min::T = x.Intv.lo
  eps_max::T = x.Intv.hi
	midcc::T,cc_id::Int64 = mid3(x.cc,x.cv,eps_max)
	midcv::T,cv_id::Int64 = mid3(x.cc,x.cv,eps_min)
  if (MC_param.mu >= 1)
     #println("ran me max 1")
	   cc::T,dcc::T = cc_max(midcc,x.Intv.lo,x.Intv.hi,c)
     cv::T,dcv::T = cv_max(midcv,x.Intv.lo,x.Intv.hi,c)
     gcc1::T,gdcc1::T = cc_max(x.cv,x.Intv.lo,x.Intv.hi,c)
 	   gcv1::T,gdcv1::T = cv_max(x.cv,x.Intv.lo,x.Intv.hi,c)
 	   gcc2::T,gdcc2::T = cc_max(x.cc,x.Intv.lo,x.Intv.hi,c)
 	   gcv2::T,gdcv2::T = cv_max(x.cc,x.Intv.lo,x.Intv.hi,c)
 	   cv_grad::SVector{N,T} = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
 	   cc_grad::SVector{N,T} = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
  else
    #println("ran me max 2")
    cc,dcc = cc_max_NS(midcc,x.Intv.lo,x.Intv.hi,c)
    cv,dcv = cv_max_NS(midcv,x.Intv.lo,x.Intv.hi,c)
    cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
    cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
  end

  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, max(x.Intv),x.cnst,x.IntvBox,x.xref)
end
@inline max(c::$i,x::SMCg{N,T}) where {N,T<:AbstractFloat} = max(x,c)
@inline min(c::$i,x::SMCg{N,T}) where {N,T<:AbstractFloat} = -max(-x,-c)
@inline min(x::SMCg{N,T},c::$i) where {N,T<:AbstractFloat} = -max(-x,-c)
         end )
 end
# defines functions on which bivariant maximum mapping from Khan 2016
@inline function psil_max(x::T,y::T,lambda::Interval{T},nu::Interval{T},
                          f1::SMCg{N,T},f2::SMCg{N,T}) where {N,T<:AbstractFloat}
   if (nu.hi<=lambda.lo)
     val::T = x
   elseif (lambda.hi<=nu.lo)
     val = y
   elseif ((nu.lo<=lambda.lo)&&(lambda.lo<nu.hi))
     val = x+(nu.hi-lambda.lo)*max(zero(T),((y-x)/(nu.hi-lambda.lo)))^(MC_param.mu+1)
   else
     val =  y + (lambda.hi-nu.lo)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^(MC_param.mu+1)
   end
   if (nu.hi <= lambda.lo)
     grad_val::SVector{N,T} = f1.cv_grad
   elseif (lambda.hi <= nu.lo)
     grad_val = f1.cc_grad
   else
     grad_val = max(zero(T),psil_max_dx(x,y,lambda,nu))*f1.cv_grad +
                min(zero(T),psil_max_dx(x,y,lambda,nu))*f1.cc_grad +
                max(zero(T),psil_max_dy(x,y,lambda,nu))*f2.cv_grad +
                min(zero(T),psil_max_dy(x,y,lambda,nu))*f2.cc_grad
   end
   return val,grad_val
end
@inline function thetar(x::T,y::T,lambda::Interval{T},nu::Interval{T}) where {T<:AbstractFloat}
    return (max(lambda.lo,nu.lo) + max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi) +
    max(lambda.hi,nu.lo))*max(zero(T),((lambda.hi-x)/(lambda.hi-lambda.lo)-(y-nu.lo)/(nu.hi-nu.lo)))^(MC_param.mu+1)
end
function psil_max_dx(x::T,y::T,lambda::Interval{T},nu::Interval{T}) where {T<:AbstractFloat}
  if (nu.lo <= lambda.lo < nu.hi)
    return one(T)-(MC_param.mu+1)*max(zero(T),(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return (MC_param.mu+1)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end
function psil_max_dy(x::T,y::T,lambda::Interval{T},nu::Interval{T}) where {T<:AbstractFloat}
  if (nu.lo <= lambda.lo < nu.hi)
    return (MC_param.mu+1)*max(zero(T),(y-x)/(nu.hi-lambda.lo))^MC_param.mu
  else
    return one(T)-(MC_param.mu+1)*max(zero(T),(x-y)/(lambda.hi-nu.lo))^MC_param.mu
  end
end

@inline function psir_max(x::T,y::T,xgrad::SVector{N,T},ygrad::SVector{N,T},
                          lambda::Interval{T},nu::Interval{T}) where {N,T<:AbstractFloat}
    if (nu.hi<=lambda.lo)
      return x,xgrad
    elseif (lambda.hi<=nu.lo)
      return y,ygrad
    else
      val::T = max(lambda.hi,nu.hi)-(max(lambda.hi,nu.hi)-max(lambda.lo,nu.hi))*
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
      grad_val::SVector{N,T} = coeff[1]*xgrad + coeff[2]*ygrad +
                 coeff[3]*coeff[4]*(coeff[5]*xgrad+coeff[6]*ygrad)
      return val,grad_val
    end
end

@inline function max(x::SMCg{N,T},y::SMCg{N,T}) where {N,T<:AbstractFloat}
    if (MC_param.mu >= 1)
      cc::T = zero(T)
      cv::T = zero(T)
      temp_mid::T = zero(T)
      if ((y.Intv.hi<=x.Intv.lo)||(x.Intv.hi<=y.Intv.lo))
        cv,cv_grad::SVector{N,T} = psil_max(x.cv,y.cv,x.Intv,y.Intv,x,y)
      elseif ((y.Intv.lo<=x.Intv.lo) & (x.Intv.lo<y.Intv.hi))
        temp_mid,blank = mid3(x.cv,x.cc,y.cv-(y.Intv.hi-x.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(temp_mid,y.cv,x.Intv,y.Intv,x,y)
      elseif ((x.Intv.lo<y.Intv.lo) & (y.Intv.lo<x.Intv.hi))
        temp_mid,blank = mid3(y.cv,y.cc,x.cv-(x.Intv.hi-y.Intv.lo)*(MC_param.mu+1)^(-1/MC_param.mu))
        cv,cv_grad = psil_max(x.cv,temp_mid,x.Intv,y.Intv,x,y)
      end
      cc,cc_grad::SVector{N,T} = psir_max(x.cc,y.cc,x.cc_grad,y.cv_grad,x.Intv,y.Intv)
      return SMCg{N,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),(x.cnst && y.cnst),x.IntvBox,x.xref)
    elseif (x.Intv.hi <= y.Intv.lo)
      cc = y.cc
      cc_grad = y.cnst ? zeros(y.cc_grad) : y.cc_grad
    elseif (x.Intv.lo >= y.Intv.hi)
      cc = x.cc
      cc_grad = x.cnst ? zeros(x.cc_grad) : x.cc_grad
    elseif (MC_param.multivar_refine)
      maxLL::T = max(x.Intv.lo,y.Intv.lo)
      maxLU::T = max(x.Intv.lo,y.Intv.hi)
      maxUL::T = max(x.Intv.hi,y.Intv.lo)
      maxUU::T = max(x.Intv.hi,y.Intv.hi)
      thin1::Bool = (diam(x.Intv) == zero(T))
      thin2::Bool = (diam(y.Intv) == zero(T))
      r11::T = thin1 ? zero(T) : (maxUL-maxLL)/diam(x.Intv)
      r21::T = thin1 ? zero(T) : (maxLU-maxUU)/diam(x.Intv)
      r12::T = thin2 ? zero(T) : (maxLU-maxLL)/diam(y.Intv)
      r22::T = thin2 ? zero(T) : (maxUL-maxUU)/diam(y.Intv)
      cc1::T = maxLL + r11*(x.cc-x.Intv.lo) + r12*(y.cc-y.Intv.lo)
      cc2::T = maxUU - r21*(x.cc-x.Intv.hi) - r22*(y.cc-y.Intv.hi)
      if (cc1 <= cc2)
        cc = cc1
        cc_grad =  (x.cnst ? zeros(y.cc_grad) : r11*x.cc_grad) + (y.cnst ? zeros(x.cc_grad) : r12*y.cc_grad)
      else
        cc = cc2
        cc_grad = -(x.cnst ? zeros(y.cc_grad) : r21*x.cc_grad) - (y.cnst ? zeros(x.cc_grad) : r22*y.cc_grad)
      end
    else
      ccMC::SMCg{N,T} = (x+y+abs(x-y))/2
      cc = ccMC.cc
      cc_grad = ccMC.cc_grad
    end
    cv = max(x.cv,y.cv)
    cv_grad = (x.cv > y.cv) ? (x.cnst ? zeros(x.cv_grad): x.cv_grad) :
                              (y.cnst ? zeros(y.cv_grad): y.cv_grad)
    cnst = y.cnst ? x.cnst : (x.cnst ? y.cnst : (x.cnst || y.cnst) )

    return SMCg{N,T}(cc, cv, cc_grad, cv_grad, max(x.Intv,y.Intv),cnst,x.IntvBox,x.xref)
end

@inline function maxcv(x::SMCg{N,T},y::SMCg{N,T}) where {N,T<:AbstractFloat}
        cv::T = zero(T)
        temp_mid::T = zero(T)
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
@inline function maxcc(x,y::SMCg{N,T}) where {N,T<:AbstractFloat}
        cc::T = psir_max(x.cc,y.cc,x.Intv,y.Intv)
end
@inline mincv(x,y::SMCg{N,T}) where {N,T<:AbstractFloat} = - maxcc(-x,-y)
@inline min(x::SMCg{N,T},y::SMCg{N,T}) where {N,T<:AbstractFloat} = -max(-x,-y)
@inline max(x::SMCg{N,T},y::Interval{T}) where {N,T<:AbstractFloat} = max(x,SMCg{N,T}(y))
@inline max(y::Interval{T},x::SMCg{N,T}) where {N,T<:AbstractFloat} = max(x,SMCg{N,T}(y))
@inline min(x::SMCg{N,T},y::Interval{T}) where {N,T<:AbstractFloat} = min(x,SMCg{N,T}(y))
@inline min(y::Interval{T},x::SMCg{N,T}) where {N,T<:AbstractFloat} = min(x,SMCg{N,T}(y))
