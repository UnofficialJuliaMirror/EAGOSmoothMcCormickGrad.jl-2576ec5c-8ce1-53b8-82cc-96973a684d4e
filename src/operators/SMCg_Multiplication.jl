########### Defines multiplicaiton of SMC and SMC
function sigu(x)
	 (zero(x)<=x) ? exp(log(x)/MC_param.mu) : -exp(log(abs(x))/MC_param.mu)
 end
function xstar(y,lambda::Interval,nu::Interval)
	return lambda.lo+(lambda.hi-lambda.lo)*(((nu.hi-y)/(nu.hi-nu.lo))
				 -sigu((nu.lo+nu.hi)/((MC_param.mu+1.0)*(nu.hi-nu.lo))))
end

function ystar(x,lambda::Interval,nu::Interval)
	return nu.lo+(nu.hi-nu.lo)*(((lambda.hi-x)/(lambda.hi-lambda.lo))
				 -sigu((lambda.lo+lambda.hi)/((MC_param.mu+1.0)*(lambda.hi-lambda.lo))))
end
@inline function dxA(x,y,lambda::Interval,nu::Interval)
	      return (lambda.hi-lambda.lo)*(nu.hi-nu.lo)*
        abs((y-nu.lo)/(nu.hi-nu.lo)-(lambda.hi-x)/(lambda.hi-lambda.lo))^(MC_param.mu+1.0)
end
@inline function GxA(x,y,lambda::Interval,nu::Interval)
	      return 0.5*(x*(nu.lo+nu.hi)+y*(lambda.lo+lambda.hi)-
        (lambda.lo*nu.lo+lambda.hi*nu.hi)+dxA(x,y,lambda,nu))
end
@inline function dxB(x,y,lambda::Interval,nu::Interval)
	      return (lambda.hi-lambda.lo)*(nu.hi-nu.lo)*max(zero(y),(((y-nu.lo)/
        (nu.hi-nu.lo))-(lambda.hi-x)/(lambda.hi-lambda.lo)))^(Int64(MC_param.mu+1))
end
@inline function GxB(x,y,lambda::Interval,nu::Interval)
	      return x*nu.lo+y*lambda.lo-lambda.lo*nu.lo+dxB(x,y,lambda,nu)
end
@inline function gCxA(alpha::Interval,beta::Interval,lambda::Interval,nu::Interval,x1,x2)
        term1,blank = mid3(alpha.lo,alpha.hi,xstar(beta.lo,lambda,nu))
	      term2,blank = mid3(alpha.lo,alpha.hi,xstar(beta.hi,lambda,nu))
	      term3,blank = mid3(beta.lo,beta.hi,ystar(alpha.lo,lambda,nu))
	      term4,blank = mid3(beta.lo,beta.hi,ystar(alpha.hi,lambda,nu))
	      a,b = findmin([GxA(term1,beta.lo,lambda,nu),GxA(term2,beta.hi,lambda,nu),
                           GxA(alpha.lo,term3,lambda,nu),GxA(alpha.hi,term4,lambda,nu)])
        if (b == 1)
			  	grad = max(0,psi_mlt_Ax(term1,beta.lo,lambda,nu))*x1.cv_grad+
						     min(0,psi_mlt_Ax(term1,beta.lo,lambda,nu))*x1.cc_grad+
					 		   max(0,psi_mlt_Ay(term1,beta.lo,lambda,nu))*x2.cv_grad+
				 			   min(0,psi_mlt_Ay(term1,beta.lo,lambda,nu))*x2.cc_grad
				elseif (b == 2)
					grad = max(0,psi_mlt_Ax(term2,beta.hi,lambda,nu))*x1.cv_grad+
						     min(0,psi_mlt_Ax(term2,beta.hi,lambda,nu))*x1.cc_grad+
					 		   max(0,psi_mlt_Ay(term2,beta.hi,lambda,nu))*x2.cv_grad+
				 			   min(0,psi_mlt_Ay(term2,beta.hi,lambda,nu))*x2.cc_grad
				elseif (b == 3)
					grad = max(0,psi_mlt_Ax(alpha.lo,term3,lambda,nu))*x1.cv_grad+
						     min(0,psi_mlt_Ax(alpha.lo,term3,lambda,nu))*x1.cc_grad+
					 		   max(0,psi_mlt_Ay(alpha.lo,term3,lambda,nu))*x2.cv_grad+
				 			   min(0,psi_mlt_Ay(alpha.lo,term3,lambda,nu))*x2.cc_grad
				else
					grad = max(0,psi_mlt_Ax(alpha.hi,term4,lambda,nu))*x1.cv_grad+
						     min(0,psi_mlt_Ax(alpha.hi,term4,lambda,nu))*x1.cc_grad+
					 		   max(0,psi_mlt_Ay(alpha.hi,term4,lambda,nu))*x2.cv_grad+
				 			   min(0,psi_mlt_Ay(alpha.hi,term4,lambda,nu))*x2.cc_grad
				end

        return a,grad
end

function psi_mlt_Ax(x,y,lambda,nu)
	term = [(y-nu.lo)/(nu.hi-nu.lo)
	        (lambda.hi-x)/(lambda.hi-lambda.lo)]
	return 0.5*(nu.lo+nu.hi+(MC_param.mu+1)*(nu.hi-nu.lo)*(term[1]-term[2])*abs(term[1]-term[2])^(MC_param.mu-1))
end
function psi_mlt_Ay(x,y,lambda,nu)
	term = [(y-nu.lo)/(nu.hi-nu.lo)
	        (lambda.hi-x)/(lambda.hi-lambda.lo)]
	return 0.5*(lambda.lo+lambda.hi+(MC_param.mu+1)*(lambda.hi-lambda.lo)*(term[1]-term[2])*abs(term[1]-term[2])^(MC_param.mu-1))
end
function psi_mlt_Bx(x,y,lambda,nu)
	term = [(y - nu.lo)/(nu.hi-nu.lo)
	        (lambda.hi - x)/(lambda.hi-lambda.lo)]
	return nu.lo+(MC_param.mu+1)*(nu.hi-nu.lo)*max(0,term[1]-term[2])^MC_param.mu
end
function psi_mlt_By(x,y,lambda,nu)
	term = [(y - nu.lo)/(nu.hi-nu.lo)
	        (lambda.hi - x)/(lambda.hi-lambda.lo)]
	return lambda.lo+(MC_param.mu+1)*(lambda.hi-lambda.lo)*max(0,term[1]-term[2])^MC_param.mu
end

function multiply_MV(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	cv = zero(x1.cc)
	cc = zero(x1.cc)
	alpha0 = Interval(x1.cv,x1.cc)
	beta0 =  Interval(x2.cv,x2.cc)
	#println("x1: ", x1)
	#println("x2: ", x2)
	if (zero(x1.Intv.lo)<=x1.Intv.lo) && (zero(x2.Intv.lo)<=x2.Intv.lo)
		cv = GxB(x1.cv,x2.cv,x1.Intv,x2.Intv)
		cv_grad = x1.cv_grad*psi_mlt_Bx(x1.cv,x2.cv,x1.Intv,x2.Intv)
						+ x2.cv_grad*psi_mlt_By(x1.cv,x2.cv,x1.Intv,x2.Intv)
	elseif ((x1.Intv.hi<=zero(x1.Intv.hi))) && (x2.Intv.hi<=zero(x2.Intv.hi))
		cv = GxB(-x1.cc,-x2.cc,-x1.Intv,-x2.Intv)
		cv_grad = - x1.cc_grad*psi_mlt_Bx(-x1.cc,-x2.cc,-x1.Intv,-x2.Intv)
						- x2.cc_grad*psi_mlt_By(-x1.cc,-x2.cc,-x1.Intv,-x2.Intv)
	else
		cv,cv_grad = gCxA(alpha0,beta0,x1.Intv,x2.Intv,x1,x2)
	end
	if ((x1.Intv.hi<=zero(x1.Intv.hi))) && ((zero(x2.Intv.lo))<=x2.Intv.lo)
		cc = -GxB(-x1.cc,x2.cv,-x1.Intv,x2.Intv)
		cc_grad = x1.cc_grad*psi_mlt_Bx(-x1.cc,x2.cv,-x1.Intv,x2.Intv)
						- x2.cv_grad*psi_mlt_By(-x1.cc,x2.cv,-x1.Intv,x2.Intv)
	elseif ((zero(x1.Intv.lo))<=x1.Intv.lo) && (x2.Intv.hi<=zero(x2.Intv.hi))
		cc = -GxB(x1.cv,-x2.cc,x1.Intv,-x2.Intv)
		cc_grad = - x1.cv_grad*psi_mlt_Bx(x1.cv,-x2.cc,x1.Intv,-x2.Intv)
							+ x2.cc_grad*psi_mlt_By(x1.cv,-x2.cc,x1.Intv,-x2.Intv)
	else
		cct,cc_gradt = gCxA(-alpha0,beta0,-x1.Intv,x2.Intv,x1,x2)
		cc = -cct
		cc_grad = -cc_gradt
	end
	if (min(x1.Intv.lo,x2.Intv.lo)<zero(x1.Intv.lo)<max(x1.Intv.hi,x2.Intv.hi))
		#println("interval calc: ")
		lo_Intv_calc,blank = gCxA(x1.Intv,x2.Intv,x1.Intv,x2.Intv,x1,x2)
		hi_Intv_calct,blankt = gCxA(-x1.Intv,x2.Intv,-x1.Intv,x2.Intv,x1,x2)
		hi_Intv_calc = -hi_Intv_calct
		Intv_calc = @interval(lo_Intv_calc,hi_Intv_calc)
	else
		#println("alt interval calc: ")
		Intv_calc = x1.Intv*x2.Intv
	end
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	#println("cc: ", cc)
	#println("cv: ", cv)
	#println("Intv_calc: ", Intv_calc)
	return SMCg{N,T}(cc, cv, cc_grad, cv_grad, Intv_calc, cnst, x1.IntvBox,x1.xref)
end

function mul1_u1pos_u2pos(x1::SMCg{N,T},x2::SMCg{N,T},cnst) where {N,T}
  Intv = x1.Intv*x2.Intv
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cv_grad
  end

  cc1 = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cc_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul1_u1pos_u2mix(x1::SMCg{N,T},x2::SMCg{N,T},cnst) where {N,T}
  Intv = x1.Intv*x2.Intv
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end

  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
# ERROR KrawczykCW, NewtonGS (2x)
function mul1_u1mix_u2mix(x1::SMCg{N,T},x2::SMCg{N,T},cnst) where {N,T}
  Intv = x1.Intv*x2.Intv
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
  cv = cv1
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x2.Intv.hi*x1.cv_grad
  else
    cv = cv2
    cv_grad = x2.Intv.lo*x1.cc_grad
  end
  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
  cc = cc1
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x2.Intv.lo*x1.cv_grad
  else
    cc = cc2
    cc_grad = x2.Intv.hi*x1.cc_grad
  end
  return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul2_u1pos_u2pos(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	#println("x1:   ", x1)
	#println("x2:   ", x2)
	Intv = x1.Intv*x2.Intv
	cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.lo*x2.cv_grad
	end

	cc1 = x2.Intv.lo*x1.cc + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end
	#println("cv1:   ", cv)
	#println("cc:   ", cc)
	return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv, cnst,x1.IntvBox,x1.xref)
end
function mul2_u1pos_u2mix(x1::SMCg{N,T},x2::SMCg{N,T},cnst) where {N,T}
  Intv = x1.Intv*x2.Intv
  cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
  cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
  if (cv1 > cv2)
    cv = cv1
    cv_grad = x1.Intv.hi*x2.cv_grad
  else
    cv = cv2
    cv_grad = x1.Intv.lo*x2.cc_grad
  end

  cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
  cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
  if (cc1 < cc2)
    cc = cc1
    cc_grad = x1.Intv.hi*x2.cc_grad
  else
    cc = cc2
    cc_grad = x1.Intv.lo*x2.cc_grad
  end
  return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul2_u1mix_u2mix(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	Intv = x1.Intv*x2.Intv
  cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end

	cc1 = x2.Intv.lo*x1.cv + x1.Intv.lo*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end

	return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end
function mul3_u1pos_u2mix(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	Intv = x1.Intv*x2.Intv
  cnst = x2.cnst ? x1.cnst : (x1.cnst ? x2.cnst : x1.cnst || x2.cnst)
	cv1 = x2.Intv.hi*x1.cv + x1.Intv.hi*x2.cv - x1.Intv.hi*x2.Intv.hi
	cv2 = x2.Intv.lo*x1.cc + x1.Intv.lo*x2.cv - x1.Intv.lo*x2.Intv.lo
	if (cv1 > cv2)
		cv = cv1
		cv_grad = x2.Intv.hi*x1.cv_grad + x1.Intv.hi*x2.cv_grad
	else
		cv = cv2
		cv_grad = x2.Intv.lo*x1.cc_grad + x1.Intv.lo*x2.cv_grad
	end

	cc1 = x2.Intv.lo*x1.cv + x1.Intv.hi*x2.cc - x1.Intv.hi*x2.Intv.lo
	cc2 = x2.Intv.hi*x1.cc + x1.Intv.lo*x2.cc - x1.Intv.lo*x2.Intv.hi
	if (cc1 < cc2)
		cc = cc1
		cc_grad = x2.Intv.lo*x1.cv_grad + x1.Intv.hi*x2.cc_grad
	else
		cc = cc2
		cc_grad = x2.Intv.hi*x1.cc_grad + x1.Intv.lo*x2.cc_grad
	end

	return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x1.IntvBox,x1.xref)
end

function mul_MV_ns1cv(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return MC2.Intv.hi*x1+MC1.Intv.hi*x2-MC2.Intv.hi*MC1.Intv.hi
end
function mul_MV_ns2cv(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return MC2.Intv.lo*x1+MC1.Intv.lo*x2-MC2.Intv.lo*MC1.Intv.lo
end
function mul_MV_ns3cv(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return max(mul_MV_ns1cv(x1,x2,MC1,MC2),mul_MV_ns2cv(x1,x2,MC1,MC2))
end
function mul_MV_ns1cc(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return MC2.Intv.lo*x1+MC1.Intv.hi*x2-MC2.Intv.lo*MC1.Intv.hi
end
function mul_MV_ns2cc(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return MC2.Intv.hi*x1+MC1.Intv.lo*x2-MC2.Intv.hi*MC1.Intv.lo
end
function mul_MV_ns3cc(x1,x2,MC1::SMCg{N,T},MC2::SMCg{N,T}) where {N,T}
	return min(mul_MV_ns1cc(x1,x2,MC1,MC2),mul_MV_ns2cc(x1,x2,MC1,MC2))
end
function multiply_MV_NS(x1::SMCg{N,T},x2::SMCg{N,T},ngrad,cnst) where {N,T}
 # CV calculation
 k = diam(x2.Intv)/diam(x1.Intv)
 z = (x1.Intv.hi*x2.Intv.hi - x1.Intv.lo*x2.Intv.lo)/diam(x1.Intv)
 x1vta,blank = mid3(x1.cv,x1.cc,(x2.cv-z)/k)
 x1vtb,blank = mid3(x1.cv,x1.cc,(x2.cc-z)/k)
 x2vta,blank = mid3(x2.cv,x2.cc, k*x1.cv+z)
 x2vtb,blank = mid3(x2.cv,x2.cc, k*x1.cc+z)
 x1vt = [x1.cv, x1.cc, x1vta, x1vtb, x1.cv, x1.cc]
 x2vt = [x2vta, x2vtb, x2.cv, x2.cc, x2.cv, x2.cc]
 vt  = [mul_MV_ns3cv(x1vt[1],x2vt[1],x1,x2), mul_MV_ns3cv(x1vt[2],x2vt[2],x1,x2),
        mul_MV_ns3cv(x1vt[3],x2vt[3],x1,x2), mul_MV_ns3cv(x1vt[4],x2vt[4],x1,x2),
				mul_MV_ns3cv(x1vt[5],x2vt[5],x1,x2), mul_MV_ns3cv(x1vt[6],x2vt[6],x1,x2)]
 cv,cvind = findmax(vt)

 if (ngrad>0)
	if isequal(mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
 		alph = [0.0,1.0]

 		MC1thin = isequal(x1.cv,x1.cc,MC_param.mv_tol,MC_param.mv_tol)
 		if ((~MC1thin) && (x1vt[cvind] > x1.cv))
	 		if (~isequal(x1vt[cvind],x1.cv,MC_param.mv_tol,MC_param.mv_tol))
		 		alph[2] = min(alph[2],-x2.Intv.lo/diam(x2.Intv))
	 		end
 		end
 		if ((~MC1thin) && (x1vt[cvind] < x1.cc))
			if (~isequal(x1vt[cvind],x1.cc,MC_param.mv_tol,MC_param.mv_tol))
				alph[1] = max(alph[1],-x2.Intv.lo/diam(x2.Intv))
			end
 		end

 		MC2thin = isequal(x2.cv,x2.cc,MC_param.mv_tol,MC_param.mv_tol)
 		if ((~MC2thin) && (x2vt[cvind] > x2.cv))
			if (~isequal(x2vt[cvind],x2.cv,MC_param.mv_tol,MC_param.mv_tol))
				alph[2] = min(alph[2],-x1.Intv.lo/diam(x1.Intv))
			end
 		end
 		if ((~MC2thin) && (x2vt[cvind] < x2.cc))
 			if (~isequal(x2vt[cvind],x2.cc,MC_param.mv_tol,MC_param.mv_tol))
	 			alph[1] = max(alph[1],-x1.Intv.lo/diam(x1.Intv))
 			end
 		end

 		alphthin = isequal(alph[1],alph[2],MC_param.mv_tol,MC_param.mv_tol)
 		if (~alphthin && (alph[1]>alph[2]))
	 		error("Multivariant mult error alphaL = alphaU")
 		end
 		myalph = 0.5*(alph[1]+alph[2])
 	elseif (mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2) >
			    mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2))
		myalph = 1.0
	else
		myalph = 0.0
	end
	sigma_cv1 = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cv2 = x1.Intv.lo + myalph*diam(x1.Intv)
	if (x1.cnst)
		term1 = zero(N)
	elseif (sigma_cv1>=0.0)
		term1 = x1.cv_grad
	else
		term1 = x1.cc_grad
	end
	if (x2.cnst)
		term2 = zero(N)
	elseif (sigma_cv1>=0.0)
		term2 = x2.cv_grad
	else
		term2 = x2.cc_grad
	end
	cv_grad = term1*sigma_cv1 + term2*sigma_cv2
 end

 # CC Calculation
 z = (x1.Intv.hi*x2.Intv.lo - x1.Intv.lo*x2.Intv.hi)/diam(x1.Intv)
 x1vta,blank = mid3(x1.cv,x1.cc,(x2.cv-z)/k)
 x1vtb,blank = mid3(x1.cv,x1.cc,(x2.cc-z)/k)
 x2vta,blank = mid3(x2.cv,x2.cc, k*x1.cv+z)
 x2vtb,blank = mid3(x2.cv,x2.cc, k*x1.cc+z)
 x1vt = [x1.cv, x1.cc, x1vta, x1vtb, x1.cv, x1.cc]
 x2vt = [x2vta, x2vtb, x2.cv, x2.cc, x2.cc, x2.cv]
 vt  = [mul_MV_ns3cc(x1vt[1],x2vt[1],x1,x2), mul_MV_ns3cc(x1vt[2],x2vt[2],x1,x2),
				mul_MV_ns3cc(x1vt[3],x2vt[3],x1,x2), mul_MV_ns3cc(x1vt[4],x2vt[4],x1,x2),
			  mul_MV_ns3cc(x1vt[5],x2vt[5],x1,x2), mul_MV_ns3cc(x1vt[6],x2vt[6],x1,x2)]
 cc,ccind = findmax(vt)

 if (ngrad>0)
 	if isequal(mul_MV_ns1cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 mul_MV_ns2cc(x1vt[cvind],x2vt[cvind],x1,x2),
						 MC_param.mv_tol,MC_param.mv_tol)
		 alph = [0.0,1.0]

		 MC1thin = isequal(x1.cv,x1.cc,MC_param.mv_tol,MC_param.mv_tol)
		 if ((~MC1thin) && (x1vt[cvind] > x1.cv))
		 	if (~isequal(x1vt[cvind],x1.cv,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[1] = max(alph[1],-x2.Intv.lo/diam(x2.Intv))
		 	end
		 end
		 if ((~MC1thin) && (x1vt[cvind] < x1.cc))
		 	if (~isequal(x1vt[cvind],x1.cc,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[2] = min(alph[2],-x2.Intv.lo/diam(x2.Intv))
		 	end
		 end

		 MC2thin = isequal(x2.cv,x2.cc,MC_param.mv_tol,MC_param.mv_tol)
		 if ((~MC2thin) && (x2vt[cvind] > x2.cv))
		 	if (~isequal(x2vt[cvind],x2.cv,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[2] = min(alph[2],x1.Intv.hi/diam(x1.Intv))
		 	end
		 end
		 if ((~MC2thin) && (x2vt[cvind] < x2.cc))
			if (~isequal(x2vt[cvind],x2.cc,MC_param.mv_tol,MC_param.mv_tol))
			 	alph[1] = max(alph[1],x1.Intv.hi/diam(x1.Intv))
			end
		 end

		 alphthin = isequal(alph[1],alph[2],MC_param.mv_tol,MC_param.mv_tol)
		 if (~alphthin && (alph[1]>alph[2]))
		 	error("Multivariant mult error alphaL = alphaU")
		 end
		 myalph = 0.5*(alph[1]+alph[2])
	elseif (mul_MV_ns1cv(x1vt[cvind],x2vt[cvind],x1,x2) >
				  mul_MV_ns2cv(x1vt[cvind],x2vt[cvind],x1,x2))
	   myalph = 1.0
  else
	 myalph = 0.0
  end
	sigma_cc1 = x2.Intv.lo + myalph*diam(x2.Intv)
	sigma_cc2 = x1.Intv.hi - myalph*diam(x1.Intv)
	if (x1.cnst)
		term1 = zero(N)
	elseif (sigma_cc1>=0.0)
		term1 = x1.cc_grad
	else
		term1 = x1.cv_grad
	end
	if (x2.cnst)
		term2 = zero(N)
	elseif (sigma_cc1>=0.0)
		term2 = x2.cc_grad
	else
		term2 = x2.cv_grad
	end
	cv_grad = term1*sigma_cc1 + term2*sigma_cc2
 end
 return SMCg{N,T}(cc,cv,cc_grad,cv_grad,x1.Intv*x2.Intv,cnst,x1.IntvBox,x1.xref)
end

function multiply_STD_NS(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	if (x2.Intv.lo >=0)
    	if (x2.cnst)
      		return mul1_u1pos_u2pos(x1,x2,x1.cnst) # cnst to do
    	elseif (x1.cnst)
      		return mul1_u1pos_u2pos(x2,x1,x2.cnst) # cnst to do
    	else
      		return mul2_u1pos_u2pos(x1,x2) # DONE
    	end
	elseif (x2.Intv.hi <=0)
		return -(x1*(-x2))
	else
    	if (x2.cnst)
      		return mul1_u1pos_u2mix(x1,x2,x1.cnst) # cnst to do
    	elseif (x1.cnst)
      		return mul2_u1pos_u2mix(x1,x2,x2.cnst) # cnst to do
    	else
	  		return mul3_u1pos_u2mix(x1,x2) # cnst to do
    	end
	end
end

function STD_NS_ALT(x::SMCg{N,T},y::SMCg{N,T}) where {N,T}
	alpha1 = min( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	alpha2 = min( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )
	beta1  = min( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	beta2  = min( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	gamma1 = max( y.Intv.lo*x.cv,  y.Intv.lo*x.cc )
	gamma2 = max( x.Intv.hi*y.cv,  x.Intv.hi*y.cc )
	delta1 = max( y.Intv.hi*x.cv,  y.Intv.hi*x.cc )
	delta2 = max( x.Intv.lo*y.cv,  x.Intv.lo*y.cc )

	cv1 = alpha1 + alpha2 - x.Intv.lo*y.Intv.lo
	cv2 = beta1  + beta2  - x.Intv.hi*y.Intv.hi
	cc1 = gamma1 + gamma2 - x.Intv.hi*y.Intv.lo
	cc2 = delta1 + delta2 - x.Intv.lo*y.Intv.hi
	#println("cv1: ", cv1)
	#println("cv2: ", cv2)

	s_alpha1 = (y.Intv.lo >= 0.0) ? y.Intv.lo*x.cv_grad : y.Intv.lo*x.cc_grad
	s_alpha2 = (x.Intv.lo >= 0.0) ? x.Intv.lo*x.cv_grad : x.Intv.lo*x.cc_grad
	s_beta1  = (y.Intv.hi >= 0.0) ? y.Intv.hi*x.cv_grad : y.Intv.hi*x.cc_grad
	s_beta2  = (x.Intv.hi >= 0.0) ? x.Intv.hi*x.cv_grad : x.Intv.hi*x.cc_grad
	s_gamma1 = (y.Intv.lo >= 0.0) ? y.Intv.lo*x.cc_grad : y.Intv.lo*x.cv_grad
	s_gamma2 = (x.Intv.hi >= 0.0) ? x.Intv.hi*x.cc_grad : x.Intv.hi*x.cv_grad
	s_delta1 = (y.Intv.hi >= 0.0) ? y.Intv.hi*x.cc_grad : y.Intv.hi*x.cv_grad
	s_delta2 = (x.Intv.lo >= 0.0) ? x.Intv.lo*x.cc_grad : x.Intv.lo*x.cv_grad

	if (cv1 >= cv2)
	    cv = cv1
	    cv_grad = s_alpha1 + s_alpha2
	else
	    cv = cv2
	    cv_grad = s_beta1 + s_beta2
	end
	#println("cv: ", cv)

	if (cc1 <= cc2)
	    cc = cc1
	    cc_grad = s_gamma1 + s_gamma2
	else
	    cc = cc2
	    cc_grad = s_delta1 + s_delta2
	end
	Intv = x.Intv*y.Intv
	return SMCg{N,T}(cc,cv,cc_grad,cv_grad,Intv,(x.cnst && y.cnst),x.IntvBox,x.xref)
end

@inline function *(x1::SMCg{N,T},x2::SMCg{N,T}) where {N,T}
	if x1 == x2
		return sqr(x1)
	end

	degen1::Bool = ((x1.Intv.hi - x1.Intv.lo) == 0.0)
	degen2::Bool = ((x2.Intv.hi - x2.Intv.lo) == 0.0)

	if (MC_param.mu >= 1 && ~(degen1||degen2))
		return multiply_MV(x1,x2)
	elseif (MC_param.multivar_refine && ~(degen1||degen2))
		if (x2.cnst)
			ngrad,cnst = length(x1.cc_grad),x1.cnst
		elseif (x1.cnst)
			ngrad,cnst = length(x2.cc_grad),x2.cnst
		elseif (length(x1.cc_grad) != length(x2.cc_grad))
			error("Unequal subgradients")
		else
			ngrad,cnst = length(x1.cc_grad),(x1.cnst||x2.cnst)
		end
		return multiply_MV_NS(x1,x2,ngrad,cnst) # DONE (minus gradients & case handling?)
	elseif (x1.Intv.lo >= 0)
		#return multiply_STD_NS(x1,x2)
		return STD_NS_ALT(x1,x2)
	elseif (x1.Intv.hi <= 0)
		if (x2.Intv.lo >= 0)
			return -((-x1)*x2)
		elseif (x2.Intv.hi <= 0)
			return (-x1)*(-x2)
		else
			return -(x2*(-x1))
		end
	elseif (x2.Intv.lo >= 0)
		return x2*x1
	elseif (x2.Intv.hi <= 0)
		return -((-x2)*x1)
	else
    if (x2.cnst)
	  return STD_NS_ALT(x1,x2)
      #return mul1_u1mix_u2mix(x1,x2,x1.cnst)
    elseif (x1.cnst)
	  return STD_NS_ALT(x1,x2)
      #return mul1_u1mix_u2mix(x2,x1,x2.cnst)
    else
	  return STD_NS_ALT(x1,x2)
	  #return mul2_u1mix_u2mix(x1,x2)
    end
	end
end
