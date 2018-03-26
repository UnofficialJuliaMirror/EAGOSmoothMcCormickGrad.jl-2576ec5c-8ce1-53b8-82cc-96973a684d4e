function +(x::SMCg{N,T},y::SMCg{N,T}) where {N,T<:AbstractFloat}
	return SMCg{N,T}(x.cc+y.cc, x.cv+y.cv, x.cc_grad+y.cc_grad, x.cv_grad+y.cv_grad,
						 (x.Intv+y.Intv),(x.cnst && y.cnst),x.IntvBox,x.xref)
end
function -(x::SMCg{N,T},c::SMCg{N,T}) where {N,T<:AbstractFloat}
	return x + (-c)
end
function -(x::SMCg{N,T}) where {N,T<:AbstractFloat}
	return SMCg{N,T}(-x.cv, -x.cc, -x.cv_grad, -x.cc_grad, -x.Intv,
							 x.cnst, x.IntvBox, x.xref)
end
function /(x::SMCg{N,T},y::SMCg{N,T}) where {N,T<:AbstractFloat}
	pos_orth::Bool = (x.Intv.lo >= 0) && (y.Intv.lo >= 0)
	if (x==y)
		#println("ran me to 1")
		return one(T)
	elseif (MC_param.multivar_refine) && (~ (MC_param.mu >= 1)) && (pos_orth)
		#println("ran me to 2")
		Intv::Interval{T} = x.Intv/y.Intv
		cv1::T,pos1::Int64 = mid3(x.cv,x.cc,x.Intv.lo)
		cv1t::T = (cv1+sqrt(x.intv.lo*x.intv.hi))
				  /(sqrt(x.intv.lo)+sqrt(x.intv.hi))
		cv2t::T,pos2::Int64 = mid3(y.cv,y.cc,y.Intv.hi)
		cv::T = sqr(cv1t)/cv2t
		cv_grad::SVector{N,T} = 2.0*(cv1t/cv2t)/(sqrt(x.Intv.lo)+sqrt(x.Intv.hi))*
								  mid_grad(x.cc_grad, x.cv_grad, pos1) - ((cv1t/cv2t)^2)*
									mid_grad(y.cc_grad, y.cv_grad, pos2)
		cc1::T,pos1 = mid3(x.cv,x.cc,x.Intv.hi)
		cc2::T,pos2 = mid3(y.cv,y.cc,y.Intv.lo)
		gcc1::T = y.Intv.hi*cc1 - x.Intv.lo*cc2 + x.Intv.lo*y.Intv.lo
		gcc2::T = y.Intv.lo*cc1 - x.Intv.hi*cc2 + x.Intv.hi*y.Intv.hi
		if gcc1 <= gcc2
			cc::T = gcc1/(y.Intv.lo*y.Intv.hi)
			cc_grad::SVector{N,T} = one(T)/y.Intv.lo*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.Intv.lo/(y.Intv.lo*y.Intv.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		else
			cc = gcc2/(y.Intv.lo*y.Intv.hi)
			cc_grad = one(T)/y.Intv.hi*mid_grad(x.cc_grad, x.cv_grad, pos1)
							- x.Intv.hi/(y.Intv.lo*y.Intv.hi)*mid_grad(y.cc_grad, y.cv_grad, pos2)
		end
		cnst::Bool = y.cnst ? x.cnst : (x.cnst ? y.cnst : x.cnst || y.cnst)
		return SMCg{T}(cc,cv,cc_grad,cv_grad,Intv,cnst,x.IntvBox,x.xref)
	else
		#println("ran me to 3")
		return x*inv(y)
	end
end

for numtype in union(int_list,float_list)
eval(quote
function +(x::SMCg{N,T},y::$numtype) where {N,T<:AbstractFloat}
	return SMCg{N,T}(x.cc+y, x.cv+y, x.cc_grad, x.cv_grad, (x.Intv+y),
	            x.cnst, x.IntvBox,x.xref)
end
function +(x::$numtype,y::SMCg{N,T}) where {N,T<:AbstractFloat}
	return SMCg{N,T}(x+y.cc, x+y.cv, y.cc_grad, y.cv_grad, (x+y.Intv),
	            y.cnst, y.IntvBox,y.xref)
end
function -(x::SMCg{N,T},c::$numtype) where {N,T<:AbstractFloat}
	return x + (-c)
end
function -(c::$numtype,x::SMCg{N,T}) where {N,T<:AbstractFloat}
	return c + (-x)
end
function *(x::SMCg{N,T},c::$numtype) where {N,T<:AbstractFloat}
	if (c>=zero(c))
		temp1::T = convert(T,c*x.cc)
		temp2::T = convert(T,c*x.cv)
		return SMCg{N,T}(temp1,temp2,c*x.cc_grad,c*x.cv_grad,c*x.Intv,x.cnst,x.IntvBox,x.xref)
	elseif (c<zero(c))
		temp1 = convert(T,c*x.cv)
		temp2 = convert(T,c*x.cc)
		return SMCg{N,T}(temp1,temp2,c*x.cv_grad,c*x.cc_grad,c*x.Intv,x.cnst,x.IntvBox,x.xref)
	end
end
function *(c::$numtype,x::SMCg{N,T}) where {N,T<:AbstractFloat}
	return x*c
end
@inline sqr(x::$numtype) = x*x
function /(x::SMCg{N,T},y::$numtype) where {N,T<:AbstractFloat}
	x*inv(y)
end
function /(x::$numtype,y::SMCg{N,T}) where {N,T<:AbstractFloat}
	x*inv(y)
end
end)
end
# convex relaxation of square (Khan2016)
function cv_sqr(xMC2::T,xL::T,xU::T) where {T<:AbstractFloat}
  if (zero(xMC2)<=xL||xU<=zero(xMC2))
		  return xMC2^2,2*xMC2
  elseif ((zero(xMC2)<=xMC2)&&((xL<zero(xL))&&(zero(xU)<xU)))
      return (xMC2^3)/xU,(3*xMC2^2)/xU
  elseif (zero(xMC2)>xMC2)&&((xL<zero(xL))&&(zero(xU)<xU))
      return (xMC2^3)/xL,(3*xMC2^2)/xL
  end
end
# concave relaxation of square (Khan2016)
function cc_sqr(xMC1::T,xL::T,xU::T) where {T<:AbstractFloat}
  return xL^2+(xL+xU)*(xMC1-xL),(xL+xU)
end

function sqr_cc_NS(x::T,lo::T,hi::T) where {T<:AbstractFloat}
  return line_seg(x,lo,lo^2,hi,hi^2),dline_seg(x,lo,lo^2,hi,hi^2,2*x)
end
function sqr_cv_NS(x::T,lo::T,hi::T) where {T<:AbstractFloat}
  return x^2,2*x
end
function sqr(x::SMCg{N,T}) where {N,T<:AbstractFloat}
  	eps_min,blank = mid3(x.Intv.lo,x.Intv.hi,zero(x.Intv.lo))
  	eps_max = ifelse((abs(x.Intv.lo)>=abs(x.Intv.hi)),x.Intv.lo,x.Intv.hi)
	midcc,cc_id = mid3(x.cc,x.cv,eps_max)
	midcv,cv_id = mid3(x.cc,x.cv,eps_min)
	if (MC_param.mu >= 1)
  	cc,dcc = cc_sqr(midcc,x.Intv.lo,x.Intv.hi)
  	cv,dcv = cv_sqr(midcv,x.Intv.lo,x.Intv.hi)
		gcc1,gdcc1 = cc_sqr(x.cv,x.Intv.lo,x.Intv.hi)
		gcv1,gdcv1 = cv_sqr(x.cv,x.Intv.lo,x.Intv.hi)
		gcc2,gdcc2 = cc_sqr(x.cc,x.Intv.lo,x.Intv.hi)
		gcv2,gdcv2 = cv_sqr(x.cc,x.Intv.lo,x.Intv.hi)
		cv_grad = max(zero(T),gdcv1)*x.cv_grad + min(zero(T),gdcv2)*x.cc_grad
		cc_grad = min(zero(T),gdcc1)*x.cv_grad + max(zero(T),gdcc2)*x.cc_grad
	else
		cc,dcc = sqr_cc_NS(midcc,x.Intv.lo,x.Intv.hi)
		cv,dcv = sqr_cv_NS(midcv,x.Intv.lo,x.Intv.hi)
		cc_grad = mid_grad(x.cc_grad, x.cv_grad, cc_id)*dcc
		cv_grad = mid_grad(x.cc_grad, x.cv_grad, cv_id)*dcv
	end
  return SMCg{N,T}(cc, cv, cc_grad, cv_grad, x.Intv^2, x.cnst, x.IntvBox,x.xref)
end

########### Defines functions required for linear algebra packages

one(x::SMCg{N,T}) where {N,T<:AbstractFloat} = SMCg{N,T}(one(Float64),one(Float64),zero.(x.cc_grad),zero.(x.cv_grad),Interval(one(Float64)),x.cnst,x.IntvBox,x.xref)

zero(x::SMCg{N,T}) where {N,T<:AbstractFloat} = SMCg{N,T}(zero(x.cc),zero(x.cv),zeros(SVector{N}),zeros(SVector{N}),Interval(zero(x.Intv.lo)),x.cnst,x.IntvBox,x.xref)

real(x::SMCg{N,T}) where {N,T<:AbstractFloat} = x

dist(x1::SMCg{N,T}, x2::SMCg{N,T}) where {N,T<:AbstractFloat} = max(abs(x1.cc-x2.cc), abs(x1.cv-x2.cv))

eps(x::SMCg{N,T}) where {N,T<:AbstractFloat} = max(eps(x.cc), eps(x.cv))

#=
###### Defines boolean operators

"""==(x::SMC,y::SMC) defines == for SMC type
"""
function ==(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv == y.Intv && x.cv == y.cv && x.cc == y.cc
end
"""!=(x::SMC,y::SMC) defines != for SMC type
"""
function !=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv != y.Intv || x.cv != y.cv || x.cc != y.cc
end
"""<=(x::SMC,y::SMC) defines <= for SMC type
"""
function <=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv <= y.Intv && x.cv <= y.cv && x.cc <= y.cc
end
"""=>(x::SMC,y::SMC) defines => for SMC type
"""
function >=(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv >= y.Intv && x.cv >= y.cv && x.cc >= y.cc
end
""">(x::SMC,y::SMC) defines > for SMC type
"""
function >(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv > y.Intv && x.cv > y.cv && x.cc > y.cc
end
"""<(x::SMC,y::SMC) defines < for SMC type
"""
function <(x::SMCg{N,T},y::SMCg{N,T})
	x.Intv < y.Intv && x.cv < y.cv && x.cc < y.cc
end
=#
