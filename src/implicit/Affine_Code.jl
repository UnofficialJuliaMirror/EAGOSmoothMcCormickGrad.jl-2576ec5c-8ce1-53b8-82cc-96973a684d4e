"""
    Affine_Exp!(x::Vector{SMCg{N,T}},p::Vector{SMCg{N,T}},p_ref::Vector{SMCg{N,T}},
               xa::Vector{SMCg{N,T}},xA::Vector{SMCg{N,T}},z::Vector{SMCg{N,T}},
               opt::Array{Any})

Computates the affine relaxations of the state variable. Inputs are:
* `x::Vector{SMCg{N,T}}`: State variable relaxation
* `p::Vector{SMCg{N,T}}`: Decision variable relaxation
* `p_ref::Vector{SMCg{N,T}}`: Reference variable relaxation
* `xa::Vector{SMCg{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{SMCg{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{SMCg{N,T}}`: Affine function in `X`
* `opt::Array{Any,1}`: `[np,nx,lambda]` values for relaxation
Returns the tuple `(xa,xA,z)`:
* `xa::Vector{SMCg{N,T}}`: Lower affine relaxation of the state variable
* `xA::Vector{SMCg{N,T}}`: Upper affine relaxation of the state variable
* `z::Vector{SMCg{N,T}}`: Affine function in X
--------------------------------------------------------------------------------
"""
function Affine_Exp!(x::Vector{SMCg{N,T}},p::Vector{SMCg{N,T}},p_ref::Vector{SMCg{N,T}},
                    xa::Vector{SMCg{N,T}},xA::Vector{SMCg{N,T}},z::Vector{SMCg{N,T}},
                    opt::Array{Any,1}) where {N,T}
  nx::Int64,np::Int64,lambda::T = opt
  S1::T,S2::T,S3::T = 0.0,0.0,0.0

  for i = 1:nx
   S1 = 0.0
   S2 = 0.0
   S3 = 0.0
   for j = 1:np
      S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      S3 = S3 + (lambda*x[i].cv_grad[j]+(1.0-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   temp1::T = x[i].cv + S1
   temp2::T = x[i].cc + S2
   temp3::SVector{N,T} = x[i].cv_grad
   temp4::SVector{N,T} = x[i].cc_grad
   temp5::T = lambda*x[i].cv+(1.0-lambda)*x[i].cc+S3
   temp6::SVector{N,T} = lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad
   cnst_xa::Bool = S1.cnst
   cnst_xA::Bool = S2.cnst
   cnst_z::Bool = S3.cnst
   xa[i] = SMCg{N,T}(temp1.cc,temp1.cv,temp3,temp3,Interval(temp1.cv,temp1.cc),cnst_xa,x[i].IntvBox,x[i].xref)
   xA[i] = SMCg{N,T}(temp2.cc,temp2.cv,temp4,temp4,Interval(temp2.cv,temp2.cc),cnst_xA,x[i].IntvBox,x[i].xref)
   z[i] = SMCg{N,T}(temp5.cc,temp5.cv,temp6,temp6,Interval(temp5.cv,temp5.cc),cnst_z,x[i].IntvBox,x[i].xref)
  end
  return xa,xA,z
end

"""
    Correct_Exp!(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},xL::Vector{Float64},
                 xU::Vector{Float64},nx::Int64,np::Int64,epsv::Float64)

Corrects the relaxation of the state variable `x_mc` if the affine relaxation,
'z_mc', exceeds the interval bounds `xL` or `xU`.
* `z_mc::Vector{SMCg{N,T}}`: Affine relaxation
* `x_mc::Vector{SMCg{N,T}}`: Relaxation of state variable
* `xL::Vector{Float64}`: Lower bound on state vector
* `xU::Vector{Float64}`: Upper bound on state vector
* `nx::Int64`: Size of the state vector
* `np::Int64`: Size of the decision vector
* `epsv::Float64`: Tolerance for checking that subgradient exceeds bound
"""
function Correct_Exp!(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                     xL::Vector{Float64},xU::Vector{Float64},
                     nx::Int64,np::Int64,epsv::Float64) where {N,T}
  zero_grad::SVector{N,T} = @SVector zeros(T,N)
  for i = 1:nx
    if (z_mc[i].Intv.lo-epsv < xL[i])
      x_mc[i] = SMCg{N,T}(x_mc[i].cc,xL[i],x_mc[i].cc_grad,zero_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
    if (z_mc[i].Intv.hi+epsv > xU[i])
      x_mc[i] = SMCg(xU[i],x_mc[i].cv,zero_grad,x_mc[i].cv_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
  end
end
