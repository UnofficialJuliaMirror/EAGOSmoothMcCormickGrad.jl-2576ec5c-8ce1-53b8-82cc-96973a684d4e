"""
--------------------------------------------------------------------------------
Function: Affine_Exp
--------------------------------------------------------------------------------
Description: Computates the affine relaxations of the state variable
--------------------------------------------------------------------------------
Inputs:
x       Vector{SMCg} - State variable relaxation
p       Vector{SMCg} - Decision variable relaxation
p_ref   Vector{SMCg} - Reference variable relaxation
xa      Vector{SMCg} - Lower affine relaxation of the state variable
xA      Vector{SMCg} - Upper affine relaxation of the state variable
z       Vector{SMCg} - Affine function in X
opt     Array - [np,nx,lambda] values for relaxation
--------------------------------------------------------------------------------
Returns:
The tuple (xa,xA,z):
xa - Lower affine relaxation of the state variable
xA - Upper affine relaxation of the state variable
z  - Affine function in X
--------------------------------------------------------------------------------
"""
function Affine_Exp(x,p,p_ref,xa,xA,z,opt)
  nx,np,lambda = opt
  S1,S2,S3 = 0.0,0.0,0.0

  for i = 1:Int(nx)
   S1 = 0.0
   S2 = 0.0
   S3 = 0.0
   for j = 1:Int(np)
      S1 = S1 + (p[j]-p_ref[j])*x[i].cv_grad[j]
      S2 = S2 + (p[j]-p_ref[j])*x[i].cc_grad[j]
      S3 = S3 + (lambda*x[i].cv_grad[j]+(1.0-lambda)*x[i].cc_grad[j])*(p[j]-p_ref[j])
   end
   temp1 = x[i].cv + S1
   temp2 = x[i].cc + S2
   temp3 = x[i].cv_grad
   temp4 = x[i].cc_grad
   temp5 = lambda*x[i].cv+(1.0-lambda)*x[i].cc+S3
   temp6 = lambda*x[i].cv_grad+(1.0-lambda)*x[i].cc_grad
   cnst_xa = S1.cnst
   cnst_xA = S2.cnst
   cnst_z = S3.cnst
   xa[i] = SMCg(temp1.cc,temp1.cv,temp3,temp3,Interval(temp1.cv,temp1.cc),cnst_xa,x[i].IntvBox,x[i].xref)
   xA[i] = SMCg(temp2.cc,temp2.cv,temp4,temp4,Interval(temp2.cv,temp2.cc),cnst_xA,x[i].IntvBox,x[i].xref)
   z[i] = SMCg(temp5.cc,temp5.cv,temp6,temp6,Interval(temp5.cv,temp5.cc),cnst_z,x[i].IntvBox,x[i].xref)
  end
  return xa,xA,z
end

"""
--------------------------------------------------------------------------------
Function: Correct_Exp!
--------------------------------------------------------------------------------
Description: Corrects the relaxation of the state variable if the affine
             relaxation exceeds the interval bounds.
--------------------------------------------------------------------------------
Inputs:
z_mc    Vector{SMCg} - Affine relaxation
x_mc    Vector{SMCg} - Relaxation of state variable
xL      Vector{Float64} - Lower bound on state vector
xU      Vector{Float64} - Upper bound on state vector
nx      Int64 - Size of the state vector
np      Int64 - Size of the decision vector
epsv    Float64 - Tolerance for checking that subgradient exceeds bound
--------------------------------------------------------------------------------
Returns: Mutates x_mc in place.
--------------------------------------------------------------------------------
"""
function Correct_Exp!(z_mc,x_mc,xL,xU,nx::Int64,np::Int64,epsv::Float64)
  for i = 1:nx
    if (z_mc[i].Intv.lo-epsv < xL[i])
      for j = 1:np
        x_mc[i].cv_grad[j] = 0.0
      end
      x_mc[i] = SMCg(x_mc[i].cc,xL[i],x_mc[i].cc_grad,x_mc[i].cv_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
    if (z_mc[i].Intv.hi+epsv > xU[i])
      for j = 1:np
          x_mc[i].cc_grad[j] = 0.0
      end
    x_mc[i] = SMCg(xU[i],x_mc[i].cv,x_mc[i].cc_grad,x_mc[i].cv_grad,x_mc[i].Intv,true,x_mc[i].IntvBox,x_mc[i].xref)
    end
  end
end
