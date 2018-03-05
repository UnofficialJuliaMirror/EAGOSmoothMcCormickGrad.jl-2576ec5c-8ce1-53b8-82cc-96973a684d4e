"""
--------------------------------------------------------------------------------
Function: GenExpansionParams
--------------------------------------------------------------------------------
Description: Generates relaxation of state variable at pmid in P
--------------------------------------------------------------------------------
Inputs:
h        function - h(z,p) which implicit defines state function
hj       function - hj(z,p) w.r.t z which implicit defines state function
X        IntervalBox - State variable bounds
P        ntervalBox - Decision variable bounds
pmid     Array - Point at which to generate relaxation
mc_opts  - Options for generating implicit function relaxation
--------------------------------------------------------------------------------
Returns:
The tuple (xa,xA,z):
xa - Lower affine relaxation of the state variable
xA - Upper affine relaxation of the state variable
z  - Affine function in X
--------------------------------------------------------------------------------
"""
function GenExpansionParams(h,hj,X,P,pmid,mc_opts)

  nx = length(X)
  np = length(P)
  szero = @SVector zeros(np)
  sone = @SVector ones(np)
  exp_opt = [nx,np,mc_opts.lambda]
  sto_out = []
  sto_z = []
  sto_x = []

  xL = [X[i].lo for i=1:nx]
  xU = [X[i].hi for i=1:nx]
  #x_mc = [SMCg(xU[i],xL[i],ones(np,1),ones(np,1),@interval(xL[i],xU[i]),true,[],[]) for i=1:nx]
  x_mc =[SMCg{np,Float64}(xU[i],xL[i],szero,szero,@interval(xL[i],xU[i]),false,[∅],[1.0]) for i=1:nx]
  xa_mc = [SMCg{np,Float64}(xL[i],xL[i],szero,szero,@interval(xL[i],xL[i]),false,[∅],[1.0]) for i=1:nx]
  xA_mc = [SMCg{np,Float64}(xU[i],xU[i],szero,szero,@interval(xU[i],xU[i]),false,[∅],[1.0]) for i=1:nx]
  z_mct = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc # should be const initially
  if mc_opts.z_rnd_all == true
    z_mc = Rnd_Out_Z_All(z_mct,mc_opts.z_rnd_all_eps)
  elseif mc_opts.z_rnd == true
    z_mc = Rnd_Out_Z_Intv(z_mct,mc_opts.z_rnd_eps)
  else
    z_mc = z_mct
  end
  push!(sto_out,copy(x_mc))
  p_mc = [SMCg(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  pref_mc = copy(p_mc)

  # setsup aff_mc, H_mc storage
  aff_mc = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
  aff_mct = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]

  H_mc = h(z_mc,p_mc)
  YH_mc = copy(H_mc)
  dH_mc = hj(aff_mc,p_mc)
  YdH_mc = copy(dH_mc)
  J_int = [dH_mc[i,j].Intv for i=1:nx,j=1:nx]
  Jm = J_int
  Y = (size(Jm) == 1) ? [1.0/Jm[1]] : inv(Jm)

  # stores some things for repeated us
  optc = [szero,sone]

  # Begins loop to generate parameters
  push!(sto_z,copy(z_mc))
  push!(sto_x,copy(x_mc))
  for k=1:mc_opts.kmax
    aff_mc = [SMCg{np,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                   @interval(x_mc[i].cv,x_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    if mc_opts.aff_rnd_all == true
      aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
    elseif mc_opts.aff_rnd == true
      aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
    end
    H_mc = h(z_mc,p_mc)
    dH_mc = hj(aff_mc,p_mc)
    if mc_opts.hhj_rnd_all == true
      H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
    elseif mc_opts.hhj_rnd == true
      H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
    end
    J_int = [dH_mc[i,j].Intv for i=1:nx,j=1:nx]
    Jm = mid.(J_int)
    Y = (size(Jm) == 1) ? [1.0/Jm[1]] : inv(Jm)

    # applies preconditioner
    H_mc,dH_mc = Precondition(H_mc,dH_mc,Y,nx)

    # applies parametric iteration
    if (mc_opts.style == "NewtonGS")
      x_mc = MC_NewtonGS(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
    elseif (mc_opts.style == "KrawczykCW")
      x_mc = MC_KrawczykCW(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
    else
      error("Unsupported Style of Implicit Relaxation")
    end
    # update affine relaxations & correct
    push!(sto_x,copy(x_mc))
    xa,xA,z_mc = Affine_Exp(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,exp_opt)
    Correct_Exp!(z_mc,x_mc,xL,xU,nx,np,mc_opts.aff_correct_eps)
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
    end

    # store relaxation
    push!(sto_out,copy(x_mc))
    push!(sto_z,copy(z_mc))
  end
  return sto_out,sto_z,sto_x
end
