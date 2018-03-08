"""
    MC_impRelax(h,hj,p,pmid,X,P,mc_opts,param)

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are param and the basic parameters of the
fixed point method are `mc_opt`.
"""
function MC_impRelax(h,hj,p,pmid,X,P,mc_opts,param)

    nx::Int64 = length(X)
    np::Int64 = length(P)
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    exp_opt = [nx,np,mc_opts.lambda]
    sto_out = []

    xL = [X[i].lo for i=1:nx]
    xU = [X[i].hi for i=1:nx]
    x_mc = [SMCg(xU[i],xL[i],szero,szero,@interval(xL[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
    #x_mc = param[1]
    xa_mc = [SMCg{np,Float64}(xL[i],xL[i],szero,szero,@interval(xL[i],xL[i]),true,[∅],[1.0]) for i=1:nx]
    xA_mc = [SMCg{np,Float64}(xU[i],xU[i],szero,szero,@interval(xU[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
    z_mct = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc # should be const initially
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mct,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mct,mc_opts.z_rnd_eps)
    else
      z_mc = z_mct
    end
    p_mc = copy(p)
    pref_mc = [SMCg{np,Float64}(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]

    #println("p_mc  :",p_mc)
    # setsup aff_mc, H_mc storage
    aff_mc = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    aff_mct = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    # h DAG graph

    H_mc = h(z_mc,p_mc)
    YH_mc = copy(H_mc)
    #println("isolation point #1")
    dH_mc = hj(aff_mc,p_mc)
    #println("isolation point #2")
    YdH_mc = copy(dH_mc)
    J_int = [dH_mc[i,j].Intv for i=1:nx,j=1:nx]
    Jm = J_int
    Y = (size(Jm) == 1) ? [1.0/Jm[1]] : inv(Jm)

    # stores some things for repeated us
    optc = [szero,sone]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      xa_mc,xA_mc,z_mc = Affine_Exp(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,exp_opt)
      if mc_opts.z_rnd_all == true
        z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
      elseif mc_opts.z_rnd == true
        z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
      end

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     @interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
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
      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
      end
      # applies parametric iteration
      if (mc_opts.style == "NewtonGS")
        x_mc = MC_NewtonGS(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
      elseif (mc_opts.style == "KrawczykCW")
        x_mc = MC_KrawczykCW(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
      else
        error("Unsupported Style of Implicit Relaxation")
      end
    end
    return x_mc
end

"""
    MC_impRelax!(h,hj,p,pmid,X,P,mc_opts,param)

Relaxes the implicit function determined by `h(x,p)` with `x` in `X` and `p` in
`P`. The reference point for the affine relaxations is `pmid`. The parameters
generated from the relaxation at `pmid` are param and the basic parameters of the
fixed point method are `mc_opt`. Performs operations in place for jacobian and
preconditioning calculations.
"""
function MC_impRelax!(h,hj,p,pmid,X,P,mc_opts,param)

    nx::Int64 = length(X)
    np::Int64 = length(P)
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    exp_opt = [nx,np,mc_opts.lambda]
    sto_out = []

    xL = [X[i].lo for i=1:nx]
    xU = [X[i].hi for i=1:nx]
    x_mc = [SMCg(xU[i],xL[i],szero,szero,@interval(xL[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
    #x_mc = param[1]
    xa_mc = [SMCg{np,Float64}(xL[i],xL[i],szero,szero,@interval(xL[i],xL[i]),true,[∅],[1.0]) for i=1:nx]
    xA_mc = [SMCg{np,Float64}(xU[i],xU[i],szero,szero,@interval(xU[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
    z_mct = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc # should be const initially
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mct,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mct,mc_opts.z_rnd_eps)
    else
      z_mc = z_mct
    end
    p_mc = copy(p)
    pref_mc = [SMCg{np,Float64}(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]

    #println("p_mc  :",p_mc)
    # setsup aff_mc, H_mc storage
    aff_mc = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    aff_mct = [SMCg{np,Float64}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    # h DAG graph

    H_mc = h(z_mc,p_mc)
    YH_mc = copy(H_mc)
    #println("isolation point #1")
    dH_mc = hj(aff_mc,p_mc)
    #println("isolation point #2")
    YdH_mc = copy(dH_mc)
    J_int = [dH_mc[i,j].Intv for i=1:nx,j=1:nx]
    Jm = J_int
    Y = (size(Jm) == 1) ? [1.0/Jm[1]] : inv(Jm)

    # stores some things for repeated us
    optc = [szero,sone]

    # Begins loop to generate parameters
    for k=1:mc_opts.kmax
      xa_mc,xA_mc,z_mc = Affine_Exp(param[k],p_mc,pref_mc,xa_mc,xA_mc,z_mc,exp_opt)
      if mc_opts.z_rnd_all == true
        z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
      elseif mc_opts.z_rnd == true
        z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
      end

      # sets up affine relaxations and preconditioning
      aff_mc = [SMCg(xA_mc[i].cc,xa_mc[i].cv,xA_mc[i].cc_grad,xa_mc[i].cv_grad,
                     @interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
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
      if mc_opts.hhj_rnd_all == true
        H_mc,dH_mc = Rnd_Out_H_All(H_mc,dH_mc,mc_opts.hhj_rnd_all_eps)
      elseif mc_opts.hhj_rnd == true
        H_mc,dH_mc = Rnd_Out_H_Intv(H_mc,dH_mc,mc_opts.hhj_rnd_eps)
      end
      # applies parametric iteration
      if (mc_opts.style == "NewtonGS")
        x_mc = MC_NewtonGS(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
      elseif (mc_opts.style == "KrawczykCW")
        x_mc = MC_KrawczykCW(z_mc,x_mc,dH_mc,H_mc,nx,np,optc)
      else
        error("Unsupported Style of Implicit Relaxation")
      end
    end
    return x_mc
end
