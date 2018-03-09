"""
    GenExpansionParams(h::Function, hj::Function, X::Vector{Interval{T}},
                       P::Vector{Interval{T}},pmid::Vector{T},mc_opts::mc_opts)

Generates relaxation of state variable at pmid in P. Inputs are:
* `h::Function`: h(z,p) which implicit defines state function
* `hj::Function`: hj(z,p) w.r.t z which implicit defines state function
* `P::Vector{Interval}`: State variable bounds
* `X::Vector{Interval{T}}`: Decision variable bounds
* `pmid::Vector{T}`: Point at which to generate relaxation
* `mc_opts`: Options for generating implicit function relaxation
Outs the following the tuple `sto_out`:
* `sto_out::Vector{SMCg{np,T}}` - McCormick parameter object stack
--------------------------------------------------------------------------------
"""
function GenExpansionParams(h::Function, hj::Function,
                            X::Vector{Interval{T}},
                            P::Vector{Interval{T}},
                            pmid::Vector{T},mc_opts::mc_opts{T}) where {T<:AbstractFloat}

  nx::Int64 = length(X)
  np::Int64 = length(P)
  szero::SVector{np,T} = @SVector zeros(np)
  sone::SVector{np,T} = @SVector ones(np)
  exp_opt::Vector{Any} = Any[nx,np,mc_opts.lambda]
  sto_out::Vector{Vector{SMCg{np,T}}} = [[zero(SMCg{np,T}) for i=1:nx] for j=1:(mc_opts.kmax+1)]

  x_mc::Vector{SMCg{np,T}} = [SMCg{np,T}(X[i].hi,X[i].lo,szero,szero,@interval(X[i].lo,X[i].hi),false,[∅],[1.0]) for i=1:nx]
  xa_mc::Vector{SMCg{np,T}} = [SMCg{np,T}(X[i].lo,X[i].lo,szero,szero,@interval(X[i].lo,X[i].lo),false,[∅],[1.0]) for i=1:nx]
  xA_mc::Vector{SMCg{np,T}} = [SMCg{np,T}(X[i].hi,X[i].hi,szero,szero,@interval(X[i].hi,X[i].hi),false,[∅],[1.0]) for i=1:nx]
  if mc_opts.z_rnd_all == true
    z_mc::Vector{SMCg{np,T}} = Rnd_Out_Z_All(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_all_eps)
  elseif mc_opts.z_rnd == true
    z_mc = Rnd_Out_Z_Intv(mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc,mc_opts.z_rnd_eps)
  else
    z_mc = mc_opts.lambda*xa_mc+(1.0-mc_opts.lambda)*xA_mc
  end

  p_mc::Vector{SMCg{np,T}} = [SMCg(pmid[i],pmid[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  pref_mc::Vector{SMCg{np,T}} = copy(p_mc)
  aff_mc::Vector{SMCg{np,T}} = [SMCg{np,T}(xA_mc[i].cc,xa_mc[i].cv,szero,szero,@interval(xa_mc[i].cv,xA_mc[i].cc),false,[∅],[1.0]) for i=1:nx]

  H_mc::Vector{SMCg{np,T}} = h(z_mc,p_mc)
  dH_mc::Union{Vector{SMCg{np,T}},Array{SMCg{np,T},2}} = hj(aff_mc,p_mc)
  Y::Union{Vector{T},Array{T,2}} = (nx == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))
  # stores some things for repeated us
  optc = Any[szero,sone]

  # Begins loop to generate parameters
  sto_out[1] = copy(x_mc)

  for k=1:mc_opts.kmax
    aff_mc = [SMCg{np,Float64}(x_mc[i].cc,x_mc[i].cv,szero,szero,
                   @interval(x_mc[i].cv,x_mc[i].cc),false,[∅],[1.0]) for i=1:nx]
    if mc_opts.aff_rnd_all == true
      aff_mc = Rnd_Out_Z_All(aff_mc,mc_opts.aff_rnd_all_eps)
    elseif mc_opts.aff_rnd == true
      aff_mc = Rnd_Out_Z_Intv(aff_mc,mc_opts.aff_rnd_eps)
    end

    if mc_opts.hhj_rnd_all == true
      H_mc,dH_mc = Rnd_Out_H_All(h(z_mc,p_mc),hj(aff_mc,p_mc),mc_opts.hhj_rnd_all_eps)
    elseif mc_opts.hhj_rnd == true
      H_mc,dH_mc = Rnd_Out_H_Intv(h(z_mc,p_mc),hj(aff_mc,p_mc),mc_opts.hhj_rnd_eps)
    else
      H_mc,dH_mc = h(z_mc,p_mc),hj(aff_mc,p_mc)
    end
    Y = (nx == 1) ? [one(T)/mid(dH_mc[1].Intv)] : mid.(Intv.(dH_mc))

    # applies preconditioner
    Precondition!(H_mc,dH_mc,Y,nx)

    # applies parametric iteration
    if (mc_opts.style == "NewtonGS")
      MC_NewtonGS!(z_mc,x_mc,dH_mc,H_mc,nx)
    elseif (mc_opts.style == "KrawczykCW")
      MC_KrawczykCW!(z_mc,x_mc,dH_mc,H_mc,nx)
    else
      error("Unsupported Style of Implicit Relaxation")
    end

    # update affine relaxations & correct
    Affine_Exp!(x_mc,p_mc,p_mc,xa_mc,xA_mc,z_mc,exp_opt)
    Correct_Exp!(z_mc,x_mc,X,nx,np,mc_opts.aff_correct_eps)
    if mc_opts.z_rnd_all == true
      z_mc = Rnd_Out_Z_All(z_mc,mc_opts.z_rnd_all_eps)
    elseif mc_opts.z_rnd == true
      z_mc = Rnd_Out_Z_Intv(z_mc,mc_opts.z_rnd_eps)
    end

    # store relaxation
    sto_out[k+1] = copy(x_mc)
    end
  #return sto_out,sto_z,sto_x
  return sto_out
end
