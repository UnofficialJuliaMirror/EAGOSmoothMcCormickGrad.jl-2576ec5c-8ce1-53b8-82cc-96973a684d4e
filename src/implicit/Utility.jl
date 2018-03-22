"""
    Smooth_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable fashion.
"""
function Smooth_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T}) where {N,T<:AbstractFloat}
  t_cv::SMCg{N,T} = max(x_mc,x_mc_int.Intv.lo)
  t_cc::SMCg{N,T} = min(x_mc,x_mc_int.Intv.hi)
  return SMCg{N,T}(t_cc.cc,t_cv.cv,t_cc.cc_grad,t_cv.cv_grad,
                   (x_mc.Intv ∩ x_mc_int.Intv),(t_cv.cnst && t_cc.cnst),
                   x_mc.IntvBox,x_mc.xref)
end

"""
    Final_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T})

An operator that cuts the `x_mc` object using the `x_mc_int bounds` in a
differentiable or nonsmooth fashion as specified by the `MC_param.mu flag`.
"""
function Final_Cut(x_mc::SMCg{N,T},x_mc_int::SMCg{N,T}) where {N,T<:AbstractFloat}
  if (MC_param.mu < 1)
    Intv::Interval = x_mc.Intv ∩ x_mc_int.Intv
    if (x_mc.cc <= x_mc_int.cc)
      cc::T = x_mc.cc
      cc_grad::SVector{N,T} = x_mc.cc_grad
    else
      cc = x_mc_int.cc
      cc_grad = x_mc_int.cc_grad
    end
    if (x_mc.cv >= x_mc_int.cv)
      cv::T = x_mc.cv
      cv_grad::SVector{N,T} = x_mc.cv_grad
    else
      cv = x_mc_int.cv
      cv_grad = x_mc_int.cv_grad
    end
    x_mc::SMCg{N,T} = SMCg{N,T}(cc,cv,cc_grad,cv_grad,(x_mc.Intv ∩ x_mc_int.Intv),x_mc.cnst,x_mc.IntvBox,x_mc.xref)
  else
    x_mc = Smooth_Cut(x_mc,x_mc_int)
  end
  return x_mc
end

"""
    Rnd_Out_Z_Intv(z_mct::SMCg{N,T},epsvi::Float64)

Rounds the interval of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_Intv(z_mct::Vector{SMCg{N,T}},epsvi::S) where {N,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  return [SMCg{N,T}(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             Interval{T}(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
end

"""
    Rnd_Out_Z_All(z_mct::Vector{SMCg{N,T}},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` vector elements out by `epsvi`.
"""
function Rnd_Out_Z_All(z_mct::Vector{SMCg{N,T}},epsvi::S) where {N,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  return [SMCg{N,T}(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             Interval{T}(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
end

"""
    Rnd_Out_H_Intv(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S)

Rounds the interval bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_All(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S) where {N,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  temp1::Vector{SMCg{N,T}} = [SMCg{N,T}(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
                                        z_mct[i].cc_grad, z_mct[i].cv_grad,
                                        Interval{T}(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
                                        z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
  temp2::Array{SMCg{N,T},2} = [SMCg{N,T}(Y_mct[i,j].cc+epsv,Y_mct[i,j].cv-epsv,
                                        Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
                                        Interval{T}(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
                                        Y_mct[i,j].cnst, Y_mct[i,j].IntvBox,Y_mct[i,j].xref) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

"""
    Rnd_Out_H_All(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S)

Rounds the interval and relaxation bounds of the `z_mct` and `Y_mct` elements out by `epsvi`.
"""
function Rnd_Out_H_Intv(z_mct::Vector{SMCg{N,T}},Y_mct::Array{SMCg{N,T},2},epsvi::S) where {N,S<:AbstractFloat,T<:AbstractFloat}
  epsv::T = convert(T,epsvi)
  temp1::Vector{SMCg{N,T}} = [SMCg(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             @interval(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
  temp2::Array{SMCg{N,T},2} = [SMCg(Y_mct[i,j].cc,Y_mct[i,j].cv,
             Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
             @interval(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
             Y_mct[i,j].cnst, Y_mct[i,j].IntvBox,Y_mct[i,j].xref) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

"""
    Precondition(hm::Vector{SMCg{N,T}},hJm::Union{Vector{SMCg{N,T}},Array{SMCg{N,T},2}},
                 Y::Union{Vector{T},Array{T,2}},nx::Int64)

Preconditions `hm` and `hJm` by `Y` in place where all dimensions are `nx`.
"""
function Precondition!(hm::Vector{SMCg{N,T}},hJm::Union{Vector{SMCg{N,T}},Array{SMCg{N,T},2}},
                      Y::Union{Vector{T},Array{T,2}},nx::Int64) where {N,T<:AbstractFloat}
  S1::SMCg{N,T},S2::SMCg{N,T} = zero(SMCg{N,T}),zero(SMCg{N,T})
  for i=1:nx
    S2 = zero(SMCg{N,T})
    for j=1:nx
      S1 = zero(SMCg{N,T})
      for k=1:nx
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end
