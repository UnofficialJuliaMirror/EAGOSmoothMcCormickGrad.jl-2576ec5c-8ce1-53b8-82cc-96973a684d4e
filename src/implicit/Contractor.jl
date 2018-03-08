"""
    MC_NewtonGS(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                         YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                         nx::Int64,np::Int64,optc::Vector{Any})

Performs a Newton Gauss-Seidel bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MC_NewtonGS(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                       YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                       nx::Int64,np::Int64,optc::Vector{Any}) where {N,T}
  S1::SMCg{N,T} = SMCg{N,T}(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
  x_mc_int = copy(x_mc)
  for i=1:nx
    S1 = SMCg(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
    for j=1:nx
      if (i<j)
        S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
      elseif (j>i)
        S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
      end
    end
    x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
  return x_mc
end

"""
    MC_KrawczykCW(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                           YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                           nx::Int64,np::Int64,optc::Vector{Any})

Performs a componentwise Krawczyk bounding iteration using McCormick relaxations
where `x_mc` is the initial state variable, `z_mc` are relaxations of the affine
function in `x_mc`, `YdH_mc` is the preconditioned jacobian, `YH_mc` is the
equality constraints jacobian, `nx` size of the state variables, `np` is the
size of the decision space, and `optc` is a vector with stored values.
"""
function MC_KrawczykCW(z_mc::Vector{SMCg{N,T}},x_mc::Vector{SMCg{N,T}},
                       YdH_mc::Array{SMCg{N,T},2},YH_mc::Array{SMCg{N,T},2},
                       nx::Int64,np::Int64,optc::Vector{Any}) where {N,T}
  S1::SMCg{N,T} = SMCg{N,T}(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
  x_mc_int::Vector{SMCg{N,T}} = copy(x_mc)
  for i=1:nx
    S1 = SMCg{N,T}(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
    for j=1:nx
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j>i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (1.0-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
    end
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
  return x_mc
end
