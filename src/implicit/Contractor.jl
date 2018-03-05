function MC_NewtonGS(z_mc,x_mc,YdH_mc,YH_mc,nx,np,optc)
  S1::SMCg = SMCg(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
  x_mc_int = copy(x_mc)
  for i=1:nx
    S1 = SMCg(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
    for j=1:nx
      #println("i: ", i)
      #println("j: ", j)
      #println("YdH_mc[i,j]: ", YdH_mc[i,j])
      #println("x_mc[j]-z_mc[j]: ", x_mc[j]-z_mc[j])
      #println("z_mc[j]: ", z_mc[j])
      #println("x_mc[j]: ", x_mc[j])
      if (i<j)
        S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
      elseif (j>i)
        S1 = S1 + YdH_mc[i,j]*(x_mc[j]-z_mc[j])
      end
      #println("S1: ",S1)
    end
    #println("z_mc[i]: ",z_mc[i])
    #println("YH_mc[i]: ",YH_mc[i])
    #println("YdH_mc[i,i]: ",YdH_mc[i,i])
    #println("(YH_mc[i]+S1)/YdH_mc[i,i]: ",(YH_mc[i]+S1)/YdH_mc[i,i])
    println("YdH_mc[i,i]: ", YdH_mc[i,i])
    x_mc[i] = z_mc[i] - (YH_mc[i]+S1)/YdH_mc[i,i]
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
  return x_mc
end

function MC_KrawczykCW(z_mc,x_mc,YdH_mc,YH_mc,nx,np,optc)
  S1::SMCg = SMCg(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
  x_mc_int = copy(x_mc)
  for i=1:nx
    S1 = SMCg(zero(x_mc[1].cc),zero(x_mc[1].cv),optc[1],optc[1],zero(x_mc[1].Intv),true,x_mc[1].IntvBox,x_mc[1].xref)
    for j=1:nx
      #println("i: ", i)
      #println("j: ", j)
      #println("YdH_mc[i,j]: ", YdH_mc[i,j])
      #println("x_mc[j]-z_mc[j]: ", x_mc[j]-z_mc[j])
      #println("z_mc[j]: ", z_mc[j])
      #println("x_mc[j]: ", x_mc[j])
      #println("x_mc[j].cnst: ", x_mc[j].cnst)
      #println("z_mc[j].cnst: ", z_mc[j].cnst)
      if (i<j)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      elseif (j>i)
        S1 = S1 - (YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      else
        S1 = S1 + (1.0-YdH_mc[i,j])*(x_mc[j]-z_mc[j])
      end
      #println("S1: ",S1)
    end
    #println("z_mc[i]: ",z_mc[i])
    #println("YH_mc[i]: ",YH_mc[i])
    x_mc[i] =  z_mc[i] - YH_mc[i] + S1
    x_mc[i] = Final_Cut(x_mc[i],x_mc_int[i])
  end
  return x_mc
end
