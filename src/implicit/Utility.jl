"""
"""
function Smooth_Cut(x_mc::SMCg,x_mc_int::SMCg)
  t_cv = max(x_mc,x_mc_int.Intv.lo)
  t_cc = min(x_mc,x_mc_int.Intv.hi)
  cnst = t_cv.cnst && t_cc.cnst
  Intv::Interval = x_mc.Intv ∩ x_mc_int.Intv
  x_mc = SMCg(t_cc.cc,t_cv.cv,t_cc.cc_grad,t_cv.cv_grad,Intv,cnst,x_mc.IntvBox,x_mc.xref)
  return x_mc
end

function Final_Cut(x_mc::SMCg,x_mc_int::SMCg)
  if (MC_param.mu < 1)
    Intv::Interval = x_mc.Intv ∩ x_mc_int.Intv
    println("Final Cut #1 Intv: ", Intv)
    if (x_mc.cc <= x_mc_int.cc)
      cc = x_mc.cc
      cc_grad = x_mc.cc_grad
    else
      cc = x_mc_int.cc
      cc_grad = x_mc_int.cc_grad
    end
    if (x_mc.cv >= x_mc_int.cv)
      cv = x_mc.cv
      cv_grad = x_mc.cv_grad
    else
      cv = x_mc_int.cv
      cv_grad = x_mc_int.cv_grad
    end
    println("Final Cut #2 cc: ", cc)
    println("Final Cut #2 cv: ", cv)
    println("Final Cut #2 cc_grad: ", cc_grad)
    println("Final Cut #2 cv_grad: ", cv_grad)
    x_mc = SMCg(cc,cv,cc_grad,cv_grad,Intv,x_mc.cnst,x_mc.IntvBox,x_mc.xref)
  else
    println("Smooth_Cut Branch")
    x_mc = Smooth_Cut(x_mc,x_mc_int)
  end
  return x_mc
end

function Rnd_Out_Z_Intv(z_mct,epsv)
  return [SMCg(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             @interval(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
end

function Rnd_Out_Z_All(z_mct,epsv)
  return [SMCg(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             @interval(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
end

function Rnd_Out_H_Intv(z_mct,Y_mct,epsv)
  temp1 = [SMCg(z_mct[i].cc+epsv,z_mct[i].cv-epsv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             @interval(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
  temp2 = [SMCg(Y_mct[i,j].cc+epsv,Y_mct[i,j].cv-epsv,
             Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
             @interval(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
             Y_mct[i,j].cnst, Y_mct[i,j].IntvBox,Y_mct[i,j].xref) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

function Rnd_Out_H_All(z_mct,Y_mct,epsv)
  temp1 = [SMCg(z_mct[i].cc,z_mct[i].cv,
             z_mct[i].cc_grad, z_mct[i].cv_grad,
             @interval(z_mct[i].Intv.lo-epsv, z_mct[i].Intv.hi+epsv),
             z_mct[i].cnst, z_mct[i].IntvBox,z_mct[i].xref) for i=1:length(z_mct)]
  temp2 = [SMCg(Y_mct[i,j].cc,Y_mct[i,j].cv,
             Y_mct[i,j].cc_grad, Y_mct[i,j].cv_grad,
             @interval(Y_mct[i,j].Intv.lo-epsv, Y_mct[i,j].Intv.hi+epsv),
             Y_mct[i,j].cnst, Y_mct[i,j].IntvBox,Y_mct[i,j].xref) for i=1:length(z_mct), j=1:length(z_mct)]
  return temp1,temp2
end

"""
    Precondition(hm,hJm,Y,nx)

Preconditions `hm` and `hJm` by `Y` where all dimensions are `nx`. Returns the
tuple `(Y*hm,Y*hJm)`.
"""
function Precondition(hm,hJm,Y,nx)
  S1,S2 = 0.0,0.0
  for i=1:nx
    S2 = 0.0
    for j=1:nx
      S1 = 0.0
      for k=1:nx
        temp1 = Y[i,k]
        temp2 = hJm[k,j]
        S1 = S1 + Y[i,k]*hJm[k,j]
      end
      hJm[i,j] = S1
      S2 += Y[i,j]*hm[j]
    end
    hm[i] = S2
  end
  return hm,hJm
end
