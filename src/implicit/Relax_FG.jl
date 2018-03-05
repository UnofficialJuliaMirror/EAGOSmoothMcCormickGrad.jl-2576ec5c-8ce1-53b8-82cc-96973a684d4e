function impRelax_f(f,h,hj,X,P,p,pmid,mc_opt,param)
  nx::Int64 = length(X)
  np::Int64 = length(P)
  szero = @SVector zeros(np)
  sone = @SVector ones(np)
  #exp_opt = [nx,np,lambda]
  xL = [X[i].lo for i=1:nx]
  xU = [X[i].hi for i=1:nx]
  xa,xA,z = copy(param[1]),copy(param[1]),copy(param[1])
  XP = copy(X)
  x = [SMCg{np,Float64}(xU[i],xL[i],szero,szero,@interval(xL[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
  p_mc = [SMCg{np,Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC = MC_impRelax(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc)
end

function impRelax_fg(f,g,h,hj,X,P,p,pmid,mc_opt,param)
  nx::Int64 = length(X)
  np::Int64 = length(P)
  szero = @SVector zeros(np)
  sone = @SVector ones(np)
  #exp_opt = [nx,np,lambda]
  xL = [X[i].lo for i=1:nx]
  xU = [X[i].hi for i=1:nx]
  xa,xA,z = copy(param[1]),copy(param[1]),copy(param[1])
  XP = copy(X)
  x = [SMCg{np,Float64}(xU[i],xL[i],szero,szero,@interval(xL[i],xU[i]),true,[∅],[1.0]) for i=1:nx]
  p_mc = [SMCg{np,Float64}(p[i],p[i],sone,sone,@interval(P[i].lo,P[i].hi),false,[∅],[1.0]) for i=1:np]
  xpMC = MC_impRelax(h,hj,p_mc,pmid,X,P,mc_opt,param)
  return f(xpMC,p_mc),g(xpMC,p_mc)
end
