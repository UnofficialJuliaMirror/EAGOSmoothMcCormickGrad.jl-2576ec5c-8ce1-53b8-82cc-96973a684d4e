########### Defines root-finding functions
# newton's method root-finding algorithm
@inline function newton(x0,xL,xU,f,df,rinp,iinp)
  dfk = zero(x0)
  xk = zero(x0)
  fk = zero(x0)

  xk = max(xL,min(x0,xU))
  fk= f(xk,rinp,iinp)

  for i=1:MC_param.env_max_int
    dfk = df(xk,rinp,iinp)
    if (abs(fk)<MC_param.env_tol)
      return xk
    end
    if (dfk == zero(x0))
      error("NEWTON EXCEPTION")
    elseif (xk==xL && fk/dfk>zero(x0))
      return xk
    elseif (xk==xU && fk/dfk<zero(x0))
      return xk
    end
    xk = max(xL,min(xU,xk-fk/dfk))
    fk = f(xk,rinp,iinp)
  end
  error("NEWTON EXCEPTION")
end
# secant method root-finding algorithm
@inline function secant(x0,x1,xL,xU,f,rinp,iinp)
  xkm = max(xL,min(xU,x0))
  xk = max(xL,min(xU,x1))
  fkm = f(xkm,rinp,iinp)

  for i=1:MC_param.env_max_int
    fk = f(xk,rinp,iinp)
    Bk = (fk-fkm)/(xk-xkm)
    if (abs(fk)<MC_param.env_tol)
      return xk
    end
    if (Bk == zero(x0))
      error("SECANT EXCEPTION")
    elseif ((xk==xL)&(fk/Bk>zero(x0)))
      return xk
    elseif ((xk==xU)&(fk/Bkzero(x0)))
      return xk
    end
    xkm = xk
    fkm = fk
    xk = max(xL,min(xU,xk-fk/Bk))
  end
  error("SECANT EXCEPTION")
end
# golden section root-finding algorithm
@inline function golden_section(xL,xU,f,rinp,iinp)
  fL = f(xL,rinp,iinp)
  fU = f(xU,rinp,iinp)

  if (fL*fU > 0.0)
    error("GOLDEN EXCEPTION")
  end
  xm = xU-(2.0-golden)*(xU-xL)
  fm = f(xm,rinp,iinp)
  return golden_section_it(1,xL,fL,xm,fm,xU,fU,f,rinp,iinp)
end
@inline function golden_section_it(init::Integer,a,fa,b,fb,c,fc,f,rinp,iinp)
  b_t_x = (c-b > b-a)
  if (b_t_x)
    x = b + (2.0-golden)*(c-b)
  else
    x = b - (2.0-golden)*(b-a)
  end
  itr = init
  if (abs(c-a)<MC_param.env_tol*(abs(b)+abs(x))||(itr>MC_param.env_max_int))
    return (c+a)/2.0
  end
  itr += 1
  fx = f(x,rinp,iinp)
  if (b_t_x)
    if (fa*fx<zero(fa))
      golden_section_it(itr,a,fa,b,fb,x,fx,f,rinp,iinp)
    else
      golden_section_it(itr,b,fb,x,fx,c,fc,f,rinp,iinp)
    end
  else
    if (fa*fb<(fa))
      golden_section_it(itr,a,fa,x,fx,b,fb,f,rinp,iinp)
    else
      golden_section_it(itr,x,fx,b,fb,c,fc,f,rinp,iinp)
    end
  end
end
