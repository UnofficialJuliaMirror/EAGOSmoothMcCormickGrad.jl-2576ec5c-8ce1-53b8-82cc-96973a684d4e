#=
abstract type abstractInterval <: Real end
struct Interval <: abstractInterval
  lo
  hi
end
∞ = Inf
Interval(x::S) where {S<:Integer}  = Interval(x,x)
Interval(x::S) where {S<:AbstractFloat}  = Interval(x,x)
Interval(x::Interval) = x
Interval(a::Tuple) = Interval(a...)
#Interval(x) = Interval(x,x)

diam(x::Interval) = x.hi- x.lo

emptyinterval() = Interval(Inf, -Inf)
isempty(x::Interval) = x.lo == Inf && x.hi == -Inf

entireinterval() = Interval(-Inf, Inf)
isentire(x::Interval) = x.lo == -Inf && x.hi == Inf

isunbounded(x::Interval) = x.lo == -Inf || x.hi == Inf

zero(x::Interval) = Interval(zero(x.lo),zero(x.lo))
one(x::Interval) = Interval(one(x.lo),one(x.lo))
Base.iszero(x::Interval) = iszero(x.lo) && iszero(x.hi)

real(x::Interval) = Interval(real(x.lo),real(x.hi))

promote_rule(::Type{Interval}, ::Type{S}) where {S<:AbstractFloat} =
    Interval
promote_rule(::Type{Interval}, ::Type{S}) where {S<:Integer} =
    Interval

function convert(::Type{Interval},x::S) where {S<:Integer}
          Interval(x,x)
end
function convert(::Type{Interval},x::S) where {S<:AbstractFloat}
          Interval(x,x)
end
function convert(::Type{Interval},x::Interval)
          x
end

function ==(a::Interval, b::Interval)
    isempty(a) && isempty(b) && return true
    a.lo == b.lo && a.hi == b.hi
end
!=(a::Interval, b::Interval) = !(a==b)

function islessprime(a, b)
    (isinf(a) || isinf(b)) && a==b && return true
    a < b
end

# Weakly less, \le, <=
function <=(a::Interval, b::Interval)
    isempty(a) && isempty(b) && return true
    (isempty(a) || isempty(b)) && return false
    (a.lo ≤ b.lo) && (a.hi ≤ b.hi)
end

# Strict less: <
function <(a::Interval, b::Interval)
    isempty(a) && isempty(b) && return true
    (isempty(a) || isempty(b)) && return false
    islessprime(a.lo, b.lo) && islessprime(a.hi, b.hi)
end

# precedes
function precedes(a::Interval, b::Interval)
    (isempty(a) || isempty(b)) && return true
    a.hi ≤ b.lo
end

# strictpreceds
function strictprecedes(a::Interval, b::Interval)
    (isempty(a) || isempty(b)) && return true
    # islessprime(a.hi, b.lo)
    a.hi < b.lo
end

∩(a::Interval,b::Interval) = Interval(max(a.lo,b.lo),min(a.hi,b.hi))
dist(a::Interval, b::Interval) = max(abs(a.lo-b.lo), abs(a.hi-b.hi))
eps(a::Interval) = max(eps(a.lo), eps(a.hi))

function mag(a::Interval)
    isempty(a) && return convert(eltype(a), NaN)
    max( abs(a.lo), abs(a.hi) )
end

function mig(a::Interval)
    isempty(a) && return convert(eltype(a), NaN)
    a.lo <= zero(a.lo) <= a.lo==hi && return zero(a.lo)
    min( abs(a.lo), abs(a.hi ))
end

min(a::Interval, b::Interval) = Interval(min(a.lo, b.lo), min(a.hi, b.hi))
max(a::Interval, b::Interval) = Interval(max(a.lo, b.lo), max(a.hi, b.hi))
hull(a::Interval, b::Interval) = Interval(min(a, b), max(a, b))

########## Defines nonvalidated interval operations
# Defines constants
pi_interval(T) = Interval(pi,pi)
half_pi(::Type{T}) where T = pi_interval(T) / 2
half_pi(x::T) where T<:AbstractFloat = half_pi(T)
two_pi(::Type{T})  where T = pi_interval(T) * 2

# Defines for on increasing univariant functions
inc_list1 = [:exp,:exp2,:exp10,:tanh,:sinh,:asinh,:atan]
inc_list2 = [:log,:log2,:log10,:sqrt]
dec_list = [:acos,:acosh]
mon_list = [:sqr,:cosh,:abs]
val1_list = [:sin,:cos]
val2_list = [:pow,:^,+,*,-,inv]
# defines nonvalidated interval arithmetic over increasing functions
for j in inc_list1
  eval(quote function ($j)(x::Interval)
              isempty(x) && return emptyinterval()
              Interval(($j)(x.lo),($j)(x.hi))
        end end)
end
for j in inc_list2
  eval(quote function ($j)(x::Interval)
              isempty(x) && return emptyinterval()
              domain = Interval(0, Inf)
              y = x ∩ domain
              Interval(($j)(y.lo),($j)(y.hi))
        end end)
end
function acosh(x::Interval)
  isempty(x) && return emptyinterval()
  domain = Interval(1, Inf)
  y = x ∩ domain
  Interval(acosh(y.lo),acosh(y.hi))
end
function atanh(x::Interval)
  isempty(x) && return emptyinterval()
  domain = Interval(-one(x.lo), one(x.lo))
  y = x ∩ domain
  Interval(atanh(y.lo),atanh(y.hi))
end
function asin(x::Interval)
  isempty(x) && return emptyinterval()
  domain = Interval(-one(x.lo), one(x.lo))
  y = x ∩ domain
  Interval(asin(y.lo),asin(y.hi))
end
function tan(x::Interval)
  isempty(x) && return emptyinterval()
  domain = Interval(-one(x.lo), one(x.lo))
  y = x ∩ domain
  Interval(tan(y.lo),tan(y.hi))
end
# defines nonvalidated interval arithmetic over decreasing functions
for j in dec_list
  eval(quote function $j(x::Interval)
              isempty(x) && return emptyinterval(x)
              domain = Interval(-one(x.lo), -one(x.lo))
              y = x ∩ domain
              Interval($j(y.hi),$j(y.lo))
            end end)
end
# defines nonvalided interval arithmetic functions with monotone derivatives & roots at zero
for j in mon_list
  eval( quote function ($j)(x::Interval)
          isempty(x) && return emptyinterval(x)
          xmin = (($j)(x.lo)<($j)(x.hi)) ? ($j)(x.lo) : ($j)(x.hi)
          xmax = (($j)(x.lo)<($j)(x.hi)) ? ($j)(x.hi) : ($j)(x.lo)
          return ((xmin<=zero(x.lo)) & (xmax>=zero(x.lo))) ? Interval(zero(x.lo),xmax) : Interval(xmin,xmax)
        end end )
end
# defines nonvalided interval arithmetic functions for step function
@inline function step(x::Interval)
           isempty(x) && return emptyinterval(x)
           xmin = ((x.lo)<zero(x.lo)) ? zero(x.lo) : one(x.lo)
           xmax = ((x.hi)>=zero(x.lo)) ? one(x.lo) : zero(x.lo)
           return Interval(xmin,xmax)
end
# defines nonvalided interval arithmetic functions for sign function
@inline function sign(x::Interval)
           isempty(x) && return emptyinterval(x)
           xmin = ((x.lo)<zero(x.lo)) ? -one(x.lo) : one(x.lo)
           xmax = ((x.hi)>=zero(x.lo)) ? one(x.lo) : -one(x.lo)
           return Interval(xmin,xmax)
end
function floor(a::Interval)
    isempty(a) && return emptyinterval(a)
    Interval(Base.floor(a.lo), Base.floor(a.hi))
end
function ceil(a::Interval)
    isempty(a) && return emptyinterval(a)
    Interval(Base.ceil(a.lo), Base.ceil(a.hi))
end
function trunc(a::Interval)
    isempty(a) && return emptyinterval(a)
    Interval(Base.trunc(a.lo), Base.trunc(a.hi))
end
function mid(a::Interval)

    isempty(a) && return NaN
    isentire(a) && return zero(a.lo)

    a.lo == -∞ && return nextfloat(a.lo)
    a.hi == +∞ && return prevfloat(a.hi)

    return 0.5 * (a.lo + a.hi)
end

# defines nonvalided interval arithmetic functions for *
@inline function (*)(x1::Interval, x2::Interval)
     (isempty(x1)|isempty(x2)) && return emptyinterval(x)
     a = x1.lo
     b = x1.hi
     c = x2.lo
     d = x2.hi
     ac = a*c
     ad = a*d
     bc = b*c
     bd = b*d
     return Interval(min(ac,ad,bc,bd), max(ac,ad,bc,bd))
end
# defines nonvalided interval arithmetic functions for /
@inline function (/)(x1::Interval, x2::Interval)
    (isempty(x1)|isempty(x2)) && return emptyinterval(x)
    a = x1.lo
    b = x1.hi
    c = x2.lo
    d = x2.hi
    ac = a/c
    ad = a/d
    bc = b/c
    bd = b/d
    return Interval(min(ac,ad,bc,bd), max(ac,ad,bc,bd))
end
# defines nonvalided interval arithmetic functions for inv
function inv(a::Interval)
  isempty(a) && return emptyinterval(a)
  if a.lo <= zero(a.lo) <= a.hi
    a.lo < zero(a.lo) == a.hi && return Interval(-Inf, inv(a.lo))
    a.lo == zero(a.lo) < a.hi && return Interval(inv(a.hi), Inf)
    a.lo < zero(a.lo) < a.hi && return entireinterval()
    a == zero(a.lo) && return emptyinterval()
  end
  return Interval(inv(a.hi), inv(a.lo))
end


+(x::Interval,y::Interval) = Interval(x.lo+y.lo,x.hi+y.hi)
-(x::Interval) = Interval(-x.hi,-x.lo)
for i in union(int_list, float_list)
  eval(quote
        +(x::Interval,y::$i) = Interval(x.lo+y,x.hi+y)
        +(y::$i,x::Interval) = Interval(x.lo+y,x.hi+y)
        -(x::Interval,y::$i) = x + (-y)
        -(x::$i,y::Interval) = x + (-y)
        *(y::Interval,x::$i) = x >= 0 ? Interval(x*y.lo,x*y.hi) : Interval(x*y.hi,x*y.lo)
        *(x::$i,y::Interval) = x >= 0 ? Interval(x*y.lo,x*y.hi) : Interval(x*y.hi,x*y.lo)
        /(x::Interval,y::$i) = x*(1/y)
        /(x::$i,y::Interval) = x*inv(y)
       end)
end

function find_quadrants(x)
    temp = x / half_pi(x)
    (Base.floor(temp.lo), Base.floor(temp.hi))
end
function cos(a::Interval)
    T = eltype(a.lo)
    isempty(a) && return a

    whole_range = Interval(-one(T), one(T))

    diam(a) > two_pi(T).lo && return whole_range

    lo_quadrant = minimum(find_quadrants(a.lo))
    hi_quadrant = maximum(find_quadrants(a.hi))

    if hi_quadrant - lo_quadrant > 4  # close to limits
        return Interval(-one(T), one(T))
    end

    lo_quadrant = mod(lo_quadrant, 4)
    hi_quadrant = mod(hi_quadrant, 4)

    # Different cases depending on the two quadrants:
    if lo_quadrant == hi_quadrant # Interval limits in the same quadrant
        a.hi - a.lo > pi_interval(T).lo && return whole_range
        lo = Interval(cos(a.lo), cos(a.lo))
        hi = Interval(cos(a.hi), cos(a.hi))
        return hull(lo, hi)

    elseif lo_quadrant == 2 && hi_quadrant==3
        return Interval(cos(a.lo), cos(a.hi))

    elseif lo_quadrant == 0 && hi_quadrant==1
        return Interval(cos(a.hi), cos(a.lo))

    elseif ( lo_quadrant == 2 || lo_quadrant==3 ) && ( hi_quadrant==0 || hi_quadrant==1 )
        return Interval(min(cos(a.lo), cos(a.hi)), 1)

    elseif ( lo_quadrant == 0 || lo_quadrant==1 ) && ( hi_quadrant==2 || hi_quadrant==3 )
        return Interval(-1, max(cos(a.lo), cos(a.hi)))

    else#if ( lo_quadrant == 3 && hi_quadrant==2 ) || ( lo_quadrant == 1 && hi_quadrant==0 )
        return whole_range
    end
end
function sin(a::Interval)
  cos(pi/2-a)
end

function pow(a::Interval,n::Integer)
  isempty(a) && return a
    n == 0 && return one(a)
    n == 1 && return a
    # n == 2 && return sqr(a)
    n < 0 && a == zero(a) && return emptyinterval(a)

    if isodd(n) # odd power
        isentire(a) && return a
        if n > 0
            a.lo == 0 && return Interval(0, a.hi^n)
            a.hi == 0 && return Interval(a.lo^n, 0)
            return Interval(a.lo^n, a.hi^n)
        else
            if a.lo ≥ 0
                a.lo == 0 && return Interval(a.hi^n, Inf)
                return Interval(a.hi^n, a.lo^n)

            elseif a.hi ≤ 0
                a.hi == 0 && return Interval(-Inf, a.lo^n)
                return Interval(a.hi^n, a.lo^n)
            else
                return entireinterval(a)
            end
        end

    else # even power
        if n > 0
            if a.lo ≥ 0
                return Interval(a.lo^n, a.hi^n)
            elseif a.hi ≤ 0
                return Interval(a.hi^n, a.lo^n)
            else
                return Interval(mig(a)^n, mag(a)^n)
            end

        else
            if a.lo ≥ 0
                return Interval(a.hi^n, a.lo^n)
            elseif a.hi ≤ 0
                return Interval(a.lo^n, a.hi^n)
            else
                return Interval(mag(a)^n, mig(a)^n)
            end
        end
    end
end

^(a::Interval,c::Integer) = pow(a,c)
=#
