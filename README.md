# EAGOSmoothMcCormickGrad.jl
A (differentiable) McCormick Relaxation Library w/Embedded (sub)gradient

[![Build Status](https://travis-ci.org/MatthewStuber/EAGOSmoothMcCormickGrad.jl.svg?branch=master)](https://travis-ci.org/MatthewStuber/EAGOSmoothMcCormickGrad.jl)
[![Coverage Status](https://coveralls.io/repos/github/MatthewStuber/EAGOSmoothMcCormickGrad.jl/badge.svg?branch=master)](https://coveralls.io/github/MatthewStuber/EAGOSmoothMcCormickGrad.jl?branch=master)
[![codecov.io](http://codecov.io/github/MatthewStuber/EAGOSmoothMcCormickGrad.jl/coverage.svg?branch=master)](http://codecov.io/github/MatthewStuber/EAGOSmoothMcCormickGrad.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://MatthewStuber.github.io/EAGO.jl/latest)

## Authors

[Matthew Wilhelm](https://psor.uconn.edu/our-team/), Department of Chemical and Biomolecular Engineering,  University of Connecticut (UCONN)

## Installation

```julia
julia> Pkg.add("EAGOSmoothMcCormickGrad.jl")
```

## Capabilities

**EAGOSmoothMcCormickGrad.jl** provides a library of smooth McCormick objects `SMCg{N,T}` with associate rules
for propagating the relaxations via operator-overloading. A McCormick relaxation framework allows for the
easy computation of convex and concave relaxations of functions without defining auxiliary variables. Techniques
using this approach have been utilized in global optimization to bound implicit function such as those described
by differential equations. The user has the option to specify whether differentiable, multivariant, and/or affine
subgradient tightening variants.

The `SMCg{N,T}` has the following fields
- `cc::T`: the concave relaxation
- `cv::T`: the convex relaxation
- `cc_grad::SVector{N,T}`: the sub-gradient of cc if nonsmooth McCormick relaxations are used and
                           the gradient if differentiable McCormick relaxations are specified
- `cv_grad::SVector{N,T}`: the sub-gradient of cv if nonsmooth McCormick relaxations are used
- `Intv::Interval{T}`: interval bounds for the object
- `cnst::Bool`: used to denote constant
- `IntvBox::Vector{Interval{T}}`: domain of initial problem used for subgradient refinement.
- `xref::Vector{T}`: reference point defined for affine interval tightening procedure

The routine are used extensively in the [`EAGO.jl`](https://github.com/MatthewStuber/EAGO.jl) solver.
Please see the example files for usage cases.

## Future Work

- We're always looking to added new relaxations.
- Additionally, we'll be updating this to provide better support for control-flow logic.
- Lastly, we're developing a lighter weight package that better integrates with
  Julia's automatic differentiation ecosystem.

## Related Packages
- [**EAGO.jl**](https://github.com/MatthewStuber/EAGO.jl): A package containing global and robust solvers based mainly on McCormick relaxations.
This package supports a JuMP and MathProgBase interface.
- [**MC++**](https://omega-icl.github.io/mcpp/): A mature McCormick relaxation package in C++ that also includes McCormick-Taylor, Chebyshev
Polyhedral and Ellipsoidal arithmetics.

## References
- Chachuat, B.: MC++: a toolkit for bounding factorable functions, v1.0. Retrieved 2 July 2014 https://
projects.coin-or.org/MCpp (2014)
- A. Mitsos, B. Chachuat, and P. I. Barton. McCormick-based relaxations of algorithms.
SIAM Journal on Optimization, 20(2):573–601, 2009.
- G. P. McCormick. Computability of global solutions to factorable nonconvex programs:
Part I-Convex underestimating problems. Mathematical Programming, 10:147–175, 1976.
- G. P. McCormick. Nonlinear programming: Theory, Algorithms, and Applications. Wi-
ley, New York, 1983.
- J. K. Scott, M. D. Stuber, and P. I. Barton. Generalized McCormick relaxations. Journal
of Global Optimization, 51(4):569–606, 2011.
- Stuber, M.D., Scott, J.K., Barton, P.I.: Convex and concave relaxations of implicit functions. Optim.
Methods Softw. 30(3), 424–460 (2015)
- A. Tsoukalas and A. Mitsos. Multivariate McCormick Relaxations. Journal of Global
Optimization, 59:633–662, 2014.
