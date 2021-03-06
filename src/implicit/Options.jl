"""
    mc_opts{T}

Storage type used for parameter options in implicit bounding routine.
* `lambda::T`: Affine weighting parameter
* `kmax::Int64`: Number of iterations run
* `style::String`: Type of contractor used
* `z_rnd::Bool`: Flag for rounding z on interval bounds
* `z_rnd_eps::T`: Amount for rounding z on interval bounds
* `z_rnd_all::Bool`: Flag for rounding z on all bounds
* `z_rnd_all_eps::T`: Amount for rounding z on all bounds
* `aff_rnd::Bool`: Flag for rounding affine relaxation on interval bounds
* `aff_rnd_eps::T`: Amount for rounding affine relaxation on interval bounds
* `aff_rnd_all::Bool`: Flag for rounding affine relaxation on all bounds
* `aff_rnd_all_eps::T`: Amount for rounding affine relaxation on all bounds
* `hhj_rnd::Bool`: Flag for rounding Jacobian on interval bounds
* `hhj_rnd_eps::T`: Amount for rounding Jacobian on interval bounds
* `hhj_rnd_all::Bool`: Flag for rounding Jacobian on all bounds
* `hhj_rnd_all_eps::T`: Amount for rounding Jacobian on all bounds
* `aff_correct_eps::T`: Affine correction tolerance
"""
type mc_opts{T<:AbstractFloat}
  lambda::T
  kmax::Int64
  style::String
  z_rnd::Bool
  z_rnd_eps::T
  z_rnd_all::Bool
  z_rnd_all_eps::T
  aff_rnd::Bool
  aff_rnd_eps::T
  aff_rnd_all::Bool
  aff_rnd_all_eps::T
  hhj_rnd::Bool
  hhj_rnd_eps::T
  hhj_rnd_all::Bool
  hhj_rnd_all_eps::T
  aff_correct_eps::T
end

"""
    mc_opts(T)

Initialization function for currently sets the weight `.lambda = 0.5`, the
contractor style `.style = KrawczykCW`, the number of iterations to `.kmax = 2`,
and the affine correction tolerance as `.aff_correct_eps = 1E-12`. Other rounding
options disabled.
"""
mc_opts(T) = mc_opts{T}(0.5,2,"KrawczykCW",
                        false,0.0,false,0.0,
                        false,0.0,false,0.0,
                        false,0.0,false,0.0,1E-12)

mc_opts() = mc_opts(Float64)
"""
    set_default!(x::mc_opts)

Restores `x` to default values. Sets the weight `.lambda = 0.5`, the
contractor style `.style = KrawczykCW`, the number of iterations to `.kmax = 2`,
and the affine correction tolerance as `.aff_correct_eps = 1E-12`. Other rounding
options disabled.
"""
function set_default!(x::mc_opts{T}) where {T<:AbstractFloat}
  x.lambda = 0.5
  x.kmax = 2
  x.style = "KrawczykCW"
  x.z_rnd = false
  x.z_rnd_eps = 0.0
  x.z_rnd_all= false
  x.z_rnd_all_eps = 0.0
  x.aff_rnd= false
  x.aff_rnd_eps = 0.0
  x.aff_rnd_all= false
  x.aff_rnd_all_eps = 0.0
  x.hhj_rnd= false
  x.hhj_rnd_eps = 0.0
  x.hhj_rnd_all= false
  x.hhj_rnd_all_eps = 0.0
  x.aff_correct_eps = 1E-12
end
