"""
--------------------------------------------------------------------------------
Function: default_options
--------------------------------------------------------------------------------
Description:
Sets McCormickParameters to default.
env_max_int = 100 # number of iterations used for envelope calculations
env_tol = 1E-10 # tolerance used for envelope calculations
mu = 0 # nonsmooth McCormick relaxations used
valid_check = true # relaxtions and intervals are checked for validity
subgrad_refine = false # don't use interval refinement by subgradient propagation
multivar_refine = false # don't use multivariant relaxations
mv_tol = 1E-15          # tolerance for multivariant relaxations.
outer_rnding = false    # outer rounding of interval and relaxation disabled
outer_param = 0.0       # amount of outer rounding
--------------------------------------------------------------------------------
Inputs: None.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function default_options()
  MC_param = McCormickParamters()
  println("Options set to default values")
end

"""
--------------------------------------------------------------------------------
Function: set_iterations
--------------------------------------------------------------------------------
Description:
Sets number of iterations.
--------------------------------------------------------------------------------
Inputs:
val - Number of iterations.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_iterations(val)
  MC_param.env_max_int = val
  println("Maximum number of iteration for root-finding algorithms used to determing convex/concave envelopes set to $(MC_param.env_max_int)")
end

"""
--------------------------------------------------------------------------------
Function: set_tolerance
--------------------------------------------------------------------------------
Description:
Sets tolerance for envelope calculations.
--------------------------------------------------------------------------------
Inputs:
val - tolerance for envelope calculations.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_tolerance(val)
  MC_param.env_tol = val
end

"""
--------------------------------------------------------------------------------
Function: set_diff_relax
--------------------------------------------------------------------------------
Description:
Sets differentiability of McCormick relaxations.
--------------------------------------------------------------------------------
Inputs:
val - differentiability of McCormick relaxations.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_diff_relax(val::Integer)
  diff_relax = val>0
  if (diff_relax>0)
    MC_param.mu = val+1
  else
    MC_param.mu = 0
  end
end

function set_validated(val)
  MC_param.valid_check = valid_check
end

"""
--------------------------------------------------------------------------------
Function: set_valid_check
--------------------------------------------------------------------------------
Description:
Sets flag for validity check used.
--------------------------------------------------------------------------------
Inputs:
val - Bool.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_valid_check(val)
  MC_param.valid_check = val
end

"""
--------------------------------------------------------------------------------
Function: set_subgrad_refine
--------------------------------------------------------------------------------
Description:
Sets flag for using subgradient refinement.
--------------------------------------------------------------------------------
Inputs:
val - Bool.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_subgrad_refine(val)
  MC_param.subgrad_refine = val
end

"""
--------------------------------------------------------------------------------
Function: set_multivar_refine
--------------------------------------------------------------------------------
Description:
Sets flag for use of multivariant relaxations.
--------------------------------------------------------------------------------
Inputs:
bool - flag for whether multivariant relaxations will be used
tol - tolerance for multivariant relaxations.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_multivar_refine(bool,tol)
  MC_param.multivar_refine = bool
  MC_param.mv_tol = tol
end

"""
--------------------------------------------------------------------------------
Function: set_outer_rnd
--------------------------------------------------------------------------------
Description:
Sets outer rounding options.
--------------------------------------------------------------------------------
Inputs:
bool - flag for whether outer rounding will be used
tol - amount to outer round.
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function set_outer_rnd(bool,tol)
  MC_param.outer_rnding = bool
  MC_param.outer_param = tol
end

"""
--------------------------------------------------------------------------------
Function: show_options
--------------------------------------------------------------------------------
Description:
Prints current options to console.
--------------------------------------------------------------------------------
Inputs: None
--------------------------------------------------------------------------------
Returns: None.
--------------------------------------------------------------------------------
"""
function show_options()
end
