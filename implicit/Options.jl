type mc_opts
  lambda
  kmax
  style
  z_rnd
  z_rnd_eps
  z_rnd_all
  z_rnd_all_eps
  aff_rnd
  aff_rnd_eps
  aff_rnd_all
  aff_rnd_all_eps
  hhj_rnd
  hhj_rnd_eps
  hhj_rnd_all
  hhj_rnd_all_eps
  aff_correct_eps
end
mc_opts() = mc_opts(0.5,2,"KrawczykCW",
                    false,0.0,false,0.0,
                    false,0.0,false,0.0,
                    false,0.0,false,0.0,
                    1E-12)

function set_default!(x::mc_opts)
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
