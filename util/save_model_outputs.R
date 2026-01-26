# save full model outputs as individual files, one for each model

require(brms)
# fitted cover error models
load("output/cover_error.rda")
# letters indicate wether model was fit for same surveyor ("s") or different suveyors ("d")
# fitted to all species ("0"), or to graminoids vs. non-graminoids ("g")
mod_c0s = update(mod_c0s, save_pars = save_pars(all = FALSE))
mod_c0d = update(mod_c0d, save_pars = save_pars(all = FALSE))
mod_cgs = update(mod_cgs, save_pars = save_pars(all = FALSE))
mod_cgd = update(mod_cgd, save_pars = save_pars(all = FALSE))

# fitted presence/absence models
load("output/id_error.rda")
mod_ps = update(mod_ps, save_pars = save_pars(all = FALSE))
mod_pd = update(mod_pd, save_pars = save_pars(all = FALSE))
mod_pgs = update(mod_pgs, save_pars = save_pars(all = FALSE))
mod_pgd = update(mod_pgd, save_pars = save_pars(all = FALSE))

save(list = c("mod_c0s","mod_c0d","mod_cgs","mod_cgd",
              "mod_ps","mod_pd","mod_pgs","mod_pgd"),
     file = "output/brms_coverid_models_small.rda")
