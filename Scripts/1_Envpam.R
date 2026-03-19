# Packages
# install.packages("devtools")
# install_github("macroecology/letsR")
library(devtools)
library(letsR)

# Data load
data("prec")
data("temp")
prec <- unwrap(prec)
temp <- unwrap(temp)
envs <- lets.addvar(PAM_crop, c(prec, temp), onlyvar = TRUE)

colnames(envs) <- c("Preciptation", "Temperature")
wrld_simpl <- get(utils::data("wrld_simpl", package = "letsR"))
PAM_crop42 <- lets.pamcrop(PAM_crop, vect(wrld_simpl))


res <- lets.envpam(
  PAM_crop42,
  envs,
  n_bins = 30,
  remove.cells = TRUE,
  remove.sp = TRUE,
  count = FALSE
)

lets.plot.envpam(
  res,
  species = NULL,
  cell_id_env = NULL,
  cell_id_geo = NULL,
  geo_plot = TRUE,
  env_plot = TRUE,
  world = TRUE,
  rast_return = FALSE,
  col_rich = NULL,
)

str(res)

out <- lets.envcells(res, perc = 0.2)


str(envs)
any(is.na(envs))
sum(is.na(envs))
colSums(is.na(envs))
