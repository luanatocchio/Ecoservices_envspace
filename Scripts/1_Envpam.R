# Code to map species across spaces

# Packages ----------------------------------------------------------------
## install.packages("devtools")
## library(devtools)
## install_github("macroecology/letsR")
library(letsR)
library(terra)
library(geodata)
library(tidyverse)

# Load Data ---------------------------------------------------------------
## Env data
bio_data <- worldclim_global(var = "bio", res = 10, 
                             path = "Data/Primary Data/ worldclim")
clim <- terra::subset(bio_data, c(1, 12))
Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
clim_br0 <- terra::crop(clim, Brazil)
clim_br <- terra::mask(clim_br0, Brazil)
### Log, because effects of preciptation tend to be log scale
v <- values(clim_br)[, 2]
logv <- sqrt(v)
values(clim_br)[, 2] <-  ifelse(is.infinite(logv), NA, logv)

## Pam data
pam_mammals <- lets.load("Data/Primary Data/PAM_crop_total.RData")

## Ecossystem service data
ES <- read.csv("Data/Primary Data/ES_11_TraitsNOVO.csv")
ES$total <- rowSums(ES[, -1])

# Environment PAM ---------------------------------------------------------
# Add variables to pam
envs <- lets.addvar(pam_mammals, clim_br, onlyvar = TRUE)
colnames(envs) <- c("Temperature", "Preciptation (Log)")

# Create Env pam
res <- lets.envpam(pam_mammals, envs,   n_bins = 20)
lets.plot.envpam(res,
                 world = TRUE)



# Map ecosystem services -------------------------------------------------
# Map Service 
res_map <- lets.maplizer.env(res, 
                             y = ES$total, 
                             z = ES$species,
                             func = sum) # Change here for 

# Relate to frequency
out <- lets.envcells(res)

abund_env   <- res_map$Matrix_env[, -(1:3), drop = FALSE]

# Mantain cells without zero
keep <- res$Presence_and_Absence_Matrix_env[, 1] & abund_env > 0
abund_env <- abund_env[keep]

tibble(abund_env, out[keep, -1]) %>%
  mutate(Frequency_log = log(Frequency)) %>%
  # select(-) %>% 
  pivot_longer(-abund_env) %>% 
  ggplot(aes(y = abund_env, x = value)) +
  geom_smooth() +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  geom_point() +
  facet_wrap(.~name, scales = "free")

