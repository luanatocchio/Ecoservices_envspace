# Code to map species across spaces

# Packages ----------------------------------------------------------------
## install.packages("devtools")
## library(devtools)
## install_github("macroecology/letsR")
library(letsR)
library(terra)
library(geodata)
library(tidyverse)
library(gam)

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
res <- lets.envpam(pam_mammals, envs,   n_bins = 30)
lets.plot.envpam(res,
                 world = TRUE)



# Map ecosystem services -------------------------------------------------
# Map Service 
es_names <- names(ES)[-1]
data_plot <- list()
gam_results <- list()
# Mantain cells without zero
keep <- res$Presence_and_Absence_Matrix_env[, 1]

for (i in seq_along(es_names)) {
  res_map <- lets.maplizer.env(
    res,
    y = ES[[es_names[i]]],
    # Change for each ecosystem service
    z = ES$species,
    func = sum
  ) # Change here for mean, sum, variability
  out <- lets.envcells(res)
  # lets.plot.envpam(res_map) # for ploting (if you want to save it)
  
  abund_env   <- res_map$Matrix_env[, -(1:3), drop = FALSE]
  data_plot[[i]] <- 
    tibble(es = es_names[i], Abundance = abund_env[keep, 1], out[keep, -1]) %>%
    drop_na() %>% 
    # mutate(Frequency_log = log(Frequency)) %>%
    select(c(1:3, 7, 15)) %>% # remove variables 
    pivot_longer(3:5)
  
  # STATISTICAL MODEL
  # gam_results[[i]] <- gam(Abundance ~ ) Dplyr roda modelo para as 3 variaveis e guardar numa (salvar numa lista)
}

# Combine all
out <- do.call(bind_rows, data_plot)

# Relate to frequency


# Plot
out %>%
  ggplot(aes(y = Abundance, x = value)) +
  geom_smooth() +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  geom_point() +
  theme_bw() +
  facet_grid(es~name, scales = "free")


ggsave("Figures/es_env_relationships.tif",
       height = 60,
       width = 20, 
       units = "cm")
