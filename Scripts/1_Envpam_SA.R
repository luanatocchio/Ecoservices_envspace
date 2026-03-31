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
library(sf)

# Load Data ---------------------------------------------------------------
## Env data
bio_data <- worldclim_global(var = "bio", res = 10, 
                             path = "Data/Primary Data/ worldclim")
clim <- terra::subset(bio_data, c(1, 12))

# Loading South America Polygon
SA <- st_read("Data\\Primary Data\\SA\\Lim_america_do_sul_2021.shp")
SA <- st_transform(SA, crs(clim))
#terra::ext(SA)

clim_SA0 <- terra::crop(clim, SA)
clim_SA <- terra::mask(clim_SA0, SA)

### Log, because effects of precipitation tend to be log scale
v <- values(clim_SA)[, 2]
logv <- sqrt(v)
values(clim_SA)[, 2] <-  ifelse(is.infinite(logv), NA, logv)

## Pam data

### Creating Pam data for South America

mammals <- st_read("Data\\Primary Data\\mammals_polygons\\data_0.shp")
mammals <- mammals %>%
  rename(SCINAME = SCI_NAME)


#pam_SA <- lets.presab(mammals, xmn = -83.9361138954408, xmx = -28.8477703530605, ymn = -55.9866984156709, ymx = 13.5854011963317, res = 0.5, count = TRUE)

#lets.save(pam_SA,file = "Data\\Secondary Data\\pam_SA.RData")
pam_mammals_SA <- lets.load(file = "Data\\Secondary Data\\pam_SA.RData") 
plot(pam_mammals_SA)

############################################################################################################

# Loading species data

vale_species <- read.csv("Data/Primary Data/vale_species_last.csv")[[1]]

#str(vale_species)

########################################################################################################################################################################

# Subset a presence-absence object based on species character vector

pam_mammals_sub_SA <- lets.subsetPAM(x = pam_mammals_SA, names = vale_species, remove.cells = FALSE)
plot(pam_mammals_sub_SA, xlab = "Longitude", ylab = "Latitude",
     main = "Mammals species richness")

str(pam_mammals_sub_SA)

# Taking the test, to check for possible discrepancies

differences_vale_pam_mammals_sub_SA <- setdiff(vale_species, pam_mammals_sub_SA$Species_name)
differences_pam_mammals_sub_SA_vale <- setdiff(pam_mammals_sub_SA$Species_name, vale_species)

# lets.save(pam_mammals_sub_SA, file = "Data\\Secondary Data\\pam_mammals_sub_SA.RData")
pam_mammals_sub_SA <- lets.load(file = "Data\\Secondary Data\\pam_mammals_sub_SA.RData") 

################################################################################################################################

## Ecossystem service data
ES_SA <- read.csv("Data/Primary Data/ES_11_TraitsNOVO.csv")
ES_SA$total <- rowSums(ES_SA[, -1])

# Environment PAM ---------------------------------------------------------
# Add variables to pam
envs_SA <- lets.addvar(pam_mammals_sub_SA, clim_SA, onlyvar = TRUE)
colnames(envs_SA) <- c("Temperature", "Precipitation (Log)")

# Create Env pam
res_SA <- lets.envpam(pam_mammals_sub_SA, envs_SA,   n_bins = 30)
lets.plot.envpam(res_SA,
                 world = TRUE)

# Map ecosystem services -------------------------------------------------
# Map Service 
es_names_SA <- names(ES_SA)[-1]
data_plot_SA <- list()
gam_results_SA <- list()

# Mantain cells without zero
keep_SA <- res_SA$Presence_and_Absence_Matrix_env[, 1]

for (i in seq_along(es_names_SA)) {
  res_map_SA <- lets.maplizer.env(
    res_SA,
    y = ES_SA[[es_names_SA[i]]],
    # Change for each ecosystem service
    z = ES_SA$species,
    func = sum
  ) # Change here for mean, sum, variability
  out_SA <- lets.envcells(res_SA)
  # lets.plot.envpam(res_map) # for ploting (if you want to save it)
  
  abund_env_SA   <- res_map_SA$Matrix_env[, -(1:3), drop = FALSE]
  data_plot_SA[[i]] <- 
    tibble(es_SA = es_names_SA[i], Abundance_SA = abund_env_SA[, 1], out_SA[keep_SA, -1]) %>%
    drop_na() %>% 
    # mutate(Frequency_log = log(Frequency)) %>%
    select(c(1:3, 7, 15)) %>% # Select variables 
    pivot_longer(3:5)
  
  # STATISTICAL MODEL
  # gam_results[[i]] <- gam(Abundance ~ ) Dplyr roda modelo para as 3 variaveis e guardar numa (salvar numa lista)
}

# Combine all
out_SA <- do.call(bind_rows, data_plot_SA)


# Plot
out_SA %>%
  ggplot(aes(y = Abundance_SA, x = value)) +
  geom_smooth(method = "gam") +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  geom_point() +
  theme_bw() +
  facet_grid(es_SA~name, scales = "free")


ggsave("Figures/es_env_relationships_SA.tif",
       height = 60,
       width = 20, 
       units = "cm")




