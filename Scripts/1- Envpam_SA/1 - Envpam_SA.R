# -------------------------------------------------------------------------
# Mapping species across environmental spaces and running GAM models
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(letsR)
library(terra)
library(geodata)
library(tidyverse)
library(sf)
library(broom)
library(mgcv)  

# Load Data ---------------------------------------------------------------
## Env data (WorldClim)
bio_data <- worldclim_global(var = "bio", res = 10, 
                             path = "Data/Primary Data/worldclim")
clim <- terra::subset(bio_data, c(1, 12))

# Loading South America Polygon
SA <- st_read("Data/Primary Data/SA/Lim_america_do_sul_2021.shp")
SA <- st_transform(SA, crs(clim))

clim_SA0 <- terra::crop(clim, SA)
clim_SA <- terra::mask(clim_SA0, SA)

### Log transformation for precipitation
v <- values(clim_SA)[, 2]
logv <- sqrt(v)
values(clim_SA)[, 2] <- ifelse(is.infinite(logv), NA, logv)

### Creating Pam data for South America
#mammals <- st_read("Data\\Primary Data\\mammals_polygons\\data_0.shp")
#mammals <- mammals %>%
  #rename(SCINAME = SCI_NAME)
#pam_SA <- lets.presab(mammals, xmn = -83.9361138954408, xmx = -28.8477703530605, ymn = -55.9866984156709, ymx = 13.5854011963317, res = 0.5, count = TRUE)
#lets.save(pam_SA,file = "Data\\Secondary Data\\pam_SA.RData")

# Loading saved PAM
pam_mammals_SA <- lets.load(file = "Data/Secondary Data/pam_SA.RData") 

# Loading and subsetting species
vale_species <- read.csv("Data/Primary Data/vale_species_last.csv")[[1]]
pam_mammals_sub_SA <- lets.subsetPAM(x = pam_mammals_SA, names = vale_species, remove.cells = FALSE)

## Ecosystem Service Data
ES_SA <- read.csv("Data/Primary Data/ES_11_TraitsNOVO.csv")
ES_SA$total <- rowSums(ES_SA[, -1])

# Environment PAM ---------------------------------------------------------
envs_SA <- lets.addvar(pam_mammals_sub_SA, clim_SA, onlyvar = TRUE)
colnames(envs_SA) <- c("Temperature", "Precipitation (Log)")

# Create Env PAM
res_SA <- lets.envpam(pam_mammals_sub_SA, envs_SA, n_bins = 30)

# Map Ecosystem Services and Run Models -----------------------------------
es_names_SA <- names(ES_SA)[-1]
data_plot_SA <- list()
gam_results_SA <- list()

# Maintain cells without zero
keep_SA <- res_SA$Presence_and_Absence_Matrix_env[, 1]

# --- RESIDUALS REPORT ---
# This opens a single PDF to store all 36 diagnostic pages
pdf("Figures/Residuals_Diagnosis_Report_SA.pdf", width = 10, height = 10)
par(mfrow = c(2, 2)) # 4 plots per page

for (i in seq_along(es_names_SA)) {
  
  # 1. Map the service
  res_map_SA <- lets.maplizer.env(
    res_SA,
    y = ES_SA[[es_names_SA[i]]],
    z = ES_SA$species,
    func = sum
  ) 
  
  out_SA_cells <- lets.envcells(res_SA)
  abund_env_SA <- res_map_SA$Matrix_env[, -(1:3), drop = FALSE]
  
  # 2. Prepare data table 
  data_plot_SA[[i]] <- tibble(
    es_SA = es_names_SA[i], 
    Abundance_SA = abund_env_SA[, 1], 
    out_SA_cells[keep_SA, -1],
    .name_repair = "unique" 
  ) %>%
    drop_na() %>% 
    select(es_SA, Abundance_SA, Frequency, `Isolation (Mean)`, `Frequency Weighted Distance`) %>% 
    pivot_longer(cols = 3:5)
  
  # 3. Run Statistical Model (GAM)
  gam_results_SA[[i]] <- data_plot_SA[[i]] %>%
    group_by(name) %>%
    summarise(
      model = list(mgcv::gam(Abundance_SA ~ s(value), data = pick(everything())))
    )
  
  # 4. Generate Residual Plots for each of the 3 variables
  for(j in 1:nrow(gam_results_SA[[i]])) {
    mod_temp <- gam_results_SA[[i]]$model[[j]]
    var_nome <- gam_results_SA[[i]]$name[j]
    
    # Check diagnostics
    mgcv::gam.check(mod_temp)
    
    # Add title to the page
    mtext(paste("Service:", es_names_SA[i], "| Variable:", var_nome), 
          side = 3, line = -1.5, outer = TRUE, font = 2, cex = 0.8)
  }
}

dev.off() # Close and save the PDF report

# Final Table of Results --------------------------------------------------
names(gam_results_SA) <- es_names_SA
all_models <- bind_rows(gam_results_SA, .id = "Ecosystem_service")

table_gam <- all_models %>%
  mutate(coefficients = map(model, broom::tidy)) %>%
  mutate(quality = map(model, broom::glance)) %>%
  unnest(c(coefficients, quality), names_sep = "_") %>%
  select(
    Ecosystem_service,
    Variable = name,
    Degrees_of_freedom = coefficients_edf,
    F_statistic = coefficients_statistic,
    P_value = coefficients_p.value,
    Adjusted_R_Squared = quality_adj.r.squared,
    Deviance_Explained = quality_deviance
  )

# Save Statistics
#write.csv(table_gam, "Results/Results_GAM_Statistics.csv", row.names = FALSE)

# Visualization -----------------------------------------------------------

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

#ggsave("Figures/Es_env_relationships_SA.tif",
# height = 60,
# width = 20, 
# units = "cm")

# Creating a folder to store the combined images for South America.
dir.create("Figures/Maps_env_SA", recursive = TRUE, showWarnings = FALSE)

# Service names (garantindo que estamos usando os dados da América do Sul)
es_names_SA <- names(ES_SA)[-1]

for (i in seq_along(es_names_SA)) {
  
  # Calculating the abundance of ecosystem services in the environmental space.
  res_map_SA <- lets.maplizer.env(
    x = res_SA,
    y = ES_SA[[es_names_SA[i]]],
    z = ES_SA$species,
    func = sum
  )
  
  # Maps
  # Note that the width is double the height to fit the two maps side by side.
  png(filename = paste0("Figures/Maps_env_SA/", es_names_SA[i], "_env_SA.png"), 
      width = 3000, height = 1500, res = 300)
  
  # Plot the double panel (Geographic + Environmental)
  lets.plot.envpam(res_map_SA, world = TRUE)
  
  # saving images to folder.
  dev.off()
}

print("South America Maps successfully generated!")