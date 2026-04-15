# -------------------------------------------------------------------------
# Mapping species across environmental spaces and running GAM models - BRAZIL
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(letsR)
library(terra)
library(geodata)
library(tidyverse)
library(sf)
library(broom)
library(mgcv)  ]

# Load Data ---------------------------------------------------------------
## Env data (WorldClim)
bio_data <- worldclim_global(var = "bio", res = 10, 
                             path = "Data/Primary Data/worldclim")
clim <- terra::subset(bio_data, c(1, 12))

# Brazil Polygon
Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ] 
clim_br0 <- terra::crop(clim, Brazil)
clim_br <- terra::mask(clim_br0, Brazil)

### Log transformation for precipitation
v <- values(clim_br)[, 2]
logv <- sqrt(v)
values(clim_br)[, 2] <- ifelse(is.infinite(logv), NA, logv)

## PAM data
pam_mammals <- lets.load("Data/Primary Data/PAM_crop_total.RData")

## Ecosystem service data
ES <- read.csv("Data/Primary Data/ES_11_TraitsNOVO.csv")
ES$total <- rowSums(ES[, -1])

# Environment PAM ---------------------------------------------------------
envs <- lets.addvar(pam_mammals, clim_br, onlyvar = TRUE)
colnames(envs) <- c("Temperature", "Precipitation (Log)")

# Create Env PAM
res <- lets.envpam(pam_mammals, envs, n_bins = 30)

# Map Ecosystem Services and Run Models -----------------------------------
es_names <- names(ES)[-1]
data_plot <- list()
gam_results <- list()

# Maintain cells without zero
keep <- res$Presence_and_Absence_Matrix_env[, 1]

# --- RESIDUALS REPORT ---
# Creates a single PDF with diagnostics for all models
pdf("Figures/Residuals_Diagnosis_Report_BRAZIL.pdf", width = 10, height = 10)
par(mfrow = c(2, 2)) # 4 plots per page

for (i in seq_along(es_names)) {
  
  # 1. Map the service
  res_map <- lets.maplizer.env(
    res,
    y = ES[[es_names[i]]],
    z = ES$species,
    func = sum
  ) 
  
  out_cells <- lets.envcells(res)
  abund_env <- res_map$Matrix_env[, -(1:3), drop = FALSE]
  
  # 2. Prepare clean data table
  data_plot[[i]] <- tibble(
    es_name = es_names[i], 
    Abundance = abund_env[, 1], 
    out_cells[keep, -1],
    .name_repair = "unique"
  ) %>%
    drop_na() %>% 
    # Selecting exact environmental space variables
    select(es_name, Abundance, Frequency, `Isolation (Mean)`, `Frequency Weighted Distance`) %>% 
    pivot_longer(cols = 3:5)
  
  # 3. Statistical Model (GAM) - Univariate
  gam_results[[i]] <- data_plot[[i]] %>%
    group_by(name) %>%
    summarise(
      model = list(mgcv::gam(Abundance ~ s(value), data = pick(everything())))
    )
  
  # 4. Generate Residual Plots for the report
  for(j in 1:nrow(gam_results[[i]])) {
    mod_temp <- gam_results[[i]]$model[[j]]
    var_name <- gam_results[[i]]$name[j]
    
    # Run diagnostic
    mgcv::gam.check(mod_temp)
    
    # Title for the specific page
    mtext(paste("Service:", es_names[i], "| Variable:", var_name), 
          side = 3, line = -1.5, outer = TRUE, font = 2, cex = 0.8)
  }
}

dev.off() # Close the PDF

# Final Table of Results --------------------------------------------------
names(gam_results) <- es_names
all_models <- bind_rows(gam_results, .id = "Ecosystem_service")

table_gam_BR <- all_models %>%
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
write.csv(table_gam_BR, "Results/Results_GAM_Statistics_BRAZIL.csv", row.names = FALSE)

# Visualization -----------------------------------------------------------

# Combine all
out_final <- do.call(bind_rows, data_plot)

# Plot
out_final %>%
  ggplot(aes(y = Abundance, x = value)) +
  geom_smooth(method = "gam") +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  geom_point() +
  theme_bw() +
  facet_grid(es_name~name, scales = "free")

#ggsave("Figures/Es_env_relationships_.tif",
#height = 60,
#width = 20, 
#units = "cm")

# Creating a folder to store the combined images.
dir.create("Figures/Maps_env", recursive = TRUE, showWarnings = FALSE)

# Service names
es_names <- names(ES)[-1]

for (i in seq_along(es_names)) {
  
  # Calculating the abundance of ecosystem services in the environmental space.
  res_map <- lets.maplizer.env(
    x = res,
    y = ES[[es_names[i]]],
    z = ES$species,
    func = sum
  )
  
  # Maps
  # Note that the width is double the height to fit the two maps side by side without creasing.
  png(filename = paste0("Figures/Maps_env/", es_names[i], "_env.png"), 
      width = 3000, height = 1500, res = 300)
  
  # Plot the double panel (Geographic + Environmental)
  lets.plot.envpam(res_map, world = TRUE)
  
  # saving images to folder.
  dev.off()
}

print("Maps successfully generated!")