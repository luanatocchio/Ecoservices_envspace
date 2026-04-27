# -------------------------------------------------------------------------
# Replication of Coelho et al. (2023) framework for ecosystem services
# Multivariate GAMs with variance partitioning
# Predictors: Frequency, Isolation (Mean), PC1, PC2
# Response: Abundance of each ecosystem service (and total) in env. space
# -------------------------------------------------------------------------
# NOTE: This script assumes the following objects are already loaded
# from Envpam_SA_final.R:
#   - res_SA
#   - ES_SA (without top_down_regulation)
#   - pam_mammals_sub_SA
#   - clim_SA
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(letsR)
library(mgcv)
library(ggplot2)
library(tidyverse)
library(MetBrewer)
library(cowplot)

# -------------------------------------------------------------------------
# STEP 1: Prepare climate data and run PCA
# -------------------------------------------------------------------------

# Extract environmental variables per bin
out_SA_cells <- lets.envcells(res_SA)

# Extract temperature and precipitation directly from the environmental matrix
# This matrix already has 449 rows (one per non-empty bin)
env_matrix <- res_SA$Presence_and_Absence_Matrix_env

Temperature   <- env_matrix[, "Temperature"]
Precipitation <- env_matrix[, "Precipitation (Log)"]

envs_bins_keep <- data.frame(
  Temperature   = Temperature,
  Precipitation = Precipitation
)

# Check for NAs
cat("NAs - Temperature:", sum(is.na(envs_bins_keep[, 1])), "\n")
cat("NAs - Precipitation:", sum(is.na(envs_bins_keep[, 2])), "\n")

# Run PCA on scaled temperature and precipitation
pca_clim <- prcomp(envs_bins_keep, scale. = TRUE, center = TRUE)
summary(pca_clim)
cat("\nPCA loadings:\n")
print(pca_clim$rotation)

# Extract PC scores
pc_scores <- as.data.frame(pca_clim$x)
colnames(pc_scores) <- c("PC1", "PC2")

# -------------------------------------------------------------------------
# STEP 2: Build master data frame for models
# -------------------------------------------------------------------------

# Geography of climate metrics
# Filtrar out_SA_cells para manter apenas os bins presentes em env_matrix
cell_ids_env <- env_matrix[, "Cell_env"]

geo_clim <- out_SA_cells %>%
  as.data.frame() %>%
  filter(Cell_env %in% cell_ids_env) %>%
  select(Frequency, `Isolation (Mean)`) %>%
  rename(Isolation_Mean = `Isolation (Mean)`)

# Verificar dimensões antes de bind_cols
cat("Dimensões geo_clim:", nrow(geo_clim), "\n")
cat("Dimensões pc_scores:", nrow(pc_scores), "\n")
# -------------------------------------------------------------------------
# STEP 3: Multivariate GAMs and variance partitioning
# -------------------------------------------------------------------------

# Following Coelho et al. (2023) exactly:
# - Complete model with all 4 predictors
# - Reduced models removing one predictor at a time (fixing smoothing params)
# - Variance partitioning into unique and joint contributions

results_partition <- list()

for (es in es_names_model) {
  
  cat("\nRunning models for:", es, "\n")
  
  # Complete model (Poisson family for count data)
  b <- gam(
    as.formula(paste(es, "~ s(Frequency, k = 4) + s(Isolation_Mean, k = 4) + s(PC1, k = 4) + s(PC2, k = 4)")),
    data   = model_data,
    family = poisson()
  )
  
  # Null model
  b0 <- gam(
    as.formula(paste(es, "~ 1")),
    data   = model_data,
    family = poisson()
  )
  
  # Reduced models (fixing smoothing parameters from complete model)
  b1 <- gam(  # without Frequency
    as.formula(paste(es, "~ s(Isolation_Mean, k = 4) + s(PC1, k = 4) + s(PC2, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[2], b$sp[3], b$sp[4])
  )
  
  b2 <- gam(  # without Isolation
    as.formula(paste(es, "~ s(Frequency, k = 4) + s(PC1, k = 4) + s(PC2, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[1], b$sp[3], b$sp[4])
  )
  
  b3 <- gam(  # without PC1
    as.formula(paste(es, "~ s(Frequency, k = 4) + s(Isolation_Mean, k = 4) + s(PC2, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[1], b$sp[2], b$sp[4])
  )
  
  b4 <- gam(  # without PC2
    as.formula(paste(es, "~ s(Frequency, k = 4) + s(Isolation_Mean, k = 4) + s(PC1, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[1], b$sp[2], b$sp[3])
  )
  
  b5 <- gam(  # climate itself only (PC1 + PC2)
    as.formula(paste(es, "~ s(PC1, k = 4) + s(PC2, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[3], b$sp[4])
  )
  
  b6 <- gam(  # geography of climate only (Frequency + Isolation)
    as.formula(paste(es, "~ s(Frequency, k = 4) + s(Isolation_Mean, k = 4)")),
    data   = model_data,
    family = poisson(),
    sp     = c(b$sp[1], b$sp[2])
  )
  
  # Variance partitioning (proportion of null deviance explained)
  Total          <- (deviance(b0) - deviance(b))  / deviance(b0)
  FrequencyAlone <- (deviance(b1) - deviance(b))  / deviance(b0)
  IsolationAlone <- (deviance(b2) - deviance(b))  / deviance(b0)
  PC1Alone       <- (deviance(b3) - deviance(b))  / deviance(b0)
  PC2Alone       <- (deviance(b4) - deviance(b))  / deviance(b0)
  JointEffect    <- Total - (FrequencyAlone + IsolationAlone + PC1Alone + PC2Alone)
  
  GeoCLimPlusJoint  <- (deviance(b5) - deviance(b)) / deviance(b0)
  ClimSelfPlusJoint <- (deviance(b6) - deviance(b)) / deviance(b0)
  JointGeoClim      <- abs(GeoCLimPlusJoint  - (FrequencyAlone + IsolationAlone))
  JointClimSelf     <- abs(ClimSelfPlusJoint - (PC1Alone + PC2Alone))
  JointGeoxClimSelf <- abs(Total - JointGeoClim - JointClimSelf -
                             FrequencyAlone - IsolationAlone - PC1Alone - PC2Alone)
  
  results_partition[[es]] <- data.frame(
    Ecosystem_service           = es,
    Total_explained             = Total,
    Frequency_alone             = FrequencyAlone,
    Isolation_alone             = IsolationAlone,
    PC1_alone                   = PC1Alone,
    PC2_alone                   = PC2Alone,
    Joint_effect_all            = JointEffect,
    Joint_geography_of_climate  = JointGeoClim,
    Joint_climate_itself        = JointClimSelf,
    Joint_GeoClim_x_ClimSelf    = JointGeoxClimSelf,
    AdjR2_complete_model        = summary(b)$r.sq,
    Deviance_explained_complete = summary(b)$dev.expl
  )
  
  # Store complete model for partial residual plots
  assign(paste0("model_", es), b)
}

# Combine and save results
table_partition <- bind_rows(results_partition)
print(table_partition)

write.csv(table_partition,
          "Results/Results_Variance_Partitioning.csv",
          row.names = FALSE)

# -------------------------------------------------------------------------
# STEP 4: Partial residual plots (following Coelho et al. style)
# -------------------------------------------------------------------------

plot_partial_residuals <- function(model, data, es_name) {
  
  vars <- c("Frequency", "Isolation_Mean", "PC1", "PC2")
  
  # Build reduced models by removing one term at a time from the full model
  reduced_models <- list(
    update(model, . ~ . - s(Frequency,      k = 4)),
    update(model, . ~ . - s(Isolation_Mean, k = 4)),
    update(model, . ~ . - s(PC1,            k = 4)),
    update(model, . ~ . - s(PC2,            k = 4))
  )
  
  plot_list <- list()
  
  for (i in seq_along(vars)) {
    
    resid_data <- data.frame(
      x         = data[[vars[i]]],
      residuals = residuals(reduced_models[[i]])
    )
    
    p <- ggplot(resid_data, aes(x = x, y = residuals)) +
      geom_point(size  = 2,
                 col   = met.brewer("Egypt")[i],
                 alpha = 0.7) +
      stat_smooth(method  = "gam",
                  formula = y ~ s(x, k = 4),
                  col     = met.brewer("Hiroshige")[7],
                  alpha   = 0.10) +
      labs(x     = vars[i],
           y     = "Partial Residuals",
           title = paste(es_name, "|", vars[i])) +
      theme(
        panel.border     = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title       = element_text(size = 12),
        axis.text        = element_text(color = "black", size = 10),
        axis.line        = element_line(colour = "black")
      )
    
    plot_list[[i]] <- p
  }
  
  cowplot::plot_grid(plotlist = plot_list, nrow = 2, ncol = 2)
}

# Generate and save partial residual plots for each service
dir.create("Figures/Partial_Residuals", recursive = TRUE, showWarnings = FALSE)

for (es in es_names_model) {
  
  p <- plot_partial_residuals(
    model   = get(paste0("model_", es)),
    data    = model_data,
    es_name = es
  )
  
  ggsave(
    filename = paste0("Figures/Partial_Residuals/PartialResiduals_", es, ".tif"),
    plot     = p,
    width    = 20,
    height   = 15,
    units    = "cm",
    dpi      = 300
  )
  
  message(paste("Saved partial residuals for:", es))
}

# -------------------------------------------------------------------------
# STEP 5: Summary bar plot of variance partitioning (all services)
# -------------------------------------------------------------------------

# Clean service names for plotting
service_labels <- c(
  "cultural_service_charismatic_species"          = "Cultural Service",
  "ecotourism"                                    = "Ecotourism",
  "carrion_control"                               = "Carrion Control",
  "pollination"                                   = "Pollination",
  "seed_dispersal"                                = "Seed Dispersal",
  "pest_and_disease_control"                      = "Pest & Disease Control",
  "nutrient_transporting_horizontal_and_vertical" = "Nutrient Transporting",
  "ecosystem_engineering"                         = "Ecosystem Engineering",
  "rodent_control"                                = "Rodent Control",
  "disease_sentinelling"                          = "Disease Sentinelling",
  "total"                                         = "Total"
)

table_partition_long <- table_partition %>%
  select(Ecosystem_service,
         Frequency_alone,
         Isolation_alone,
         PC1_alone,
         PC2_alone,
         Joint_geography_of_climate,
         Joint_climate_itself,
         Joint_GeoClim_x_ClimSelf) %>%
  pivot_longer(cols      = -Ecosystem_service,
               names_to  = "Component",
               values_to = "Variance_explained") %>%
  mutate(
    Component = recode(Component,
                       "Frequency_alone"            = "Frequency",
                       "Isolation_alone"            = "Isolation",
                       "PC1_alone"                  = "PC1",
                       "PC2_alone"                  = "PC2",
                       "Joint_geography_of_climate" = "Joint: Geography of Climate",
                       "Joint_climate_itself"       = "Joint: Climate Itself",
                       "Joint_GeoClim_x_ClimSelf"   = "Joint: GeoClim x ClimSelf"
    ),
    Ecosystem_service = recode(Ecosystem_service, !!!service_labels)
  )

ggplot(table_partition_long,
       aes(x    = Ecosystem_service,
           y    = Variance_explained,
           fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = met.brewer("Hiroshige", n = 7)) +
  coord_flip() +
  labs(x    = "Ecosystem Service",
       y    = "Proportion of Deviance Explained",
       fill = "Component") +
  theme_bw() +
  theme(
    axis.text       = element_text(size = 10),
    axis.title      = element_text(size = 12),
    legend.position = "bottom"
  )

ggsave(
  "Figures/Variance_Partitioning_Summary.tif",
  width     = 25,
  height    = 15,
  units     = "cm",
  dpi       = 300,
  limitsize = FALSE
)

print("All Coelho replication analyses successfully completed!")