#Packages

library(letsR)
library(dplyr)
library(purrr)
library(tidyr)
library(gam)
library(readr)

# Creating an empty list for the data.
data_wide_list_SA <- list()

for (i in seq_along(es_names_SA)) {
  res_map_SA <- lets.maplizer.env(res_SA, y = ES_SA[[es_names_SA[i]]], z = ES_SA$species, func = sum)
  out_SA <- lets.envcells(res_SA)
  
  # Retrieving Temperature (column 2), Precipitation (column 3), and Abundance (column 4) directly from the environmental matrix.
  temp_climate_SA <- tibble(
    es_SA = es_names_SA[i],
    Temperature = res_map_SA$Matrix_env[, 2],
    `Precipitation (Log)` = res_map_SA$Matrix_env[, 3],
    Abundance_SA = res_map_SA$Matrix_env[, 4]
  )
  
  # Retrieving the geographic variables (Frequency, Isolation, etc.) from the out object.
  temp_geo_SA <- as_tibble(out_SA[keep_SA, -1])
  
  # Putting everything side-by-side in a single table and removing NAs.
  temp_wide_SA <- bind_cols(temp_climate_SA, temp_geo_SA) %>% drop_na()
  
  # Saving the completed table in the list.
  data_wide_list_SA[[i]] <- temp_wide_SA
}

# Combining the tables from the 11 services into a single data frame.
out_wide_SA <- bind_rows(data_wide_list_SA)

# -------------------------------------------------------------------------
# STATISTICAL MODEL
# -------------------------------------------------------------------------
results_gam_SA <- out_wide_SA %>%
  group_by(es_SA) %>%
  nest() %>%
  mutate(
    # Run the multivariate GAM for the 3 independent variables.
    model_SA = map(data, ~ gam(Abundance_SA ~ s(Temperature) + s(`Precipitation (Log)`) + s(Frequency), data = .x)),
    
    # Extracting the Pseudo-R2 (Deviance Explained) to assess the power of the model
    deviance_explained_SA = map_dbl(model_SA, ~ 1 - (.x$deviance / .x$null.deviance))
  )

# Final table with results
print(results_gam_SA %>% select(es_SA, deviance_explained_SA))

# Now there is an object called 'results_gam_SA'.
# To see the statistical summary of the model for the first service on the list, simply run:
# summary(results_gam_SA$model[[1]])

#############################################################################################################

# Extracting p values

# Extracting the p-values and organizing the table.
final_table_results_SA <- results_gam_SA %>%
  mutate(
    # Creating a new column that will extract the p-values from each model.
    p_values = map(model_SA, ~ {
      # Pull the parametric ANOVA table from the model.
      anova_table <- summary(.x)$parametric.anova
      
      # Creating a mini-table with the exact p-values of the 3 rows (Temp, Prec, Freq).
      tibble(
        p_Temperature = anova_table[1, "Pr(>F)"],
        p_Precipitation = anova_table[2, "Pr(>F)"],
        p_Frequency = anova_table[3, "Pr(>F)"]
      )
    })
  ) %>%
  # The `unnest` command unpacks the mini-table, placing the columns side by side.
  unnest(p_values) %>%
  # Selecting only the columns of interest, removing the column containing the raw model data.
  select(es_SA, deviance_explained_SA, p_Temperature, p_Precipitation, p_Frequency)

# View the final table 
print(final_table_results_SA)

# Saving the table 
#write_csv(final_table_results_SA, "Data\\Secondary Data\\Resultados_GAM_Capitulo2_SA.csv")

###################################################################################################

# Residuals

# Creating a PDF file to store residuals plots
pdf("Figures/Residuals_diagnosis_SA.pdf", width = 8, height = 8)

# Dividing the page
par(mfrow = c(2, 2))

# Loop (11 services)
for (i in 1:nrow(results_gam_SA)) {
  
  # Plotting 4 graphs for each service
  stats:::plot.lm(results_gam_SA$model_SA[[i]])
  
  # Placing the service name at the top of the page
  mtext(results_gam_SA$es_SA[[i]], side = 3, line = -2, outer = TRUE, font = 2, cex = 1.2)
}

# Saving the PDF file
dev.off()
par(mfrow = c(1, 1)) 

print("PDF successfully generated!")
