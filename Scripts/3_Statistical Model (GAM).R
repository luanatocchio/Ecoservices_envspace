#Packages

library(dplyr)
library(purrr)
library(tidyr)
library(gam)
library(readr)

# Creating an empty list for the data.
data_wide_list <- list()

for (i in seq_along(es_names)) {
  res_map <- lets.maplizer.env(res, y = ES[[es_names[i]]], z = ES$species, func = sum)
  out <- lets.envcells(res)
  
  # Retrieving Temperature (column 2), Precipitation (column 3), and Abundance (column 4) directly from the environmental matrix.
  temp_climate <- tibble(
    es = es_names[i],
    Temperature = res_map$Matrix_env[, 2],
    `Preciptation (Log)` = res_map$Matrix_env[, 3],
    Abundance = res_map$Matrix_env[, 4]
  )
  
  # Retrieving the geographic variables (Frequency, Isolation, etc.) from the out object.
  temp_geo <- as_tibble(out[keep, -1])
  
  # Putting everything side-by-side in a single table and removing NAs.
  temp_wide <- bind_cols(temp_climate, temp_geo) %>% drop_na()
  
  # Saving the completed table in the list.
  data_wide_list[[i]] <- temp_wide
}

# Combining the tables from the 11 services into a single data frame.
out_wide <- bind_rows(data_wide_list)

# -------------------------------------------------------------------------
# STATISTICAL MODEL
# -------------------------------------------------------------------------
results_gam <- out_wide %>%
  group_by(es) %>%
  nest() %>%
  mutate(
    # Run the multivariate GAM for the 3 independent variables.
    model = map(data, ~ gam(Abundance ~ s(Temperature) + s(`Preciptation (Log)`) + s(Frequency), data = .x)),
    
    # Extracting the Pseudo-R2 (Deviance Explained) to assess the power of the model
    deviance_explained = map_dbl(model, ~ 1 - (.x$deviance / .x$null.deviance))
  )

# Final table with results
print(results_gam %>% select(es, deviance_explained))

# Now there is an object called 'gam_results'.
# To see the statistical summary of the model for the first service on the list, simply run:
# summary(results_gam$model[[1]])

#############################################################################################################

# Extracting p values

# Extracting the p-values and organizing the table.
final_table_results <- results_gam %>%
  mutate(
    # Creating a new column that will extract the p-values from each model.
    p_values = map(model, ~ {
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
  select(es, deviance_explained, p_Temperature, p_Precipitation, p_Frequency)

# View the final table 
print(final_table_results)

# Saving the table 
#write_csv(final_table_results, "Data\\Secondary Data\\Resultados_GAM_Capitulo2.csv")

