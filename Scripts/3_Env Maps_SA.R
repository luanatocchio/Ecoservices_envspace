#Packages
library(letsR)

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