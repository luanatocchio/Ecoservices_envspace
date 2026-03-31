#Packages
library(letsR)

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