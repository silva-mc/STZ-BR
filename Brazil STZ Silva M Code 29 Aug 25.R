# Supplementary code to the paper 
# "Brazil Transfer Seed Zones: Supporting Seed Sourcing for Ecological Restoration under Current and Future Climates"
# by Silva et al. (mateuscardosobio@gmail.com)

# Number of clusters/STZs
k = 19

# 1. Initialization ----

# Required packages
library(terra)
library(factoextra)
library(cluster)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)

# Set working directory
setwd(paste("C:/Users",
            Sys.info()[["user"]],
            "OneDrive - University of Exeter (1)/Consultancy/ISA (Jan-Sep 25)/Submission/Version 1/Code",
            sep = "/"))

# Read raster layers
envi.pres = rast("wc2.1_sg2.0_30s_1970_2000.tif") # Climatic, edaphic and topographic variables for the near present (1970-2000)

# Colours (add more if k > 25)
col = c(
  "deepskyblue", "darkred", "gold", "lightgreen", "magenta", 
  "floralwhite", "coral1", "darkgreen", "orange", "wheat", 
  "purple", "slategray", "darkolivegreen1", "peru", "lavender", 
  "limegreen", "yellow", "royalblue", "lightskyblue1", "lightpink", 
  "darkturquoise", "saddlebrown", "mediumslateblue", "orchid", "tomato3"
)

# 2. Data processing ----

# Transform to data frame
envi.tabl = terra::as.data.frame(envi.pres, xy = T, na.rm = T)

# Create an PCA
pca = prcomp(envi.tabl[3:27], scale. = T) # PCA excluding lat and long

# Get PCA scores
pca.scor = pca$x

# Get eigenvalues
pca.eige = get_eig(pca)

# List PCs with eigenvalue > 1
numb = which(pca.eige$eigenvalue > 1)

# Get the number of PCs
numb = numb[length(numb)]

# Get the first n PCs
scor.numb = data.frame(pca.scor[, 1:numb])

# Scale coordinates so they have a variance equal to average variance of the PCs
coor = data.frame(scale(envi.tabl[, c("x", "y")]))

# Rename
coor = coor %>%
  rename(long = x,
         lat = y)

# Merge the scores and the rescaled coordinates
scor = bind_cols(scor.numb, coor)

# Create an empty raster
scor.spat = rast()

# Make it spatial
for(z in names(scor)) {
  
  # Filter only xyz
  scor.z = bind_cols(scor, envi.tabl[, c("x", "y")]) 
  scor.z = scor.z[, c("x", "y", z)]
  scor.z = rast(scor.z,
                type = "xyz")
  
  # Add to the empty raster
  scor.spat[[z]] = scor.z
}

# 3. Cluster analysis ----

clus = clara(x = scor,
             k = k,
             metric = "euclidean",
             samples = 50,
             sampsize = 2000,
             stand = T,
             pamLike = T,
             correct.d = T)

# Add coordinate to the cluster
clus.tabl = envi.tabl %>%
  dplyr::select(x, y) %>%
  bind_cols(clus$cluster) %>%
  rename("cluster" = "...3")

# Make it spatial
clus.spat = rast(clus.tabl,
                 type = "xyz")

# Assign a CRS
crs(clus.spat) = crs(envi.pres)

# Moving window analysis to remove isolated pixels
clus.spat = terra::focal(x = clus.spat, w = matrix(1, 11, 11), fun = "modal", na.rm = TRUE)

# Resample it
clus.spat = resample(clus.spat, envi.pres[[1]])

# Crop it
clus.spat = crop(clus.spat, envi.pres[[1]])

# Mask it
clus.spat = mask(clus.spat, envi.pres[[1]])

# Round
clus.spat = terra::round(clus.spat)

# Save files
writeRaster(clus.spat, "STZ_Present.tif", overwrite = T) # STZ present
writeRaster(scor.spat, "PCA_Scores.tif", overwrite = T) # PCA scores
write.csv(data.frame(pca[["rotation"]]), "PCA_Loadings.csv") # PCA loadings
write.csv(data.frame(summary(pca)$importance), "PCA_Importance.csv") # PCA importance

# 4. Future projection ----

# Loop for SSP1 and SSP5 scenarios
for(ssp in c("SSP1", "SSP5")) {
  
  # Read future climate
  envi.futu = rast(paste0("cmip6_sg2.0_30s_2081_2100_", ssp, ".tif")) # Ensemble climatic, edaphic and topographic variables for the far future (2081-2100)
  
  # Crop it
  envi.futu = terra::crop(envi.futu, clus.spat)
  
  # Transform to data frame
  futu.tabl = terra::as.data.frame(envi.futu, xy = T, na.rm = T)
  
  # Calculate PCA scores for the future
  pca.futu = predict(pca, newdata = futu.tabl)
  
  # Get the first n PCs
  numb.futu = data.frame(pca.futu[, 1:numb])
  
  # Get raw coordinates
  coor.raw = envi.tabl[, c("x", "y")]
  
  # Merge raw and scaled coordinates
  coor.merg = bind_cols(coor.raw, coor) %>%
    unite("key", c(x, y), remove = T)
  
  # Bind future coordinates with PCA scores
  coor.scor = bind_cols(futu.tabl, numb.futu)
  
  # Do the same for future
  futu.merg = coor.scor %>%
    unite("key", c(x, y), remove = F)
  
  # Merge datasets
  scor.futu = left_join(futu.merg %>%
                          select(key, starts_with("PC")), 
                        coor.merg,
                        by = "key") %>%
    drop_na() %>%
    dplyr::select(-key)
  
  # # Merge the scores and the rescaled coordinates
  # scor.futu = bind_cols(numb.futu, coor)
  
  # Create an empty raster
  futu.spat = rast()
  
  # Make it spatial
  for(z in names(scor.futu)) {
    
    # Filter only xyz
    scor.futu.z = bind_cols(scor.futu, futu.tabl[, c("x", "y")]) 
    scor.futu.z = scor.futu.z[, c("x", "y", z)]
    scor.futu.z = rast(scor.futu.z,
                       type = "xyz")
    
    # Add to the empty raster
    futu.spat[[z]] = scor.futu.z
  }
  
  # Assign clusters
  assign_clusters = function(scores, medoids) {
    apply(scores, 1, function(row) {
      which.min(colSums((t(medoids) - row)^2))
    })
  }
  
  # Future clusters
  futu.clus = assign_clusters(scor.futu, clus$medoids)
  
  # Get coordinates
  futu.clus = cbind(futu.tabl[, c("x", "y")], cluster = futu.clus)
  
  # Get future clusters
  clus.spat.futu = rast(futu.clus,
                        type = "xyz")
  
  # Moving window analysis to remove isolated pixels
  clus.spat.futu = terra::focal(x = clus.spat.futu, w = matrix(1, 11, 11), fun = "modal", na.rm = TRUE)
  
  # Resample it
  clus.spat.futu = resample(clus.spat.futu, envi.pres[[1]])
  
  # Crop it
  clus.spat.futu = crop(clus.spat.futu, envi.pres[[1]])
  
  # Mask it
  clus.spat.futu = mask(clus.spat.futu, envi.pres[[1]])
  
  # Round it
  clus.spat.futu = terra::round(clus.spat.futu, digits = 0) # SSP1
  
  # Save raster
  writeRaster(clus.spat.futu, paste0("STZ_2100_", ssp, ".tif"), overwrite = T)
  
}

