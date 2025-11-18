# Supplementary code to the paper 
# "Brazil Seed Transfer Zones: Supporting Seed Sourcing for Climate-Resilient Ecosystem Restoration"
# by Silva et al. (mateuscardosobio@gmail.com)

# Number of clusters/STZs
k = 48

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
setwd("path/folder") # CHANGE HERE

# Read raster layers
envi.pres = rast("bioclim_soil_1970_2000.tif") # Climatic and edaphic variables for the near present (1970-2000)

# Colours (manual specification)
col = c(
  "#83EE9A","#A5DA5B","#E245A5","#DFD433","#E0CC89","#DD9F64","#9B6CDD","#8D745F","#E7C1A4","#B7E5EB",
  "#E69092","#E33C40","#B2E7D1","#9FA1A5","#D9B0DF","#ED87BF","#74A6E5","#6CE5E8","#CBE990","#A7537B",
  "#E8F0DC","#B19DE3","#C0F039","#CA39EB","#45B07D","#60E5C8","#FFD700","#64BBB7","#6772A1","#468CA2",
  "#E0AE3B","#7A9442","#6547E4","#E37949","#DEADBF","#A6EAB5","#DD84DF","#E35775","#57EB7B","#EADDE7",
  "#64E840","#C7C5B6","#A8BF92","#E5E878","#E7EDB8","#B4C4E5","#6CCAEA","#5B74D7"
)
names(col) = as.character(1:48)

# GCM list
gcm.list = unique(c("GISS-E2-1-G", "FIO-ESM-2-0", "MPI-ESM1-2-HR", "ACCESS-CM2", "AWI-CM-1-1-MR", # Global https://esd.copernicus.org/articles/11/995/2020/
                    "EC-Earth3-Veg", "INM-CM4-8", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "CMCC-ESM2")) # Regional https://iopscience.iop.org/article/10.1088/2752-5295/ad3fdb

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
clus.spat = terra::focal(x = clus.spat, w = matrix(1, 3, 3), fun = "modal", na.rm = TRUE)

# Resample it
clus.spat = resample(clus.spat, envi.pres[[1]])

# Crop it
clus.spat = crop(clus.spat, envi.pres[[1]])

# Mask it
clus.spat = mask(clus.spat, envi.pres[[1]])

# Round
clus.spat = terra::round(clus.spat)

# Save files
writeRaster(clus.spat, "STZ_1970-2000.tif", overwrite = T) # STZ present
writeRaster(scor.spat, "PCA_Scores.tif", overwrite = T) # PCA scores
write.csv(data.frame(pca[["rotation"]]), "PCA_Loadings.csv") # PCA loadings
write.csv(data.frame(summary(pca)$importance), "PCA_Importance.csv") # PCA importance

# 4. Future projection ----

for(ssp in c("ssp126", "ssp585")) {
  
  for(time in c("2041-2060", "2081-2100")) {
    
    for(gcm in gcm.list) {
      
      # Read future climate
      envi.futu = rast(paste("wc2.1_2.5m_bioc", gcm, ssp, paste0(time, ".tif"), sep = "_")) # Ensemble
      
      # Crop
      envi.futu = terra::crop(envi.futu, envi.pres[[1]])
      
      # Mask
      envi.futu = terra::mask(envi.futu, envi.pres[[1]])
      
      # Remove BIO18 and BIO19 due to discontinuity
      envi.futu = envi.futu[[-c(18, 19)]]
      
      # Stack environmental raster layers
      envi.futu = c(envi.futu, envi.pres[[c("bdod", "cfvo", "nitrogen", "phh2o", "sand", "silt", "clay", "soc")]])
      
      # Match the names
      names(envi.futu) = names(envi.pres)
      
      # Transform to data frame
      futu.tabl = terra::as.data.frame(envi.futu, xy = T, na.rm = T)
      
      # Calculate PCA scores for the future
      pca.futu = predict(pca, newdata = futu.tabl)
      
      # Get the first n PCs
      numb.futu = data.frame(pca.futu[, 1:numb])
      
      # Scale coordinates so they have a variance equal to average variance of the PCs
      coor.futu = data.frame(scale(futu.tabl[, c("x", "y")]))
      
      # Rename
      coor.futu = coor.futu %>%
        rename(long = x,
               lat = y)
      
      # Merge the scores and the rescaled coordinates
      scor.futu = bind_cols(numb.futu, coor.futu)
      
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
      clus.spat.futu = terra::focal(x = clus.spat.futu, w = matrix(1, 3, 3), fun = "modal", na.rm = TRUE)
      
      # Resample it
      clus.spat.futu = resample(clus.spat.futu, envi.pres[[1]])
      
      # Crop it
      clus.spat.futu = crop(clus.spat.futu, envi.pres[[1]])
      
      # Mask it
      clus.spat.futu = mask(clus.spat.futu, envi.pres[[1]])
      
      # Save raster
      writeRaster(clus.spat.futu, paste("STZ", gcm, ssp, paste0(time, ".tif"), sep = "_"), overwrite = T)
      
    }
    
  }
  
}

