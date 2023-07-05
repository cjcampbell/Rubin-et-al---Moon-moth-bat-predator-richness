source("R/00_setup.R")

library(data.table)
library(sf)
library(terra)

proj4_eqe <- "+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs"
proj4_wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Load coordinates to be used. ----
if(!file.exists(file.path(wd$data, "mothOccurrences.csv"))) {
  read.csv("bin/Modified_fulldata_2-24-23.csv") %>% 
    dplyr::select(decimalLongitude, decimalLatitude, genus, species, specificEpithet) %>% 
    fwrite(file.path(wd$data, "mothOccurrences.csv"), row.names = F)
}

mothOccurrences <- fread(file.path(wd$data, "mothOccurrences.csv")) %>% 
  dplyr::mutate(id = row_number())
mothOccurrences_sf <- mothOccurrences %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = proj4_wgs) %>% 
  st_transform(proj4_eqe) %>% 
  full_join(mothOccurrences)

# Load surface ------------------------------------------------------------
# Bat richness surface.
mySurface <- rast(file.path(wd$out, "combinedContinuousSurface.tif"))
set.crs(mySurface, proj4_eqe)


# Extract values ----------------------------------------------------------

mothOccurrences$batRichness <- terra::extract(mySurface, mothOccurrences_sf, ID = F)
mothOccurrences_sf$batRichness <- mothOccurrences$batRichness

# Find nearest cells -----------------------------------------------------
# In cases where the point is quite close to an existing cell, 
# return the value of the nearest cell.

# First, fortify surface.
mySurface_df <- mySurface %>% 
  as.data.frame(xy=T) 

# For each occurence with NA value, find distance to nearest cell centroid.
for(i in 1:nrow(mothOccurrences)) {
  if(!is.na(mothOccurrences$batRichness[i])) next
  
  writeLines(paste("Working on", mothOccurrences$id[i]))
  # Find coordinates of occurrence record.
  whichPts <- unlist(mothOccurrences_sf[i,]$geometry)
  # Filter to any cells with values within 25km.
  relevantCells_df_100 <- mySurface_df %>% 
    dplyr::filter(abs(x - whichPts[1]) < 100e3 & abs(y - whichPts[2]) < 100e3) 
  relevantCells_df <- relevantCells_df_100 %>% 
    dplyr::filter(abs(x - whichPts[1]) < 25e3 & abs(y - whichPts[2]) < 25e3) 
  if(nrow(relevantCells_df) == 0) next
  relevantCells <- relevantCells_df %>% 
    st_as_sf(coords = c("x", "y"), crs = proj4_eqe) 
  closestCell <- relevantCells %>% 
    st_distance(mothOccurrences_sf[i,]) %>% 
    which.min
  value <- unlist(relevantCells[closestCell,])[1]
  
  # Optionally, plot each of these to make sure it checks out.
  # p <- ggplot() + 
  #   geom_tile(relevantCells_df_100, mapping = aes(x=x,y=y,fill=sum, color = sum)) + 
  #   geom_sf(mothOccurrences_sf[i,], mapping = aes(), color = "red") +
  #   geom_sf(relevantCells[closestCell,], mapping = aes(), color = "yellow") +
  #   coord_sf(crs = proj4_eqe, xlim = c())
  # plot(p)
  # Sys.sleep(1)
  # readline(prompt="Press [enter] to continue")
  
  # Assign value
  mothOccurrences$batRichness[i] <- value
}


# Export ------------------------------------------------------------------

fwrite(mothOccurrences,file.path(wd$out, paste0("batRichnessAtSite-", Sys.Date(), ".csv")), row.names = F)
