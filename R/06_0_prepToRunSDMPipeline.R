
source("R/00_setup.R")
# Set up pipeline ---------------------------------------------------------

## Step 1: Load in the necessary libraries ----
invisible({
  library(rnaturalearth)
  library(rgeos)
  library(dismo)
  library(ENMeval)
  library(dplyr)
  library(stringr)
  library(data.table)
})

## Step 2: Load in the pipeline scripts ----

invisible({
  list.files(wd$fun, full.names = T) %>% 
    lapply(source)
})

## Step 3: Prepare the spatial data ----

# Load predictors
mod_vars_5k <- terra::rast(file.path(wd$data, "mod_vars_5k.tif"))
names(mod_vars_5k) <- gsub("wc2.1_30s_", "", names(mod_vars_5k))
names(mod_vars_5k)[names(mod_vars_5k)=="file29dedd1476c9b5_"] <- "treeCover"
names(mod_vars_5k)[names(mod_vars_5k)=="elevation_1KMmd_GMTEDmd"] <- "elevation"
names(mod_vars_5k)[names(mod_vars_5k)=="roughness_1KMmd_GMTEDmd"] <- "roughness"
names(mod_vars_5k)[names(mod_vars_5k)=="tri_1KMmd_GMTEDmd"] <- "tri"

### 3C. load occurrence records
occs <- fread(file.path(wd$data, "occurrences_clean.csv"))

### 3D. Read in basemap for visualizing
world <- ne_countries(scale = "medium", returnclass = "sf")

### 3E. Choose projection and project data if necessary
study_proj <- proj4_eqe
world <- st_transform(world, crs = study_proj)

occs_sf <- st_as_sf(occs, 
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = proj4_wgs) 
occs_sf <- st_transform(occs_sf, crs = study_proj)
occs <- occs %>% 
  mutate(x = st_coordinates(occs_sf)[,1], y = st_coordinates(occs_sf)[,2])


rerunPipeline <- FALSE

#### Step 1: Set up species list
source("R/02b_targetSpecies.R")
spp_list <- sort(targetSpecies$canonical)
targetList <- occs %>% 
  dplyr::filter(canonical %in% spp_list) %>% 
  dplyr::select(canonical, x, y) %>% 
  distinct %>% 
  count(canonical) %>% 
  dplyr::filter(n >= 3) %>% 
  #arrange(n) %>% 
  dplyr::select(canonical) %>% 
  unlist %>% 
  as.vector()
