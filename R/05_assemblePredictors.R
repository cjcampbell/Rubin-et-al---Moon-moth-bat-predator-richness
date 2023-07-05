library(terra)
library(magrittr)

# Download MODIS data -----------------------------------------------------
library(luna)
MODIStilePath <- file.path(wd$data, "MOD44B")
if(!dir.exists(MODIStilePath)) dir.create(MODIStilePath)

# Download MODIS VCF data.
luna::getModis(
  product = "MOD44B", 
  start_date = "2020-09-01", end_date = "2020-10-01",
  aoi = c(-180,180,-90,90), 
  path = MODIStilePath, 
  version ="061", download = T, 
  username = earthdata$usr, password = earthdata$pwd
  )

# Write rasters to a tempdir with just Percent Tree Cover layer, after projecting.
# This will take awhile!
# NATIVE RESOLUTION 4885.872

TreeCoverPath <- file.path(wd$bin, "TreeCover")
if(!dir.exists(TreeCoverPath)) dir.create(TreeCoverPath)
ff <- list.files(MODIStilePath, full.names = T, pattern = "\\.061")
#bettermc::mclapply(ff, mc.cores = 15, function(i) {
for(i in ff) {
  savename <-  file.path(TreeCoverPath, paste0(basename(i), ".tif"))
  if(!file.exists(savename)) {
    try({treecover <- rast(i)$Percent_Tree_Cover})
    if(exists("treecover")) {
      try({ terra::project(treecover, terra::crs(proj4_eqe), method = "bilinear", filename=savename) })
      writeLines(paste("Done with", round(100*length(list.files(TreeCoverPath))/length(ff), 1), "%." ))
    }
  }
}
#}); print(paste("Done with", length(list.files(TreeCoverPath, full.names = T, pattern = "\\.061"))))

# Create VRT file instead of mosaic-ing.
vrtfile <- file.path(paste0(tempfile(), "_.tif"))
v <- vrt(list.files(TreeCoverPath, full.names = T, pattern = "\\.061"), vrtfile, overwrite = T)
# Write to a tempfile tif
terra::writeRaster(v, filename = file.path(wd$preds, "TreeCover.tif"), overwrite  = T)
refRast <- rast(file.path(wd$preds, "TreeCover.tif"))
plot(refRast)


# Time invariate predictors -------------------------------------------------

## Bioclimatic ----
# Native resolution 629.2727m
for(i in c("_1.tif", "_2.tif", "_4.tif", "_5.tif", "_6.tif", "_8.tif", "_9.tif","_10.tif","_11.tif", "_12.tif", "_13.tif", "_14.tif", "_15.tif", "_16.tif", "_17.tif")) {
  if(!exists("refRast")) refRast <- rast(file.path(wd$preds, "TreeCover.tif"))
  myfilename <- file.path(wd$preds, paste0("bioclim", i))
  if(file.exists(myfilename)) next
  writeLines(paste("BIOCLIM -", i))
  a <- rast(list.files(
    "/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/wc2.1_30s_bio", 
    pattern = i, 
    full.names = T
  ) )
  myfile <- paste0(tempfile(), ".tif")
  terra::project(a, proj4_eqe, filename = myfile, overwrite = T)
  a <- rast(myfile)
  terra::resample(a, refRast, filename = myfilename, overwrite = T)
  unlink(myfile)
}

## Topographic ----
# Native resolutiion 629.2727m
for(i in c("elevation_1KMmd_GMTEDmd", "roughness_1KMmd_GMTEDmd", "tri_1KMmd_GMTEDmd")[1]) {
  myfilename <- file.path(wd$preds, paste0("topo_", i, ".tif"))
  if(file.exists(myfilename)) next
  a <- list.files("/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/topographicVars", pattern = ".tif$", full.names = T) %>% 
    terra::rast(lyr = i) 
  myfile <- tempfile()
  myfile <- paste0(tempfile(), ".tif")
  terra::project(a, proj4_eqe, filename = myfile, overwrite = T)
  a <- rast(myfile)
  terra::resample(a, refRast, filename = myfilename, overwrite = T)
  unlink(myfile)
}


# # Make season-specific predictors -----------------------------------------
# makeSeasonSpecificBioClim(
#   parentDir = "/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/wc2.1_30s_prec",
#   varName = "cumPrec", 
#   fxn = "sum", 
#   outDir = wd$preds, refRast = refRast
# )
# 
# makeSeasonSpecificBioClim(
#   parentDir = "/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/wc2.1_30s_tmax",
#   varName = "maxTmax", 
#   fxn = "max", 
#   outDir = wd$preds, refRast = refRast
# )  
# 
# makeSeasonSpecificBioClim(
#   parentDir = "/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/wc2.1_30s_tmin",
#   varName = "minTmin", 
#   fxn = "min", 
#   outDir = wd$preds, refRast = refRast
# )  
# 
# makeSeasonSpecificBioClim(
#   parentDir = "/srv/duo/mbelitz/NoAm_Bat_SDM/data/preds_setup/wc2.1_30s_tavg",
#   varName = "aveTavg", 
#   fxn = "mean", 
#   outDir = wd$preds, refRast = refRast
# )  

# Assemble ----------------------------------------------------------------

if(!exists("mod_vars")) {
  
  if(file.exists(file.path(wd$data, "mod_vars.tif")) ) {
    mod_vars <- terra::rast( file.path(wd$data, "mod_vars.tif") )
  } else {
    ### 3A. Load in model variables
    varlist <- list.files(wd$preds, full.names = TRUE)
    mod_vars <- rast(varlist)
    # Check names. Tidy as needed.
    if(any(duplicated(names(mod_vars_wgs)))) stop("Names must not be duplicated")
    #c("aveTavg_1", "aveTavg_2", "aveTavg_3", "aveTavg_4", "bio_1", "bio_12", "bio_13", "bio_14", "bio_15", "bio_16", "bio_17", "bio_2", "bio_4", "bio_5", "bio_6", "bio_8", "bio_9", "cumPrec_1", "cumPrec_2", "cumPrec_3", "cumPrec_4", "Lights", "Popdensity", "maxTmax_1", "maxTmax_2", "maxTmax_3", "maxTmax_4", "minTmin_1", "minTmin_2", "minTmin_3", "minTmin_4", "elevation", "roughness", "tri", "treeCover")
    
    #mod_vars <- terra::project(mod_vars_wgs, study_proj)
    terra::writeRaster(mod_vars, file = file.path(wd$data, "mod_vars.tif"))
    
    # Aggregate
    mod_vars_5k <- terra::aggregate(mod_vars, 5000/res(mod_vars)[1], filename = file.path(wd$data, "mod_vars_5k.tif"), overwrite = T)
    #terra::writeRaster(mod_vars_5k, file = file.path(wd$data, "mod_vars_5k.tif"))
  }
}



