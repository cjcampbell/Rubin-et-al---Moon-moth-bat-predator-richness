source("R/00_setup.R")

# Find continuous surface files.
whichFiles <- list.files(wd$out, pattern = "SDM.tif$", full.names = T)
(sppNum <- length(whichFiles))

# Filter to target species.
source('/srv/duo/mbelitz/GeoBatDiv/R/02b_targetSpecies.R')
whichFiles <- basename(whichFiles) %>% 
  gsub("_SDM.tif", "", .) %>% 
  gsub("\\.", " ", .) %>% 
  {whichFiles[. %in% targetSpecies$canonical]}

# Note which species are included.
basename(whichFiles) %>% 
  gsub("_SDM.tif", "", .) %>% 
  gsub("\\.", " ", .) %>% 
  data.frame() %>% 
  fwrite(file = file.path(wd$out, "includedBatSpecies.csv"), row.names = F)

# Get target extent.
mod_vars_5k <- terra::rast(file.path(wd$data, "mod_vars_5k.tif"))
myextent <- terra::ext(mod_vars_5k)
rm(mod_vars_5k)

# Make temp dir.
mytempdir <- tempdir()
dir.exists(mytempdir)
list.files(mytempdir, full.names = T) %>% unlink() # Remove any files in tempdir (e.g., from test runs).

# Extend each tif, add to tempdir.
for(x in whichFiles) {
  myrast <- rast(x)
  terra::extend(myrast, myextent, filename = tempfile(tmpdir = mytempdir, fileext = ".tif"))
  print(paste("Done with", which(whichFiles == x), "of", length(whichFiles)))
}

# Combine.
mySDM <- rast( list.files(mytempdir, full.names = T, pattern = ".tif$") )
sum(mySDM, na.rm = T, filename = file.path(wd$out, "combinedContinuousSurface.tif"), overwrite = T)

# Examine
mySurface <- rast(file.path(wd$out, "combinedContinuousSurface.tif"))
plot(mySurface)
