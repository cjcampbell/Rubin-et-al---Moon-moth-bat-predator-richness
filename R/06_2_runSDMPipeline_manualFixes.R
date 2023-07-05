
source('R/06_0_prepToRunSDMPipeline.R')
# Check which species don't have completed SDM ----------------------------

didntRun <- lapply(make.names(targetList), function(x) {
  length( list.files(wd$out, pattern = x) ) == 0
})

checkEm <- targetList[which(unlist(didntRun))]
occs %>% 
  dplyr::filter(canonical %in% checkEm) %>% 
  dplyr::select(canonical, x, y) %>% 
  distinct %>% 
  count(canonical) %>% 
  arrange(desc(n))
length(checkEm)

## Fix alpha hull -----
for(i in which(targetList %in% c("Rhinolophus damarensis", "Rhinolophus robinsoni"))) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    ### Bug fix for specific species ###
    species_df2 <- species_df %>% 
      dplyr::mutate(dlat = floor(decimalLatitude), dlon = floor(decimalLongitude)) %>% 
      group_by(dlat, dlon) %>% 
      slice_sample(n=1)
      
    # Define accessible area.
    aa_shp <- define_accessibleArea(
      species_df = species_df2, minBuff = 200e3,
      buff_prop = 0.80, projCRS = study_proj,
      initialAlpha = 10
    )
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y
    
    ### 2B. Perform the rarefaction
    area_sqkm <- as.numeric(sf::st_area(aa_shp))/10e3
    species_df <- thinPoints(
      spp_df = species_df, 
      area_sqkm = area_sqkm, 
      bio = mymod_vars[[2]], 
      method = "simple",
      simpleMult = 5
    )
    
    # Run the rest of the pipeline. 
    
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    
    whichLayers <- which(
      names(mymod_vars) %in% c("bio_1", "bio_12", "bio_13","bio_14", "bio_15", "bio_16","bio_17","bio_2","bio_4", "bio_5", "bio_6","elevation", "roughness","tri", "treeCover" 
      )
    )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    #### Step 3: Test and fine-tune the model
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    max_model <- maxent(x = vars_selected, p = coordinates(species_df), progress = "text", silent = TRUE) 
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    bw <- make.names(targetList[i])
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df2 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df2 <- dplyr::rename(species_df2, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df2)
    
    ### 4C. Visualize the model
    ## Load in the rasters
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df2, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}



# Change method of getting accessible area -----

for(i in which(targetList ==  "Rhinolophus mossambicus")) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    # Define accessible area using sf. convex hull then buffer if I can't alpha hull, then buffer...
    aa_shp <- species_df %>% 
      st_as_sf(coords = c("x", "y"), crs = study_proj) %>% 
      st_union() %>% 
      st_convex_hull() %>% 
      st_buffer(200e3)
    st_crs(aa_shp) <- study_proj
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y
    
    ### 2B. Perform the rarefaction
    area_sqkm <- as.numeric(sf::st_area(aa_shp))/10e3
    species_df <- thinPoints(
      spp_df = species_df,
      area_sqkm = area_sqkm,
      bio = mymod_vars[[2]],
      method = "simple",
      simpleMult = 5
    )
    
    # Run the rest of the pipeline. 
    
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    
    whichLayers <- which(  names(mymod_vars) %in% c("bio_1", "bio_12", "bio_13","bio_14", "bio_15", "bio_16","bio_17","bio_2","bio_4", "bio_5", "bio_6","elevation", "roughness","tri", "treeCover" ) )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    #### Step 3: Test and fine-tune the model
    
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    max_model <- maxent(x = vars_selected, p = coordinates(species_df), progress = "text", silent = TRUE) 
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    bw <- make.names(targetList[i])
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df2 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df2 <- dplyr::rename(species_df2, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df2)
    
    ### 4C. Visualize the model
    ## Load in the rasters
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df2, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}


## Skip thinning, filter out NA-associated pts, AND specify bg pts manually -----

# Several weird examples of well-sampled species on  islands. Going to take some manual manipulation....
for(i in which(targetList %in% c("Nyctalus azoreum", "Lasiurus semotus", "Plecotus teneriffae", "Plecotus kolombatovici"))) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    # Define accessible area.
    aa_shp <- define_accessibleArea(
      species_df = species_df, minBuff = 200e3,
      buff_prop = 0.80, projCRS = study_proj
    )
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    coordinates(species_df) <- ~ x + y
    ## Filter out points with NA's
    vals <- terra::extract(mymod_vars, vect(species_df))
    species_df <- species_df[!is.na(vals$bio_1),]
    
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    
    whichLayers <- which( names(mymod_vars) %in% c( "bio_1",  "bio_12",  "bio_13", "bio_14",  "bio_15",  "bio_16", "bio_17", "bio_2", "bio_4",  "bio_5",  "bio_6", "elevation",  "roughness", "tri",  "treeCover"  ) )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    #### Step 3: Test and fine-tune the model
    
    back_pts <- randomPoints(vars_selected[[1]], n = 20) %>% 
      data.frame()

    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    max_model <- maxent(x = vars_selected, p = coordinates(species_df), a = back_pts, progress = "text", silent = TRUE) 
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      bg.coords = back_pts,
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    bw <- make.names(targetList[i])
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df2 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df2 <- dplyr::rename(species_df2, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    
    ### Manually save :(
    ENMeval_output = eval1
    AUCmin = 0.7
    spp = bw
    occ_df = species_df2
    
    spp <- stringr::str_replace(string = spp, pattern = " ", replacement = "_")
    
    ### Find the best model
    bestmod <- ENMeval_output@results %>% dplyr::filter(delta.AICc == 0) %>% slice(1)
    if(bestmod$auc.train >= AUCmin & bestmod$auc.val.avg >= AUCmin){
      bestmod <- bestmod[1,]
    } else{
      bestmod <- ENMeval_output@results %>% 
        dplyr::filter(auc.train == max(auc.train) | auc.val.avg  == max(auc.val.avg)) 
      bestmod <- bestmod[1,]
    }
    write.csv(bestmod, file = file.path(resultDir, paste0( spp, "_bestModel.csv")), row.names = F)
    
    ### Build the best model
    maxent_args <- as.character(bestmod$tune.args)
    r_best <- ENMeval_output@predictions[[maxent_args]]
    raster::writeRaster(x = r_best, filename = file.path(resultDir, paste0(spp, "_SDM.tif")), overwrite = TRUE)
    
    ### Build a presence-absence model
    # Manually classify to present in islands...
    r_best <- terra::rast(r_best)
    back_pts <- as.data.frame(spatSample(r_best, size = nrow(occ_df)*10000000, xy = T)) %>% 
      na.omit() %>% 
      #slice_sample(n = nrow(occ_df)) %>% 
      dplyr::select(x, y)
    spp_df <- occ_df %>% 
      dplyr::select(x,y)
    lpt <- terra::extract(x = r_best, y = spp_df) %>% 
      dplyr::rename(cloglog = 2) %>% dplyr::select(cloglog)
    
    pa <- c(0, 1e-10, 0, 1e-10, 1, 1)
    pa_mat <- matrix(pa, ncol = 3, byrow = TRUE)
    r_best_pa <- terra::classify(r_best, pa_mat)
    terra::writeRaster(x = r_best_pa, 
                       filename = file.path(resultDir, paste0(spp, "_SDM_PA.tif")), overwrite = TRUE)
    # Manually set value to indicate manual manipulation...
    varimp <- ENMeval_output@variable.importance[[maxent_args]] %>% 
      dplyr::mutate(percent.contribution = "999", permutation.importance = "999")
    write.csv(x = varimp,               
              file = file.path(resultDir, paste0(spp, "_variableImportance.csv")),
              row.names = F)
    
    ### 4C. Visualize the model
    ## Load in the rasters
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df2, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}

## New accessible area method using convex hull ----

for(i in which(targetList == "Mops jobensis")) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")

    # Define accessible area.
    aa_shp <- species_df %>% 
      st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84") %>% 
      st_union() %>% 
      st_convex_hull() %>% 
      st_transform(study_proj) %>% 
      st_buffer(dist = 1000e3)
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y

    # Don't thin points. 
    
    ### 2B. Perform the rarefaction
    area_sqkm <- as.numeric(sf::st_area(aa_shp))/10e3
    species_df <- thinPoints(
      spp_df = species_df, 
      area_sqkm = area_sqkm, 
      bio = mymod_vars[[2]], 
      method = "simple",
      simpleMult = 5
    )
    
    
    # Run the rest of the pipeline. 
    
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    
    whichLayers <- which(
      names(mymod_vars) %in% c( "bio_1",  "bio_12",  "bio_13", "bio_14",  "bio_15",  "bio_16", "bio_17", "bio_2", "bio_4",  "bio_5",  "bio_6", "elevation",  "roughness", "tri",  "treeCover" 
      )
    )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    #### Step 3: Test and fine-tune the model
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    max_model <- maxent(x = vars_selected, p = coordinates(species_df), progress = "text", silent = TRUE) 
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    bw <- make.names(targetList[i])
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df2 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df2 <- dplyr::rename(species_df2, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df2)
    
    ### 4C. Visualize the model
    ## Load in the rasters
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df2, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}

# Another island species fix ----------------------------------------------
# For some less-sampled island species, it makes sense to move points to nearest cell with a non-NA value...
# Only do this after visualizing points! :)

for(i in which(targetList %in% c( "Rhinolophus maendeleo", "Scotophilus marovaza", "Rhinolophus canuti", "Hipposideros kunzi", "Hipposideros papua","Hipposideros pelingensis", "Hipposideros tephrus","Hipposideros kunzi","Rhinolophus keyensis", "Rhinolophus yunanensis", "Eptesicus magellanicus", "Plecotus kolombatovici", "Mops sarasinorum", "Molossus alvarezi", "Scotophilus borbonicus"))) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    bw <- make.names(targetList[i])
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    # Define accessible area.
    aa_shp <- define_accessibleArea(
      species_df = species_df, minBuff = 200e3,
      buff_prop = 0.80, projCRS = study_proj
    )
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y
    
    ### 2B. No rarefaction.
    
    # Run the rest of the pipeline. 
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    whichLayers <- which(names(mymod_vars) %in% c("bio_1", "bio_12", "bio_13","bio_14", "bio_15", "bio_16","bio_17","bio_2","bio_4", "bio_5", "bio_6","elevation", "roughness","tri", "treeCover" ) )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    ## Move coordinate points to nearest non-NA cell. ---------
    rast_df <- vars_selected[[1]] %>% as.data.frame(xy=T) %>% na.omit()
    species_df2 <- data.frame(species_df)
    mydf <- species_df %>% st_as_sf()
    
    ggplot() + 
      geom_tile(rast_df, mapping = aes(x=x,y=y, fill = bio_1, color = bio_1)) +
      geom_sf(mydf, mapping = aes()) +
      theme_bw()
    
    st_crs(mydf) <- proj4_eqe
    for(q in 1:nrow(mydf))  {
      whichMin <- rast_df %>% 
        st_as_sf(coords = c("x", "y"), crs = proj4_eqe) %>% 
        st_distance(mydf[q,]) %>% 
        which.min
      species_df2[q,c("x2", "y2")] <- rast_df[whichMin,c("x", "y")]
    }
    coordinates(species_df2) <- ~ x2 + y2
    #terra::extract(vars_selected, species_df2)
    
    #### Step 3: Test and fine-tune the model
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    try({
      max_model <- maxent(x = vars_selected, p = coordinates(species_df2), progress = "text", silent = TRUE) 
    })

    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df2), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df3 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df3 <- dplyr::rename(species_df3, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df3)
    
    ### 4C. Visualize the model
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df3, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}

### Move THEN  increase accessible area to help background points fit -----
for(i in which(targetList %in% c("Hipposideros turpis"))) {
  
  if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    bw <- make.names(targetList[i])
    species_df = dplyr::filter(occs, canonical == targetList[i])
    stopifnot(nrow(species_df) >= 3)
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    # Define accessible area.
    aa_shp <- define_accessibleArea(
      species_df = species_df, minBuff = 500e3,
      buff_prop = 0.80, projCRS = study_proj
    )
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y
    
    ### 2B. No rarefaction.
    
    # Run the rest of the pipeline. 
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    whichLayers <- which(names(mymod_vars) %in% c("bio_1", "bio_12", "bio_13","bio_14", "bio_15", "bio_16","bio_17","bio_2","bio_4", "bio_5", "bio_6","elevation", "roughness","tri", "treeCover" ) )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    ## Move coordinate points to nearest non-NA cell. ---------
    rast_df <- vars_selected[[1]] %>% as.data.frame(xy=T) %>% na.omit()
    species_df2 <- data.frame(species_df)
    mydf <- species_df %>% st_as_sf()
    
    ggplot() + 
      geom_tile(rast_df, mapping = aes(x=x,y=y, fill = bio_1, color = bio_1)) +
      geom_sf(mydf, mapping = aes()) +
      theme_bw()
    
    st_crs(mydf) <- proj4_eqe
    for(q in 1:nrow(mydf))  {
      whichMin <- rast_df %>% 
        st_as_sf(coords = c("x", "y"), crs = proj4_eqe) %>% 
        st_distance(mydf[q,]) %>% 
        which.min
      species_df2[q,c("x2", "y2")] <- rast_df[whichMin,c("x", "y")]
    }
    coordinates(species_df2) <- ~ x2 + y2
    #terra::extract(vars_selected, species_df2)
    
    #### Step 3: Test and fine-tune the model
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    try({
      max_model <- maxent(x = vars_selected, p = coordinates(species_df2), progress = "text", silent = TRUE) 
    })
    
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df2), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df3 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df3 <- dplyr::rename(species_df3, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df3)
    
    ### 4C. Visualize the model
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df3, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}


## Helper: plot points for remaining species

for(i in which(targetList %in% checkEm)) {
  
  bw <- make.names(targetList[i])
  species_df = dplyr::filter(occs, canonical == targetList[i])
  stopifnot(nrow(species_df) >= 3)
  
  # Clip all potential predictor variables based on all cleaned occurrence
  # records from this species.
  
  # Define accessible area.
  aa_shp <- define_accessibleArea(
    species_df = species_df, minBuff = 200e3,
    buff_prop = 0.80, projCRS = study_proj
  )
  
  # Clip environmental variable layers to the defined accessible area
  mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
  
  # Thin points based on the accessible area.
  ### 2A. Prepare the coordinates for rarefaction
  coordinates(species_df) <- ~ x + y
  
  ### 2B. No rarefaction.
  
  # Run the rest of the pipeline.
  whichLayers <- which(names(mymod_vars) %in% c("bio_1", "bio_12", "bio_13","bio_14", "bio_15", "bio_16","bio_17","bio_2","bio_4", "bio_5", "bio_6","elevation", "roughness","tri", "treeCover" ) )
  vars_selected <- raster::stack(mymod_vars[[whichLayers]])
  
  ## Move coordinate points to nearest non-NA cell. ---------
  rast_df <- vars_selected[[1]] %>% as.data.frame(xy=T) %>% na.omit()
  species_df2 <- data.frame(species_df)
  mydf <- species_df %>% st_as_sf()
  
(  p <- ggplot() + 
    geom_tile(rast_df, mapping = aes(x=x,y=y, fill = bio_1, color = bio_1)) +
    geom_sf(mydf, mapping = aes()) +
    theme_bw() + 
    ggtitle(targetList[i]))
  
  Sys.sleep(1)
  
  readline(prompt="Press [enter] to continue")
}




# Rerun some species with more aggressive filtering -----------------------



for(i in which(targetList %in% c("Eptesicus nilssonii","Eptesicus serotinus", "Nyctalus lasiopterus", "Nyctalus leisleri", "Rhinolophus euryale", "Rhinolophus ferrumequinum", "Tadarida teniotis", "Vespertilio murinus")) ) {
  
  #if(length(list.files(wd$out, pattern = make.names(targetList[i]))) >= 6 & !rerunPipeline) next
  
  try({
    species_df = dplyr::filter(occs, canonical == targetList[i])
 
    # Filter down european points AGGRESSIVELY.
    euro <- species_df %>% 
      dplyr::filter(decimalLongitude < 40) %>% 
      # Filter out near-duplicates
      dplyr::mutate(
        ddx = round(x,-5),
        ddy = round(y,-5)) %>% 
      group_by(ddx, ddy) %>% 
      slice(1) %>% 
      ungroup %>%
      # And cap at 100 pts from europe.
      slice_sample(n=100)
    species_df <- bind_rows(
      dplyr::filter(species_df, decimalLongitude >= 40),
      euro
    )
    
    # Filter down korean points if needed.... (Rhinolophus ferrumequinum)
    korean_big <- species_df %>% 
      dplyr::filter(
        (decimalLatitude >= 33 & decimalLatitude <= 38.6 & decimalLongitude > 124 & decimalLongitude < 130)
      )
    if(nrow(korean_big) > 0) {
      korean <- korean_big %>% 
        # Filter out near-duplicates
        dplyr::mutate(
          ddx = round(x,-5),
          ddy = round(y,-5)) %>% 
        group_by(ddx, ddy) %>% 
        slice(1) %>% 
        ungroup %>% 
        slice_sample(n = 20)
      species_df <- species_df %>% 
        dplyr::filter(!(id %in% korean_big$id)) %>% 
        bind_rows(korean)
    }
    
    
    print("=========================================================")
    print(paste("============== Working on", targetList[i], "=============="))  
    print("=========================================================")
    
    
    # Clip all potential predictor variables based on all cleaned occurrence
    # records from this species.
    
    # Define accessible area.
    aa_shp <- define_accessibleArea(
      species_df = species_df, minBuff = 200e3,
      buff_prop = 0.80, projCRS = study_proj
    )
    
    # Clip environmental variable layers to the defined accessible area
    mymod_vars <- clip_variableLayers(rstack = mod_vars_5k, accessibleArea = aa_shp)
    
    # Thin points based on the accessible area.
    ### 2A. Prepare the coordinates for rarefaction
    coordinates(species_df) <- ~ x + y
    
    ### 2B. Perform the rarefaction
    # Rarify even more
    area_sqkm <- as.numeric(sf::st_area(aa_shp))/10e3
    species_df <- thinPoints(
      spp_df = species_df, 
      area_sqkm = area_sqkm, 
      bio = mymod_vars[[2]], 
      method = "simple",
      simpleMult = 5
    )
    
    # Run the rest of the pipeline. 
    
    print("~~~~~~~~~~ Working on model fitting ~~~~~~~~~~")
    
    whichLayers <- which( names(mymod_vars) %in% c( "bio_1",  "bio_12",  "bio_13", "bio_14",  "bio_15",  "bio_16", "bio_17", "bio_2", "bio_4",  "bio_5",  "bio_6", "elevation",  "roughness", "tri",  "treeCover"  ) )
    vars_selected <- raster::stack(mymod_vars[[whichLayers]])
    
    #### Step 3: Test and fine-tune the model
    
    ### 3A. Select top performing variables to reduce colinearity
    ## First run a test model
    print(">>>>>>>>> Running Maxent test model <<<<<<<<<")
    max_model <- maxent(x = vars_selected, p = coordinates(species_df), progress = "text", silent = TRUE) 
    ## Using the test model, iteratively test the colinearity of variables, removing highly colinear ones one at a time 
    print(">>>>>>>>> Selecting top SDM variables <<<<<<<<<")
    predictors <- select_sdmVariables(pred_vars = vars_selected, maxent_mod = max_model, maxVIF = 5)
    
    ### 3B. Evaluate various tuning parameters of the model for fine tuning to create the best performing model with best tuned parameters
    print(">>>>>>>>> Evaluating tuning variables in model <<<<<<<<<")
    
    eval1 <- ENMeval::ENMevaluate(
      occ = coordinates(species_df), env = predictors,
      method = "block", RMvalues = c(0.5, 1, 2, 3, 4),
      fc= c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
      parallel = TRUE, numCores = 20, algorithm = 'maxent.jar',
      quiet = TRUE, updateProgress = FALSE
    )
    
    #### Step 4: Create the output for the best performing model
    
    ### 4A. Prepare the output path and coordinates
    bw <- make.names(targetList[i])
    
    resultDir = wd$out
    if (!dir.exists(resultDir)) { dir.create(resultDir) }
    # Return coordinates to a data frame
    species_df2 <- as.data.frame(species_df)
    if("coords.x1" %in% names(species_df2)) {
      species_df2 <- dplyr::rename(species_df2, x = coords.x1, y = coords.x2)
    }
    
    ### 4B. Output the model
    ## Output the model evaluation
    save(eval1, file = file.path(resultDir, paste( bw, "_ENMeval", ".RData", sep = "")))
    ## Output the best performing model, including both the actual model and the presence-absence model
    print(">>>>>>>>> Saving best model <<<<<<<<<")
    
    save_SDM_results(ENMeval_output = eval1, AUCmin = 0.7, resultDir = resultDir,
                     spp = bw, occ_df = species_df2)
    
    ### 4C. Visualize the model
    ## Load in the rasters
    r    <- raster(file.path(resultDir, paste0( bw,"_SDM.tif")))
    r_pa <- raster(file.path(resultDir, paste0( bw,"_SDM_PA.tif")))
    
    create_sdmFigure(
      spp = bw, r = r, r_pa = r_pa, occ_df = species_df2, world = world
    )
    
    print(paste0("!!!!!! Model for ", gsub("_", " ", bw), " complete !!!!!"))
    
  })
  
  # Clear tempdir
  options(java.parameters = "-Xmx8000m")
  unixtools::set.tempdir("tmp/")
  tmp_dir <- unixtools::set.tempdir("tmp/")
  files <- list.files(tmp_dir, full.names = T,  all.files = T, recursive = T)
  file.remove(files)
  gc()
  
}
