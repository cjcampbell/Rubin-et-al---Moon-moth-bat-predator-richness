
source('R/06_0_prepToRunSDMPipeline.R')

############################ Execute the pipeline #############################


#### Step 2: Pass pipeline through the species list
#for(i in sample(1:length(spp_list), replace = F)) {
for(i in 1:length(targetList)) {
  
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
      names(mymod_vars) %in% c(
          "bio_1", 
          "bio_12", 
          "bio_13",
          "bio_14", 
          "bio_15", 
          "bio_16",
          "bio_17",
          "bio_2",
          "bio_4", 
          "bio_5", 
          "bio_6",
          #"bio_8",
          #"bio_9",
          "elevation", 
          "roughness",
          "tri", 
          "treeCover" 
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
