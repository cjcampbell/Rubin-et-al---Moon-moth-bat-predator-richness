
source("R/00_setup.R")
source("R/02b_targetSpecies.R")
library(tidyr)

# Find number of records per species
occs <- fread(file.path(wd$data, "occurrences_clean.csv"))

# Count records and specify if at least 3 exist.
mdf <- occs %>% 
  dplyr::filter(canonical %in% targetSpecies$canonical) %>% 
  count(canonical, name = "n_records") %>% 
  full_join(targetSpecies) %>% 
  replace_na(list(n_records = 0)) %>% 
  dplyr::mutate(atLeast3records = case_when(n_records >= 3 ~ 1, TRUE ~ 0)) %>% 
  ungroup

# Check if SDM run for each species.
sdm_run <- list.files(wd$out, pattern = "variableImportance.csv") %>% 
  gsub("_variableImportance.csv", "", .) %>% 
  gsub("[.]", " ", .) %>% 
  {.[. %in% targetSpecies$canonical]}

# Check if SDM failed completely.
varImp <- list.files(wd$out, pattern = "variableImportance.csv", full.names = T) %>% 
  lapply(function(x) {
    y <- fread(x)
    y$canonical <- basename(gsub( "[.]"," ",gsub( "_variableImportance.csv","", x)))
    return(y)
  }) %>% 
  bind_rows()

sdm_failed <- varImp %>% 
  group_by(canonical) %>% 
  summarise(varImp = sum(permutation.importance)) %>% 
  dplyr::filter(varImp == 0) %>% 
  dplyr::select(canonical) %>% 
  unlist

mdf2 <- mdf %>% 
  as.data.frame() %>% 
  dplyr::mutate(
    SDMrun = canonical %in% sdm_run,
    SDMnotFailed = !(canonical %in% sdm_failed)
  ) %>% 
  dplyr::mutate_if(is.logical, as.numeric)

# Which species don't have SDM but should?
mdf2 %>% 
  dplyr::filter(atLeast3records == 1 & SDMrun == 0) %>% 
  arrange(desc(n_records)) %>% 
  dplyr::select(canonical, starts_with("SDM"),contains("records"))

# Which SDM failed?
mdf2 %>% 
  dplyr::filter(SDMrun == 1 & SDMnotFailed == 0) %>% 
  arrange(desc(n_records)) %>% 
  dplyr::select(canonical, starts_with("SDM"),contains("records"))


# Explore variable importance ---------------------------------------------

varimp <- list.files(wd$out, pattern = "variableImportance.csv", full.names = T) %>% 
  lapply(function(x) {
    fread(x) %>% 
      dplyr::mutate(order = row_number())
  }) %>% 
  bind_rows
varimp %>% 
  group_by(variable) %>% 
  dplyr::summarise(sum = sum(percent.contribution)) %>% 
  arrange(desc(sum))



# Add manual checks for SDM performance -------------------------------------------------------
rerunCheck <- FALSE

# Run for some % of models.
set.seed(42)
while(length(list.files(wd$bin, pattern = "_SDMoutputCheck.rds")) < round(length(sdm_run)*1) ) {
  
  i <- sample(1:length(sdm_run), size = 1)
  
  spname <- gsub(" ", ".", sdm_run[i])
  
  if(!file.exists(file.path(wd$bin, paste0(spname, "_SDMoutputCheck.rds"))) | isTRUE(rerunCheck) ) {
    print( sdm_run[i])
    
    # Does file exist
    filepath <-  list.files(
      wd$out, pattern = paste0(spname, "_map.png"),
      recursive = T, full.names = T
    )
    if(length(filepath) == 0) { print("File doesn't exist; skipping"); next }
    
    # Open map
    filepath %>%
      browseURL()
    
    passfail <- readline(prompt= paste(
      "
    1 - Looks good! (pass)
    2 - Looks questionable (check)
    3 - Looks bad (fail)
    4 - Skip for now
    "))
    
    if(passfail %in% 1:3){
      data.frame(species = spname, result = passfail) %>%
        saveRDS(file = file.path(wd$bin, paste0(spname, "_SDMoutputCheck.rds")))
    }
  }
}



# Combine and check outputs of manual check. ----------------------------------------------

checkEmResults <- list.files(wd$bin, pattern = "_SDMoutputCheck.rds", full.names = T) %>% 
  lapply(readRDS) %>% 
  bind_rows()

checkEmResults %>% 
  dplyr::filter(result != 1)

checkTheseSpecies <- checkEmResults %>% 
  dplyr::filter(result != 1) %>% 
  dplyr::select(species) %>% 
  unlist %>% 
  as.character()

# Manual check reveals a few problems:
# - taxonomic splitting in Lasiurus species
# - Eurasian species sampled in Europe have waaaaay more points and habitat suitability than in asia. Try to filter more stringently.
# - a couple just look wacky due to manual bug fixes that didn't work very well...

#       Eptesicus.nilssonii      2 - resolved, filter better
#       Eptesicus.serotinus      2 - resolved, filter better
#        Hipposideros.papua      3 - continuous surface is fine, PA one is messed up (but we aren't using.)
#       Hipposideros.turpis      3 - continuous surface is fine, PA one is messed up (but we aren't using.)
#     Lasiurus.blossevillii      3 - resolved, was a taxon split / pt filtering issue
#         Lasiurus.borealis      3 - resolved, was a taxon split / pt filtering issue
#         Lasiurus.cinereus      3 - resolved, was a taxon split / pt filtering issue
#         Molossus.alvarezi      2 - resolved, jk it looks fine?
#      Nyctalus.lasiopterus      2 - resolved, filter better
#         Nyctalus.leisleri      2 - resolved, filter better
#       Rhinolophus.euryale      2 - resolved, filter better
# Rhinolophus.ferrumequinum      2 - resolved, filter better
#         Tadarida.teniotis      2 - resolved, filter better
#       Vespertilio.murinus      2 - resolved, filter better

# Rerun checks for focal species
for(spname in checkTheseSpecies) {
  print(spname)
  # Does file exist
  filepath <-  list.files(
    wd$out, pattern = paste0(spname, "_map.png"),
    recursive = T, full.names = T
  )
  if(length(filepath) == 0) { print("File doesn't exist; skipping"); next }
  
  # Open map
  filepath %>%
    browseURL()
  
  passfail <- readline(prompt= paste(
    "
    1 - Looks good! (pass)
    2 - Looks questionable (check)
    3 - Looks bad (fail)
    4 - Skip for now
    "))
  
  if(passfail %in% 1:3){
    data.frame(species = spname, result = passfail) %>%
      saveRDS(file = file.path(wd$bin, paste0(spname, "_SDMoutputCheck2.rds")))
  }
}

checkEm2Results <- list.files(wd$bin, pattern = "_SDMoutputCheck2.rds", full.names = T) %>% 
  lapply(readRDS) %>% 
  bind_rows()










