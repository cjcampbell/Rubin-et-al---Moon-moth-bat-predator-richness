
# List files.
list.files(wd$occs, recursive = T, full.names = T)

# Specify columns to read.
targetCols <-  c(
  "id", "institutionCode", "collectionCode", "basisOfRecord", "occurrenceID",
  "order", "family", "genus", "species", "specificEpithet", "scientificName", 
  "identifiedBy", "eventDate", "year", "month", "day",  "verbatimEventDate",
  "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", "locality", 
  "county", "coordinateUncertaintyInMeters", "source"
  )

# Load and tidy GBIF occurrences.
dat_gbif0 <- fread(file.path(wd$occs, "GBIF_download_20220907/0006282-220831081235567.csv")) 
dat_gbif  <- dat_gbif0 %>% 
  dplyr::mutate(source = "gbif") %>% 
  dplyr::rename(id = gbifID) %>% 
  dplyr::select( any_of(targetCols) ) %>% 
  distinct(id, .keep_all = T) 

# Load and tidy iDigBio occurrences.
dat_idb0 <- fread(file.path(wd$occs, "iDigBio_download_20220907/occurrence.csv"))
removePrefixes <- function(x) { stringr::word(x,-1, sep = ":") }
dat_idb  <- dat_idb0 %>% 
  select( "coreid", starts_with("dwc:")) %>% 
  dplyr::rename(id = coreid) %>% 
  dplyr::rename_all(removePrefixes) %>% 
  dplyr::mutate(source = "idigbio") %>% 
  dplyr::select( any_of(targetCols) ) %>% 
  dplyr::mutate_at(vars("id"), as.character())

# Combine occurrence records.
## Change classes as needed for successful binding.
dat_gbif$id <- as.character(dat_gbif$id)
dat_gbif$eventDate <- as.character(dat_gbif$eventDate)
## Bind rows
tdf <- rbind(dat_gbif, dat_idb, fill = TRUE)


# Remove obvious duplications.
tdf2 <- dplyr::distinct(tdf, occurrenceID, .keep_all = T)
tdf3 <- dplyr::distinct(tdf2, id, .keep_all = T)
tdf3 <- filter(tdf3, order == "Chiroptera")

fwrite(tdf3, file.path(wd$bin, "occurrences_tidy.csv"), row.names = F)
