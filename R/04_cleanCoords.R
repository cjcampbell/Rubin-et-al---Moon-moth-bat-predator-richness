
# Setup and functions -----------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(CoordinateCleaner)
library(rnaturalearth)
library(sf)
library(stringr)
library(dbplyr)
library(RSQLite)
library(ggpubr)

source("R/00_setup.R")
source("R/02b_targetSpecies.R")
runForAll <- FALSE # Run for all species or just target species?

targetCols <- c("canonical", "year", "month", "day", "eventDate","decimalLongitude",
                "decimalLatitude", "id")

# Functions ---------------------------------------------------------------

# Function to clean coordinates from occurrence records.
# Input is a dataframe containing occurrence records for one species only (of name "spp").
map_coord_clean <- function(dataframe, spp, out = wd$coords){

  # remove NAs n decimalLongitude & Latitude and make sure decimalLong & Lats are distinct
  cs <- dataframe %>%
    filter(decimalLongitude >= -180 & decimalLongitude <= 180) %>%
    filter(decimalLatitude >= -90 & decimalLatitude <= 90) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = T)

  cs2 <- clean_coordinates(cs, lon = "decimalLongitude", lat = "decimalLatitude",
                           species = "canonical",
                           tests = c("capitals", "centroids", "equal",
                                     "gbif", "institutions", "outliers"),
                           verbose = T) %>%
    dplyr::mutate(.summary = factor(.summary, levels = c(FALSE, TRUE)))

  mytab <- cs2 %>%
    tidyr::pivot_longer(cols =c(".val", ".equ", ".cap", ".cen", ".otl", ".gbf", ".inst") ) %>%
    group_by(name, value) %>%
    dplyr::summarise(n=n(), .groups = "drop_last") %>%
    dplyr::filter(value == FALSE | name == ".gbf")

  bw <- str_replace(spp, " ", "_")
  write.csv(x = cs2, file = file.path(out, paste0(bw, ".csv")))

  world <- ne_countries(scale = "medium", returnclass = "sf")
  #us <- ne_countries(continent = "North America", returnclass = "sf")

  lonRange <- (max(cs$decimalLongitude)- min(cs$decimalLongitude))/10
  latRange <- (max(cs$decimalLatitude)- min(cs$decimalLatitude))/10

  a0 <- ggplot() +
    geom_sf(world, mapping = aes(), fill = "grey95", color = "grey80")

  if(!exists("iucnRangemaps")) iucnRangemaps <- read_sf("/Users/cjcampbell/BigZaddyData/IUCN_BatRanges/redlist_species_data_d9b37aa6-89b6-4fe4-9ebc-8567d9a63a8c/data_0.shp") %>%
    st_simplify(10) %>%
    sf::st_transform(st_crs(world))

  if(spp %in% iucnRangemaps$BINOMIAL) {
    print("Rangemap is available")
    range <- iucnRangemaps %>% dplyr::filter( BINOMIAL == unlist(spp))
    a0 <- a0 +
      geom_sf(range, mapping = aes(), alpha = 0.8, color = NA, fill = "blue") +
          coord_sf(xlim = c(min(c(cs$decimalLongitude, sf::st_bbox(range)[1])) - lonRange, max(c(cs$decimalLongitude, sf::st_bbox(range)[3])) + lonRange),
                    ylim = c(min(c(cs$decimalLatitude, sf::st_bbox(range)[2]))- latRange, max(c(cs$decimalLatitude, sf::st_bbox(range)[4])) + latRange))

  } else {
    a0 <- a0 +
      coord_sf(xlim = c(min(cs$decimalLongitude) - lonRange, max(cs$decimalLongitude) + lonRange),
               ylim = c(min(cs$decimalLatitude) - latRange, max(cs$decimalLatitude) + latRange))

  }

  b <- a0 +
    geom_point(
      dplyr::filter(cs2, .summary == "TRUE"),
      mapping = aes(x = decimalLongitude, y = decimalLatitude),
      color = "black")  +
    geom_point(
      dplyr::filter(cs2, .summary == "FALSE"),
      mapping = aes(x = decimalLongitude, y = decimalLatitude),
      shape = 21, fill = "yellow", color = "#786213")  +
    geom_point(
      dplyr::filter(cs2, .otl == "FALSE"),
      mapping = aes(x = decimalLongitude, y = decimalLatitude),
      shape = 21, fill = "green", color = "#786213")  +
     labs(color = "Not flagged") +
    ggtitle("Cleaned") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  #cp <- ggarrange(a, b, ggtexttable(mytab, rows = NULL),  nrow = 3, ncol = 1, heights = c(2,2,1))
  cp <- ggarrange(b, ggtexttable(mytab, rows = NULL),  nrow = 2, ncol = 1, heights = c(2,1))

  ggsave(plot = cp,
         filename = file.path(out, paste0(bw, ".png")),
         width = 6, height = 8)
}

quickPlot <- function(cs) {
  if(!exists("world")) world <- ne_countries(scale = "medium", returnclass = "sf")
  lonRange <- (max(cs$decimalLongitude)- min(cs$decimalLongitude))/10
  latRange <- (max(cs$decimalLatitude)- min(cs$decimalLatitude))/10

  ggplot() +
    geom_sf(world, mapping = aes(), fill = "grey95", color = "grey80") +
    geom_point(cs, mapping = aes(x = decimalLongitude, y = decimalLatitude),
               color = "#18806d") +
    coord_sf(xlim = c(min(cs$decimalLongitude) - lonRange, max(cs$decimalLongitude) + lonRange),
             ylim = c(min(cs$decimalLatitude) - latRange, max(cs$decimalLatitude) + latRange)) +
    ggtitle("Uncleaned") +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
}


### function to pull uncleaned pts from coords/, revise, plot, and resave to coords2/----
fixPts <- function(myspp, expr, plot = FALSE) {
  spname <- gsub(" ", "_", myspp)
  print(spname)
  if(isTRUE(plot)) {
    list.files( wd$coords, pattern = paste0(spname, ".png"),
                recursive = T, full.names = T
    ) %>% browseURL()
  }

  myfile <- list.files(wd$coords2, pattern = paste0(gsub("_"," ", myspp), ".csv"), full.names = T)
  if(length(myfile) == 0) { myfile <- list.files(wd$coords2, pattern = paste0(myspp, ".csv"), full.names = T)}
  expr <- rlang::parse_expr(expr)
  out <- fread(myfile) %>%
    dplyr::mutate(id = as.numeric(id)) %>%
    dplyr::filter(!!expr)
  fwrite(out, file.path(wd$coords3, paste0(gsub("_"," ", myspp), ".csv") ))
  if(isTRUE(plot) & nrow(out) > 0) return(quickPlot(out))
  #unlink(myfile)
}

### function to identify points within a certain range of an IUCN rangemap ----
getIDsOfInRangePts <- function(myspp, dist = 200) {
  # Requires species name and distance in km from IUCN rangemap that is acceptable.
  # Returns id's of points within that distance of range.

  spname <- gsub(" ", "_", myspp)
  myfile <- list.files(wd$coords2, pattern = paste0(myspp, ".csv"), full.names = T)
  if(length(myfile) == 0) { myfile <- list.files(wd$coords2, pattern = paste0(gsub("_", " ", myspp), ".csv"), full.names = T)}
  stopifnot(length(myfile) != 0)
  cs <- fread(myfile)
  pts <- sf::st_as_sf(cs, coords = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs") %>%
    sf::st_transform(st_crs(world))

  if(!exists("iucnRangemaps")) iucnRangemaps <- read_sf("/Users/cjcampbell/BigZaddyData/IUCN_BatRanges/redlist_species_data_d9b37aa6-89b6-4fe4-9ebc-8567d9a63a8c/data_0.shp") %>%
    st_simplify(10) %>%
    sf::st_transform(st_crs(world))
  stopifnot(gsub("_", " ", myspp) %in% iucnRangemaps$BINOMIAL)
  range <- iucnRangemaps %>%
    dplyr::filter( BINOMIAL == gsub("_", " ", myspp)) %>%
    st_union()

  cs2 <- cs[ as.numeric(sf::st_distance(range, pts)) <= dist*1e3, ]

  return(as.numeric(cs2$id))
}

# Import occurrence records -----------------------------------------------
if(!exists("input_df")) {
  occs <- fread(file.path(wd$bin, "occs_syn.csv"))
  input_df <- occs %>%
    split(.$canonical)
  #rm(occs)
}


# Clean occurrence records (v1) -------------------------------------------

# order by descending nrow:
descOrder <- lapply(input_df, nrow) %>% unlist %>% data.frame(nrow = ., order = 1:length(input_df)) %>% arrange(desc(nrow)) %>% dplyr::select(order) %>% unlist()

#for(i in 1:length(input_df) ) {
for(i in descOrder ) {
  spp <- input_df[[i]][1, "canonical"]
  bw <- str_replace(spp, " ", "_")

  if(runForAll | spp %in% targetSpecies$canonical) {
    if(
      !file.exists(file.path(wd$coords, paste0(bw, ".csv"))) |
      !file.exists(file.path(wd$coords, paste0(bw, ".png")))
      ) {
        print(paste("Working on", spp))
        # Filter out coords close to 0,0
        df_exlude <- dplyr::filter(input_df[[i]], decimalLongitude > -1, decimalLongitude < 1, decimalLatitude > -1, decimalLatitude < 1 )
        df <- anti_join(input_df[[i]], df_exlude)
        if(nrow(df) == 0) { print("No points available..."); next }
        tryCatch(
          map_coord_clean( dataframe = df, spp ),
          error = function(e) print(paste(spp[i], "Error, did not run, skipping"))
        )
        Sys.sleep(10)
    } else {
      print(paste("File already exists for", spp))
    }
  }
}


# Manual check of records (v1) -------------------------------------------------
if(runForAll) {
  if(!exists("occs")) {
    occs <- fread(file.path(wd$bin, "occs_syn.csv"))
  }
  whichSpecies <- sort(unique(occs$canonical))
} else {
  whichSpecies <- sort(unique(targetSpecies$canonical))
}

rerunCheck <- FALSE
for(species in whichSpecies) {
  spname <- gsub(" ", "_", species)
  if(!file.exists(file.path(wd$coords, paste0(spname, "_out.rds"))) | isTRUE(rerunCheck) ) {
    print(species)

    # Does file exist
    filepath <-  list.files(
      wd$coords, pattern = paste0(spname, ".png"),
      recursive = T, full.names = T
    )
    if(length(filepath) == 0) { print("File doesn't exist; skipping"); next }

    # Open map
    filepath %>%
      browseURL()

    passfail <- readline(prompt= paste(
    "
    1 - Keep all points (pass)
    2 - Remove flagged only
    3 - remove flagged points except .otl (outliers)
    4 - remove / retain other points
    5 - Skip for now
    "))

    if(passfail != 5){
      data.frame(species = species, result = passfail) %>%
        saveRDS(file = file.path(wd$coords, paste0(spname, "_out.rds")))
    }
  }
}

checkOutcomes <- list.files(wd$coords, pattern = "_out.rds", full.names = T) %>%
  lapply(readRDS) %>%
  bind_rows()


# Pull out clean occurrences and prepare for revision (move to /bin/coords2) or repository (coords3) ---------------------

for(i in 1:nrow(checkOutcomes)) {

  mycoords <- fread(file.path(wd$coords, paste0(
    gsub(" ", "_", checkOutcomes[i, "species"]), ".csv"
  )))

  sppname <-  gsub("_", " ", paste0(checkOutcomes[i, "species"]))

  if( checkOutcomes[i, "result"] == 1) {
    mycoords %>%
      fwrite(., file.path(wd$coords2, paste0( sppname, ".csv")))
    mycoords %>%
      dplyr::select(all_of(targetCols)) %>%
      fwrite(., file.path(wd$coords3, paste0( sppname, ".csv")))
  } else if( checkOutcomes[i, "result"] == 2) {
    out <- mycoords %>%
      dplyr::filter(.summary == TRUE)
    fwrite(out, file.path(wd$coords2, paste0( sppname, ".csv")))
    out %>%
      dplyr::select(all_of(targetCols)) %>%
      fwrite(., file.path(wd$coords3, paste0( sppname, ".csv")))
  } else if( checkOutcomes[i, "result"] == 3) {
    out <- mycoords %>%
      dplyr::filter(.val == TRUE, .equ == TRUE,.cap == TRUE, .cen == TRUE, .inst == TRUE)
    fwrite(out, file.path(wd$coords2, paste0(sppname, ".csv")))
    out %>%
      dplyr::select(all_of(targetCols)) %>%
      fwrite(., file.path(wd$coords3, paste0( sppname, ".csv")))
    } else {
    mycoords %>%
      fwrite(.,  file.path(wd$coords2, paste0(sppname, ".csv" )))
  }
}


# Manually correct records ------------------------------------------------
## Reallocate pts from species from one to another -----
# This is typically due to taxon splits, some cases of frequent mis-IDs.

### Corynorhinus rafinesquii -> Corynorhinus townsendii ----
fileName <- list.files(wd$coords2, pattern = "Corynorhinus townsendii", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Corynorhinus rafinesquii", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLongitude < -100) %>%
  dplyr::mutate(canonical = "Corynorhinus townsendii")
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Hipposideros gentilis in southwest India -> Hipposideros pomona -----
fileName <- list.files(wd$coords2, pattern = "Hipposideros pomona.csv", full.names = T)
toMove <-  list.files(wd$coords2, pattern = "Hipposideros gentilis", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLongitude > 50 & decimalLongitude < 100) %>%
  dplyr::mutate(eventDate = as.POSIXct(eventDate)) %>%
  dplyr::mutate(canonical = "Hipposideros pomona")
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus cinereus in Hawaii -> Lasiurus semotus ----
fileName <- list.files(wd$coords2, pattern = "Lasiurus semotus", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus cinereus", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLongitude < -140) %>%
  dplyr::mutate(canonical = "Lasiurus semotus")
# quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus cinereus in south america -> Lasiurus villosissimus ----
fileName <- list.files(wd$coords2, pattern = "Lasiurus villosissimus", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus cinereus", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLatitude < 10 & decimalLongitude < -40) %>%
  dplyr::mutate(canonical = "Lasiurus villosissimus")
# quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus borealis in south america -> Lasiurus blossevillii ----
fileName <- list.files(wd$coords2, pattern = "Lasiurus blossevillii", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus borealis", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLatitude < 8 & decimalLongitude < -40) %>%
  dplyr::mutate(canonical = "Lasiurus blossevillii")
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus borealis in western and southern north america -> Lasiurus frantzii ----
states_us <- ne_states(geounit = "United States of America") %>%
  sf::st_as_sf()
uspts <-list.files(wd$coords2, pattern = "Lasiurus borealis", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLatitude > 8 & decimalLongitude < -40 ) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(states_us)) %>%
  sf::st_intersection(states_us) %>%
  dplyr::filter(name_pl %in% c("California", "Arizona", "New Mexico"))
# This might still have too many pts from NE mex...
fileName <- list.files(wd$coords2, pattern = "Lasiurus frantzii", full.names = T)
toMove1 <- input_df$`Lasiurus borealis` %>%
  dplyr::filter(decimalLatitude > 8 & decimalLongitude < -40 ) %>%
  dplyr::filter(
    decimalLatitude < 22  & decimalLongitude < -78 |
    decimalLatitude < 25  & decimalLongitude < -80 |
    decimalLatitude < 30  & decimalLongitude < -100 |
    decimalLatitude < 35  & decimalLongitude < -105
    ) %>%
  dplyr::mutate(canonical = "Lasiurus frantzii")
toMove <- input_df$`Lasiurus borealis` %>%
  dplyr::filter(id %in% c(uspts$id, toMove1$id))
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus ega to Lasiurus xanthinus ----
# Move points from western mexico and southwest us to L. xanthinus
fileName <- list.files(wd$coords2, pattern = "Lasiurus xanthinus", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus ega", full.names = T) %>%
  fread %>%
  dplyr::filter( decimalLongitude < -105 & (decimalLongitude < -100 | decimalLatitude > 20) ) %>%
  dplyr::mutate(canonical = "Lasiurus xanthinus")
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus_blossevillii in wesstern north america -> Lasiurus_frantzii ---------
fileName <- list.files(wd$coords2, pattern = "Lasiurus frantzii", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus blossevillii", full.names = T) %>%
  fread() %>%
  dplyr::filter( decimalLongitude < -80 & decimalLatitude > 10) %>%
  dplyr::mutate(canonical = "Lasiurus frantzii")
# quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Lasiurus_frantzii in south america -> Lasiurus_blossevillii ---------
fileName <- list.files(wd$coords2, pattern = "Lasiurus blossevillii", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Lasiurus frantzii", full.names = T) %>%
  fread %>%
  dplyr::filter( decimalLongitude > -80 & decimalLatitude < 10) %>%
  dplyr::mutate(canonical = "Lasiurus blossevillii")
#quickPlot(toMove)
fread(fileName) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

### Rhinolophus darlingi western southern africa -> Rhinolophus damarensis ----
fileName <- list.files(wd$coords2, pattern = "Rhinolophus damarensis.csv", full.names = T)
toMove <- list.files(wd$coords2, pattern = "Rhinolophus darlingi", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLongitude < 25 & decimalLatitude < -15) %>%
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  dplyr::mutate(canonical = "Rhinolophus damarensis")
#quickPlot(toMove)
fread(fileName) %>%
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

###  Hipposideros_pomona -> H. gentilis----
# move points from southeast asia and china to H. gentilis
fileName <- list.files(wd$coords2, pattern = "Hipposideros gentilis.csv", full.names = T)
toMove <-list.files(wd$coords2, pattern = "Hipposideros pomona", full.names = T) %>%
  fread %>%
  dplyr::filter(decimalLongitude > 90) %>%
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  dplyr::mutate(canonical = "Hipposideros gentilis")
# quickPlot(toMove)
fread(fileName) %>%
  dplyr::mutate(eventDate = as.Date(eventDate)) %>%
  rbind(toMove, fill = T) %>%
  fwrite(fileName)

## Correct the remove / retain other points codes (option 4) ----
manualSpecies <- dplyr::filter(checkOutcomes, !result %in% 1:3) %>% dplyr::select(species) %>% unlist


## Species-specific filters -----
### Afronycteris_nana ----
# remove true outlier poiints
fixPts(myspp = "Afronycteris_nana", expr = "decimalLongitude > -50 & .inst == TRUE & .cap == TRUE & (decimalLatitude > 3 | decimalLongitude > 3)")

### Corynorhinus_rafinesquii ----
# remove western points (townsendii)
fixPts("Corynorhinus_rafinesquii", expr = "decimalLongitude > -100", plot = F)

### Eptesicus_diminutus ----
# remove pt from el salvador
fixPts("Eptesicus_diminutus", expr = "decimalLatitude < 10", plot = FALSE)

### Eptesicus_gobiensis ----
# pts look fine actually.
fixPts("Eptesicus_gobiensis", expr = "decimalLatitude != 0", plot = FALSE)

### Eptesicus_nilssonii ----
# remove north american outlier but retain others.
fixPts("Eptesicus_nilssonii", expr = "decimalLongitude > -50", plot = FALSE)

### Eptesicus_serotinus ----
# Relying on rangemap from wikipedia (https://en.wikipedia.org/wiki/Serotine_bat#/media/File:Eptesicus_serotinusMap.png)
# Remove outliers in north america
# Remove Japan points
# Remove north african points
fixPts("Eptesicus_serotinus", expr = "decimalLongitude > -10 & decimalLongitude < 125 & (decimalLatitude > 35.96 | decimalLongitude > 0) & (decimalLatitude > 37.88 | decimalLongitude > 13 | decimalLongitude < 0)", plot = F)

### Eumops_bonariensi ----
## Redlist notes: "In its most restricted meaning, this species is found from northern Argentina, Paraguay, Uruguay, and Brazil (Eger 2008, Gamboa et al. 2016). This definition does not includes the taxa E. delticus or E. nanus, now treated as different species (Eger 2008)."
# remove flagged points
# remove everything north of paraguay
fixPts("Eumops_bonariensis", expr = ".cap == TRUE & .cen == TRUE & decimalLatitude < -20", plot = F)

### Eumops_nanu ----
# This is possibly low-quality considering the Colombian points, but I can't see much reason to filter them out-- they include recent points as well...
fixPts(myspp = "Eumops_nanus", expr = "decimalLatitude != 0", plot = F)

### Eumops_wilson ----
## Retain points in ecuador and peru
fixPts("Eumops_wilsoni", expr = "decimalLongitude < -74", plot = F)

### Hipposideros_alongensis ----
# This genus has had some recent research and splits within Vietnam...
# The existing pts are not particularly near the described locations in
# northeastern Vietnam, so I think best to remove them...
fixPts("Hipposideros_alongensis", expr = "decimalLongitude > 104", plot = F)

### Hipposideros_caffer ----
# This species has definitely been split along subspecies lines, see
# https://academic.oup.com/mspecies/article/doi/10.1644/845.1/2600887#114430231
# and look for who split what!
# See Vallo, Peter, et al. "Variation of mitochondrial DNA in the Hipposideros caffer complex (Chiroptera: Hipposideridae) and its taxonomic implications." Acta Chiropterologica 10.2 (2008): 193-206.
# According to Vallo et al, restrict to southeastern Africa
#
# removed outliers in asia and south america
# Removed points west of 20E
# Removed points north of 11S
# remove some points from northwestern Zambia
fixPts("Hipposideros_caffer", expr = "decimalLongitude > -50 & decimalLongitude < 75 & decimalLongitude > 20 & decimalLatitude < -11 & (decimalLatitude < -15 | decimalLongitude > 26)", plot = F)

### Hipposideros_gentilis ----
# Remove points way outside of range...
# getIDsOfInRangePts("Hipposideros_gentilis", dist = 500) %>% dput
fixPts("Hipposideros_gentilis", expr = "id %in% c(476828258, 476825710, 476864791, 476830947, 3892316235, 476830946,
476858892, 2640444661, 2640444767, 2640444784, 2636490301, 2636490303,
2636490434, 2636490438, 2636490442, 2636490444, 2636490448, 1919695068,
2640419303, 2640419305, 2640419316, 2640419365, 2640419438, 2640446304,
2640446312, 2640446354, 2640446483, 2640420336, 3892318555, 2636491301,
2636491315, 2636491326, 2636491333, 2636491337, 2636491362, 2636491377,
2636491391)", plot = F)

### Hipposideros_inornatu ----
# filter to points on same island as range map
fixPts("Hipposideros_inornatus", expr = "decimalLatitude < -11", plot = F)

### Hipposideros_kunzi ----
# See https://doi.org/10.3161/15081109ACC2018.20.1.001
# All are plausible, in penninsular malaysia...
fixPts("Hipposideros_kunzi", expr = "decimalLatitude != 0", plot = F)

### Hipposideros_pelingensi ----
# remove outlier in papau new guinea. retain one flagged outlier.
fixPts("Hipposideros_pelingensis", expr = "decimalLongitude < 140", plot = F)

### Hipposideros_pomona ----
# "This species belongs to bicolor species group. This species was listed under Hipposideros bicolor (Temminck, 1834), but is now considered distinct (Hill et al. 1986, Srinivasulu and Srinivasulu 2012). Based on morphology and bacular structure, Srinivasulu and Srinivasulu (2018) upgraded the taxon gentilis Andersen, 1918 to specific rank, and suggested that H. pomona Andersen, 1918 sensu stricto is restricted to southern India, and H. gentilis sensu lato occurs from northeast India to Southeast Asia."
# -- https://www.iucnredlist.org/species/180990825/180990948
# moved points from southeast asia and china to H. gentilis
fixPts("Hipposideros_pomona", expr = "decimalLongitude < 80", plot = F)

### Hipposideros_ruber ----
# remove points around 0,0 and south african outlier
# remove other flagged points
fread(file.path(wd$coords, paste0("Hipposideros_ruber", ".csv"))) %>%
  dplyr::mutate(id = as.numeric(id)) %>%
  dplyr::filter(`.cap` == TRUE & .equ == TRUE & (decimalLatitude > 1 | decimalLongitude > 1) & (decimalLongitude > -20)) %>%
  fwrite(file.path(wd$coords3, "Hipposideros ruber.csv") )

### Hipposideros_tephrus ----
# unusual very disjunct range...
# Remove pts outside of IUCN range, I guess.
# getIDsOfInRangePts("Hipposideros_tephrus") %>% dput
fixPts("Hipposideros_tephrus", expr = ".cap == TRUE & (id %in% c(2841673949, 3346452444, 1323016828))", plot = F)

### Lasiurus_borealis ----
# remove flagged points
# restrict to approximately North America
# Exclude south and central america. Exclude southern Mexico.
# Remove points from the west No Am
# Remove points from western mexico
# remove points from parts of bahamas
fixPts("Lasiurus_borealis", expr = ".cap == TRUE & .inst == TRUE &
       decimalLatitude > 19 &
       decimalLongitude < 0  &
       (decimalLatitude > 22 | decimalLongitude > -80) &
       decimalLongitude > -112 &
       (decimalLongitude > -105 | decimalLatitude > 35 | decimalLongitude > -95) &
       (decimalLongitude < -74 | decimalLatitude > 23)",
       plot = F)

### Lasiurus_cinereus ----
# remove outliers and points from south america, hawaii
fixPts("Lasiurus_cinereus", expr =
         "decimalLatitude > 12 &
       decimalLongitude < -50 &
       decimalLongitude > -150", plot = F)

### Lasiurus_ega ----
# Filter out flagged points
# Move western mexico and southwestern US points to Lasiurus xanthinus (western yellow bat)
fixPts("Lasiurus_ega", expr = ".cap == TRUE & .cen == TRUE & .inst == TRUE & decimalLongitude > -105 & (decimalLongitude > -100 | decimalLatitude < 2)", plot = F)

### Lasiurus_egregius ----
# my range map has since been revised; these points look good.
# "The species is known only from a few localities. It occurs in extreme southeastern Panama, but the extent of its range into northwestern Colombia is unknown; also occurs in French Guiana, northern and southern Brazil (Eisenberg 1989, Simmons 2005). Lowlands only (Reid 1997)."
# -- https://www.iucnredlist.org/species/11351/22119870
fixPts("Lasiurus_egregius", expr = "decimalLatitude != 0", plot = F)

### Lasiurus_frantzii ----
# Already moved south american records to blossevillii
# filter to north american records
# filter out other flagged points
fixPts("Lasiurus_frantzii", expr = ".cap == TRUE & .inst == TRUE & decimalLongitude < -80")

### Lasiurus_minor ----
fixPts("Lasiurus_minor", expr = "decimalLatitude < 24 & .cap == TRUE", plot = F)

### Molossus_alvarezi ----
# Batnames distribution "Yucatan Peninsula of Mexico, Belize, Guatemala, Honduras, south to Colombia, Venezuela, Surinam, French Guiana, Peru, Trinidad"
fixPts("Molossus_alvarezi", expr = "decimalLatitude > 10", plot = F)

### Molossus_aztecus ----
# Batnames distribution: "Jalisco (Mexico) to Nicaragua; Cozumel Isl (Mexico); S Venezuela, SE Brazil"
# real paucity of points for south america, and likely some spurious points in here... but none I can distinguish. Good enough?
fixPts("Molossus_aztecus", expr = "decimalLatitude != 0", plot = F)

### Molossus_bondae ----
# Distribution from batnames: "Honduras to Costa Rica; E Panama, Colombia, Ecuador, and Venezuela"
# Filter out Peru, points south...
fixPts("Molossus_bondae", expr = "decimalLatitude < 20 & decimalLatitude > -10", plot = F)

### Molossus_currentium ----
# Distribution from batnames: "Amazonian Brazil, Paraguay, Uruguay, and N Argentina"
# filtered to western points in viscinity of iucn range...
fixPts("Molossus_currentium", expr = "decimalLongitude > -65", plot = F)

### Mops_ansorgei ----
# Filter to African continent
fixPts("Mops_ansorgei", expr = "decimalLatitude < 20 & decimalLongitude > -20")

### Mops_bemmeleni ----
# Filter out point in cameroon; IUCN seems to have
# See https://www.iucnredlist.org/species/4307/22020379
fixPts("Mops_bemmeleni", expr = "decimalLongitude < 0 | decimalLongitude > 20", plot = F)

### Mops_bivittatus ----
# https://www.iucnredlist.org/species/4308/22020251
#Looks okay
fixPts("Mops_bivittatus", expr = "decimalLatitude != 0", plot = F)

### Mops_jobensis ----
# remove one outlier
fixPts("Mops_jobensis", expr = ".cen == TRUE & decimalLongitude > 0", plot = F)

### Mops_leucogaste ----
# remove outliers not in or near range
fixPts("Mops_leucogaster", expr = "decimalLatitude < -10", plot = F)

### Mops_sarasinoru ----
# remove african points...
fixPts("Mops_sarasinorum", expr = "decimalLongitude > 100", plot = F)

### Mormoops_megalophylla ----
# remove outlier in africa, other flagged points
# remove out-of-range pts in Caribbean
fixPts("Mormoops_megalophylla", expr = "decimalLongitude < -40 & .cap == TRUE & .inst == TRUE & (decimalLongitude < -85 | decimalLatitude < 15)", plot = F)

### Nyctalus_lasiopteru ----
# remove east asian points, other flagged points
fixPts("Nyctalus_lasiopterus", expr = "decimalLongitude < 100 & .cap == TRUE & .cen == TRUE", plot = F)

### Nyctalus_leisleri ----
# remove flagged points, north american outlier
fixPts("Nyctalus_leisleri", expr = "decimalLongitude > -50 & .cap == TRUE & .cen == TRUE & .inst == TRUE", plot = F)

### Nyctalus_noctula ----
# remove 2 south central asian points , other flagged points
fixPts("Nyctalus_noctula", expr = "decimalLatitude > 30 & .cap == TRUE & .cen == TRUE & .gbf == TRUE & .inst == TRUE", plot = F)

### Nycteris_thebaica ----
# remove outliers in southeast asia, americas
# remove other flagged points
# remove points near 0,0
fixPts("Nycteris_thebaica", expr = "decimalLongitude > -50 & decimalLongitude < 100 & .cap == TRUE & .cen == TRUE & .inst == TRUE & (decimalLatitude > 1 | decimalLongitude > 1)", plot = F)

### Plecotus_austriacus ----
# remove outliers from central asia
# remove outliers from middle east
# remove outliers from mid-sweden
fixPts("Plecotus_austriacus", expr = "decimalLongitude < 60 & (decimalLatitude > 38 | decimalLongitude < 31) & decimalLatitude < 57", plot = F)

### Plecotus_christie ----
# retain all points in north africa
fixPts("Plecotus_christiei", expr = "decimalLatitude != 0", plot = F)

### Plecotus_homochrous ----
# remove out-of-range pt
fixPts("Plecotus_homochrous", expr = "decimalLongitude < 100", plot = F)

### Plecotus_kolombatovic ----
# found in italy, retain all points
fixPts("Plecotus_kolombatovici", expr = "decimalLatitude != 0", plot = F)

### Pteronotus_davy ----
# remove outliers in north north america, africa
# remove other flagged points
# remove outlier from east brazil
fixPts("Pteronotus_davyi", expr = "decimalLatitude < 40 & decimalLongitude < -20 & .cap == TRUE & .cen == TRUE & decimalLongitude < -40", plot = F)

### Pteronotus_mesoamericanus ----
# remove distant outlier in venezuela
fixPts("Pteronotus_mesoamericanus", expr = "decimalLongitude < -70", plot = F)

### Pteronotus_mexicanus ----
# BatNames: "Sonora and Tamaulipas to Oaxaca and Veracruz (Mexico)"
# Looks good to me!
fixPts("Pteronotus_mexicanus", expr = "decimalLatitude != 0", plot = F)

### Rhinolophus_clivosus ----
# remove points outside of africa and middle east
fixPts("Rhinolophus_clivosus", expr = "(decimalLatitude > 1 | decimalLongitude > 1) & decimalLongitude < 60", plot = F)

### Rhinolophus_convexus ----
# Looks fine, one point in a disjunct range area...
fixPts("Rhinolophus_convexus", expr = "decimalLatitude != 0", plot = F)

### Rhinolophus_damarensis ----
# Referencing https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082614
# R. darlingi is found in east southern africa.
# $. damarensis is found in west southern africa
# After transferring pts from darlingi -> damarensis, looks good.
fixPts("Rhinolophus_damarensis", expr = "decimalLatitude < -15", plot = F)

### Rhinolophus_darlingi ----
# Referencing https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082614
# R. darlingi is found in east southern africa.
# $. damarensis is found in west southern africa
# Remove outlier points that are damarensis or not in southern africa...
fixPts("Rhinolophus_darlingi", expr = "decimalLatitude < -15 & decimalLongitude > 23", plot = F)

### Rhinolophus_denti ----
# this range gives me anxiety...
# Remove stuff way out of range of the WILDLY DISCONNECTED ranges from the range map...
fixPts("Rhinolophus_denti", expr = "decimalLongitude < 0 | (decimalLatitude < -15 & decimalLongitude < 28)", plot = F)

### Rhinolophus_euryale ----
# remove way out-of-range outlier
# remove other flagged points
fixPts("Rhinolophus_euryale", expr = "decimalLatitude > 10 & .cap == TRUE & .cen == TRUE", plot = F)

### Rhinolophus_ferrumequinum ----
# remove some outlier points
fixPts("Rhinolophus_ferrumequinum", expr = "decimalLatitude > 2 & decimalLongitude > -50", plot = F)

### Rhinolophus_hipposidero ----
# remove outlier points in south america and central africa
fixPts("Rhinolophus_hipposideros", expr = "decimalLatitude > 10 & .cap == TRUE & .cen == TRUE & .inst == TRUE", plot = F)

### Rhinolophus_landeri ----
# remove 0,0 and flagged points
# central namibia seems likely misidentification (https://biodiversity.org.na/taxondisplay.php?nr=345)
fixPts("Rhinolophus_landeri", expr = ".cap == TRUE & (decimalLatitude > 2 | decimalLongitude > 2) & decimalLatitude > -19", plot = F)

### Rhinolophus_lobatus ----
# Per https://academic.oup.com/zoolinnean/article/184/4/1249/4984486,
# this species is only found in Mozambique.
# Remove all pts, I guess?
fixPts("Rhinolophus_lobatus", expr = "decimalLongitude > 40", plot = F)

### Rhinolophus_maclaud ----
# Retain only the recent point kind of near the range, I guess?
# Used a manual check here...
fixPts("Rhinolophus_maclaudi", expr = "decimalLongitude < 0", plot = F)

### Rhinolophus_megaphyllus ----
# remove some outliers
# remove additional outliers -- range restricted to east australia and new guinea...
fixPts("Rhinolophus_megaphyllus", expr = ".cap == TRUE & .cen == TRUE & .inst == TRUE & decimalLongitude > 0 & decimalLatitude > -50 & decimalLongitude > 130", plot = F)

### Rhinolophus_mehelyi ----
# remove outliers in south america, near 0,0, in central Algeria
# remove out-of-range points in northern france, czechia
# seems like most gbif points from northern spain are since 2000 and from multiple sources (likely legit)
# remove points from france (almost exclusively from before 1950) that do not seem legit...
fixPts("Rhinolophus_mehelyi", expr = "decimalLatitude > 28 & decimalLatitude < 45 & (decimalLatitude < 43 | decimalLongitude > 10 | decimalLongitude < -2)", plot = F)

### Rhinolophus_pearsonii ----
# remove point in sri lanka
fixPts("Rhinolophus_pearsonii", expr = "decimalLatitude > 10 | decimalLongitude > 82", plot = F)

### Rhinolophus_rex ----
# Notes from IUCN:
# https://www.iucnredlist.org/species/19562/21994639
# "This species is possibly conspecific with Rhinolophus paradoxolophus (Corbet and Hill 1992). Zhang et al. (2009) also referred that R. paradoxolophus may represents a subspecies of R. rex. Until recently, phylogenies, genetic and phenotypic divergence, and species delimitation analyses supported the revised status of R. paradoxolophus as a subspecies of R. rex (Zhang et al. 2018)."
# As a chinese endemic, remove pts well outside of range...
# getIDsOfInRangePts("Rhinolophus_rex", dist = 100) %>% dput
fixPts("Rhinolophus_rex", expr = "id %in%
       c(3346733971, 3346722820, 3350146509, 2636490304, 2636490305,
        2636490433, 1919791389, 665811774, 665803221, 2640446309, 3349901979,
        3349980054, 3346883405, 3349973051)", plot = F)

### Rhinolophus_rouxi ----
## Batnames: "Sri Lanka, peninsular India to S Myanmar; possibly S China. Reports of this species from Cambodia are likely erroneous; see Kock (2000a)"
fixPts("Rhinolophus_rouxii", expr = "decimalLatitude != 0", plot = F)

### Rhinolophus_tatar ----
# To filter out points from an island neighboring Sulawese?
# The points are from 2018, unlikely mis-IDed. Keep them all!
fixPts("Rhinolophus_tatar", expr = "decimalLatitude != 0", plot = F)

### Scotophilus_dingani ----
# remove one outlier and other flagged pts.
fixPts("Scotophilus_dinganii", expr = ".cap == TRUE & decimalLatitude > -40", plot = F)

### Scotophilus_kuhli ----
# remove one outlier and other flagged pts.
fixPts("Scotophilus_kuhlii", expr = ".cap == TRUE & .cen == TRUE & decimalLongitude > 0", plot = F)

### Scotophilus_nigrita ----
# Comment: Weird disjunct range...
# Will remove points well outside of it
#
# remove cap flagged pts
# remove points far outside range
# getIDsOfInRangePts("Scotophilus_nigrita", dist = 300) %>% dput()
fixPts("Scotophilus_nigrita", expr = ".cap == TRUE & id %in% c(1065379358, 1065366168, 1044009693, 1145705789, 3025958735, 1211806532, 665898991, 665890707, 3025989896, 1322098000, 1932239984, 1932239874, 1932239932, 1932239938, 1932239934, 919418089, 3024394987, 1305201465, 686434082, 1317829889)", plot = F)

### Scotophilus_nigritellus ----
# remove extra pt near 0,0
fixPts("Scotophilus_nigritellus", expr = ".summary == TRUE & decimalLatitude > 1", plot = F)

### Scotophilus_nucella ----
# remove outlier
fixPts("Scotophilus_nucella", expr = "decimalLongitude < 20", plot = F)

### Tadarida_aegyptiaca ----
# remove flagged points excepting outliers
# remove outlier in australia
fixPts("Tadarida_aegyptiaca", expr = ".cap == TRUE & .inst == TRUE & decimalLongitude < 100", plot = F)

### Tadarida_brasiliensis ----
# remove some outliers, retain all south america.
fixPts("Tadarida_brasiliensis", expr = ".cap == TRUE & .cen == TRUE & .inst == TRUE & decimalLongitude < 0 & decimalLongitude > -140", plot = F)

### Taphozous_mauritianus ----
# remove outlier from north america
# remove non-outlier flagged points
# remove points hovering around 0,0
# remove outlier from north africa
fixPts("Taphozous_mauritianus", expr = "decimalLongitude > -50 & .cap == TRUE & .cen == TRUE & .inst == TRUE & (decimalLatitude > 2 | decimalLongitude > 2) & decimalLatitude < 20", plot = F)


# Check all transferred-to and transferred-from species again -----------------------------------------

### Corynorhinus townsendii ----
# Remove some small points out of range
fixPts("Corynorhinus townsendii", expr =
         ".cap == TRUE & .inst == TRUE & .cen == T &
       (decimalLatitude > 32 | decimalLatitude < 29 |decimalLongitude < -100) &
       (decimalLongitude < -90 | decimalLongitude > -88)", plot = F)

# fread(file.path(wd$coords2, "Corynorhinus townsendii.csv")) %>%
#   dplyr::mutate(id = as.numeric(id)) %>%
#   dplyr::filter(
#       (decimalLatitude > 32 | decimalLatitude < 29 |decimalLongitude < -100) &
#       (decimalLongitude < -90 | decimalLongitude > -88)
#     ) %>% clean_coordinates(
#       lon = "decimalLongitude", lat = "decimalLatitude",
#       species = "canonical",
#       tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers"),
#       verbose = T) %>%
#   dplyr::filter(.summary == TRUE) %>%
#   fwrite(file.path(wd$coords3, "Corynorhinus townsendii.csv") )


### Hipposideros pomona -----
fixPts("Hipposideros pomona", expr = ".cap == TRUE & .gbf == TRUE & decimalLongitude < 80", plot = F)

### Lasiurus semotus ----
fixPts("Lasiurus semotus", expr = "decimalLongitude != 0", plot = F)

### Lasiurus villosissimus ----
fixPts("Lasiurus villosissimus", expr = ".cap == TRUE & .inst == TRUE", plot = F)

### Lasiurus blossevillii ----
fixPts("Lasiurus_blossevillii", expr = ".cap == TRUE &
       decimalLongitude > -90 &
       decimalLatitude < 15 &
       (decimalLatitude < 10 | decimalLongitude > -81)",
       plot = F)

# fread(file.path(wd$coords2, "Lasiurus blossevillii.csv")) %>%
#   dplyr::mutate(id = as.numeric(id)) %>%
#   dplyr::filter(
#       decimalLongitude > -90 &
#       decimalLatitude < 15 &
#       (decimalLatitude < 10 | decimalLongitude > -81)
#   ) %>% clean_coordinates(
#   lon = "decimalLongitude", lat = "decimalLatitude",
#   species = "canonical",
#   tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers"),
#   verbose = T) %>%
#   dplyr::filter(.summary == TRUE) %>%
#   fwrite(file.path(wd$coords3, "Lasiurus blossevillii.csv") )

### Lasiurus frantzii ----
fixPts("Lasiurus frantzii", expr = ".cap == TRUE & .inst == TRUE & (decimalLongitude < -80 & decimalLatitude > 10)", plot = F)

### Lasiurus xanthinus ----
fixPts("Lasiurus xanthinus", expr = "decimalLongitude != 0", plot = F)

### Rhinolophus damarensis ----
fixPts("Rhinolophus damarensis", expr = "decimalLongitude != 0", plot = F)

###  H. gentilis ----
fixPts("Hipposideros gentilis", expr = "decimalLongitude != 0", plot = F)


# # Move all remaining approved species to coords3 --------------------------
#
# for(myspp in checkOutcomes$species) {
#   if(file.exists(file.path(wd$coords3, paste0(myspp, ".csv")))) {
#     next()
#   } else if(
#     !(file.exists(file.path(wd$coords3, paste0(myspp, ".csv")))) &
#     file.exists(file.path(wd$coords2, paste0(myspp, ".csv"))) &
#     (checkOutcomes[checkOutcomes$species == myspp, "result"] %in% 1:3)
#   ) {
#     file.copy(from=file.path(wd$coords2, paste0(myspp, ".csv")), to= file.path(wd$coords3, paste0(myspp, ".csv")),
#               overwrite = F, copy.mode = TRUE, copy.date = T)
#   } else { next() }
# }


# Synthesize and tidy dates -----------------------------------------------

dat_clean <- list.files(wd$coords3, full.names = T) %>%
  lapply(function(x) {
    fread(x) %>%
      dplyr::select(canonical, any_of(targetCols)) %>%
      dplyr::mutate(eventDate = lubridate::date(eventDate))
    }) %>%
  rbindlist(fill = T)

fwrite(dat_clean, file = file.path(wd$data, "occurrences_clean.csv"), row.names = F)

dat_clean %>%
  dplyr::filter(canonical %in% targetSpecies$canonical) %>%
  count(canonical)


# out <- list.files(wd$coords3, full.names = T) %>%
#   lapply(function(x) {
#     howMany <- fread(x) %>%
#       dplyr::filter( canonical == "Lasiurus blossevillii") %>%
#       nrow()
#     if(howMany != 0) {
#       writeLines(x)
#       return(x)
#     }
#     })
