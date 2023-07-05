source("R/00_setup.R")
# Return list of bat species included in model, with specified overlaps with some key countries.

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_transform() %>% 
  st_make_valid() %>% 
  st_simplify(dTolerance = 10e3) 

easternNoAm <- world %>% 
  dplyr::filter(continent == "North America") %>% 
  st_crop(st_bbox(c(xmin = -110, xmax = -50, ymax = 55, ymin = 0), crs = st_crs(4326)) ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 1e3) %>% 
  st_union() 

europe <- world %>% 
  dplyr::filter(name %in% c("Spain", "France", "Switzerland")) %>% 
  st_crop(st_bbox(c(xmin = -10, xmax = 40, ymax = 90, ymin = 0), crs = st_crs(4326)) ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 400e3) %>% 
  st_union() 

africa <- world %>% 
  dplyr::filter(name %in% c("Kenya", "Tanzania", "Malawi", "Mozambique", "Zimbabwe", "South Africa", "Madagascar", "Namibia")) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 300e3) %>% 
  st_union() 

india <- world %>% 
  dplyr::filter(name %in% c("India", "Bangledesh", "Nepal", "Bhutan", "Sri Lanka")) %>% 
  st_difference(st_bbox(c(xmin = 85, xmax = 100, ymax = 15, ymin = 0), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 200e3) %>% 
  st_union() 

palearctic_asia <- world %>% 
  dplyr::filter(name %in% c("China", "Korea", "Dem. Rep. Korea", "Japan", "Taiwan")) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 200e3) %>% 
  st_union() 

malaysianPenninsula <- world %>% 
  dplyr::filter(name %in% c("Malaysia", "Thailand", "Vietnam", "Cambodia", "Lao PDR", "Myanmar")) %>% 
  st_crop(st_bbox(c(xmin = 0, xmax = 110, ymax = 90, ymin = 0), crs = st_crs(4326)) ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 100e3) %>% 
  st_difference( st_bbox(c(xmin = 106, xmax = 180, ymax = 5, ymin = 0), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_union() 
ggplot() +
  geom_sf(world, mapping = aes(), fill = "white") +
  geom_sf(malaysianPenninsula, mapping = aes(), fill = "red", alpha = 0.5) +
  coord_sf(xlim = sf::st_bbox(malaysianPenninsula)[c(1,3)],
           ylim = sf::st_bbox(malaysianPenninsula)[c(2,4)])

greaterSunda <- world %>% 
  dplyr::filter(name %in% c("Malaysia","Brunei","Indonesia")) %>% 
  st_crop(st_bbox(c(xmin = 0, xmax = 130, ymax = 90, ymin = -90), crs = st_crs(4326)) ) %>% 
  st_difference( st_bbox(c(xmin = 118.8, xmax = 180, ymax = 15, ymin = -8), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_difference( st_bbox(c(xmin = 115, xmax = 180, ymax = -5, ymin = -20), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_difference( st_bbox(c(xmin = 99, xmax = 105, ymax = 90, ymin = 3), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_difference( st_bbox(c(xmin = 101, xmax = 105, ymax = 90, ymin = 2), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_difference( st_bbox(c(xmin = 102.5, xmax = 105, ymax = 90, ymin = 1), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 50e3) %>% 
  st_union() 

lesserSunda <- world %>% 
  dplyr::filter(name %in% c("Indonesia")) %>% 
  st_crop(st_bbox(c(xmin = 115, xmax = 130, ymax = 90, ymin = -90), crs = st_crs(4326)) ) %>% 
  st_difference( st_bbox(c(xmin = 0, xmax = 118, ymax = 90, ymin = -5), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_difference( st_bbox(c(xmin = 0, xmax = 120, ymax = 90, ymin = 0), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 50e3) %>% 
  st_union() 

philippines <- world %>% 
  dplyr::filter(name %in% c("Philippines")) %>% 
  #st_crop(st_bbox(c(xmin = 115, xmax = 130, ymax = 90, ymin = -90), crs = st_crs(4326)) ) %>% 
  #st_difference( st_bbox(c(xmin = 0, xmax = 118, ymax = 90, ymin = -5), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  #st_difference( st_bbox(c(xmin = 0, xmax = 120, ymax = 90, ymin = 0), crs = st_crs(4326)) %>% sf::st_as_sfc() %>% st_as_sf() ) %>% 
  st_make_valid() %>% 
  st_buffer(dist = 50e3) %>% 
  st_union() 

# Plot --------------------------------------------------------------------

p_plot <- ggplot() +
  geom_sf(world, mapping = aes(), fill = "white") +
  geom_sf(easternNoAm, mapping = aes(), fill = "blue", alpha = 0.5, color = NA) +
  geom_sf(europe, mapping = aes(), fill = "orange", alpha = 0.5, color = NA) +
  geom_sf(africa, mapping = aes(), fill = "purple", alpha = 0.5, color = NA) +
  geom_sf(india, mapping = aes(), fill = "darkgoldenrod", alpha = 0.5, color = NA) +
  geom_sf(palearctic_asia, mapping = aes(), fill = "darkgreen", alpha = 0.5, color = NA) +
  geom_sf(malaysianPenninsula, mapping = aes(), fill = "grey50", alpha = 0.5, color = NA) +
  geom_sf(greaterSunda, mapping = aes(), fill = "red", alpha = 0.5, color = NA) +
  geom_sf(lesserSunda, mapping = aes(), fill = "lightblue", alpha = 0.5, color = NA) +
  geom_sf(philippines, mapping = aes(), fill = "lightgreen", alpha = 0.5, color = NA) +
  coord_sf(crs = proj4_eqe) +
  theme_void()
ggsave(p_plot, file = file.path(wd$figs, paste0("batRegions-", Sys.Date(),".png")), width = 10, height = 10)


# Extract for species -----------------------------------------------------

source("R/00_setup.R")

world_eqe <- world %>% 
  st_transform(proj4_eqe)

# Combine regions.
regionNames <- c("easternNoAm", "europe", "africa", "india", "palearctic_asia", "malaysianPenninsula", "greaterSunda", "lesserSunda", "philippines")
myRegions <- regionNames %>% 
  lapply(function(x) {cbind(data.frame(region = x), get(x))}) %>% 
  bind_rows() %>% 
  st_as_sf() %>% 
  st_transform(proj4_eqe)

# Find continuous surface files.
whichFiles <- list.files(wd$out, pattern = "SDM.tif$", full.names = T)

# Filter to target species.
source('R/02b_targetSpecies.R')
whichFiles <- basename(whichFiles) %>% 
  gsub("_SDM.tif", "", .) %>% 
  gsub("\\.", " ", .) %>% 
  {whichFiles[. %in% targetSpecies$canonical]}
(sppNum <- length(whichFiles))

if(!file.exists(file.path(wd$out, "batRegions.csv"))) {
  
  mylist <- list()
  for(x in whichFiles) {
    print(paste(basename(x), " - file", which(whichFiles == x), "of", length(whichFiles)))
    # Convert raster to polygon
    mypoly <- rast(x) %>% 
      as.polygons() %>% 
      st_as_sf() %>% 
      st_union() 
    sf::st_crs(mypoly) <- proj4_eqe
    o <- st_intersects(myRegions, mypoly, sparse = F)
    out <- data.frame(species = basename(x), region = regionNames, intersects = as.numeric(o))
    mylist[[length(mylist) + 1]] <- out
  }
  
  # fwrite(mylist, file = file.path(wd$bin, "batRegions_notClean.csv"))
  
  mylist %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(
      species = gsub("_SDM tif","", gsub("\\.", " ", species))
    ) %>% 
    tidyr::pivot_wider(names_from = region, values_from = intersects) %>% 
    fwrite(file = file.path(wd$out, "batRegions.csv"))
  
}

# Plot richness in focal regions ------------------------------------------
# Load in richness surface.
mySurface <- rast(file.path(wd$out, "combinedContinuousSurface.tif"))
set.crs(mySurface, proj4_eqe)
focalRegions <- c( easternNoAm, europe, africa, india , palearctic_asia, malaysianPenninsula, greaterSunda, lesserSunda, philippines ) %>% 
  st_transform(proj4_eqe) %>% 
  terra::vect()
croppedSurface <- mask(mySurface, focalRegions)
croppedSurface_df <- croppedSurface %>% 
  as.data.frame(xy=T, na.rm = T) %>% 
  rename("value" = 3)
# Load worldmap.
wrld <- rnaturalearth::countries110 %>% 
  st_as_sf() %>% 
  st_transform(proj4_eqe)
# Load moth occurrences.
mothOccurrences <- fread(file.path(wd$data, "mothOccurrences.csv")) %>% 
  dplyr::mutate(id = row_number())
mothOccurrences_sf <- mothOccurrences %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = proj4_wgs) %>% 
  st_transform(proj4_eqe) %>% 
  full_join(mothOccurrences)

# Plot
p_map <- ggplot() +
  geom_tile(croppedSurface_df, mapping = aes(x=x,y=y,fill=value,color=value)) +
  geom_sf(wrld, mapping = aes(), fill = NA, color = "grey90", size = 0.25) +
  geom_sf(mothOccurrences_sf, mapping = aes(), shape = 4, color = "white") +
  coord_sf() +
  scale_fill_viridis_c( option = "turbo") +
  scale_color_viridis_c(option = "turbo") +
  theme_void() +
  theme(
    legend.text = element_text(color = "grey90"),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "black")
  )
ggsave(p_map, filename = file.path(wd$figs, paste0("richnessMap_clogLog-", Sys.Date(),".png")), width = 15, height = 10)
