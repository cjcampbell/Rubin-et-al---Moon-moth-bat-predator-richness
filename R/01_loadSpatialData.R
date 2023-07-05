library(rnaturalearth)
library(sf)

# Create full map.
map_wrld <- rnaturalearthdata::countries110 %>% 
  st_as_sf()
saveRDS(map_wrld, file = file.path(wd$bin, "map_wrld.rds"))

# Also save a lower-res version
map_wrld_smpl <- st_cast(map_wrld,  "POLYGON", do_split = FALSE) %>% 
  sf::st_simplify(dTolerance = 10000)
saveRDS(map_wrld_smpl, file = file.path(wd$bin, "map_wrld_smpl.rds"))
