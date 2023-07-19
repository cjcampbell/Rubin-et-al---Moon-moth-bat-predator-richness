source("R/00_setup.R")

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
