
# Load batnames -----------------------------------------------------------

batNames <- file.path(wd$data, "batNames", "Chiroptera2022-09-07.csv") %>%
  fread() %>%
  dplyr::mutate(canonical = stringr::str_squish(paste(Genus, Species)))


# Determine target species -----------------------------------------------------

targetGenera <- file.path(
  wd$data, "Bat_genera_for_Caitlin_tidy.csv"
) %>% fread()

# Load post-hoc sheet with removals.
toKeep <- file.path(wd$data, "canonical_include_list.csv") %>% fread

# How many species are we targeting?
targetSpecies <- batNames %>%
  dplyr::left_join(toKeep) %>% 
  dplyr::filter( Include == 1 )
nrow(targetSpecies)
# 433 as of most recent check.
# Down to 390 after some species removal (20230318)


