
library(data.table)
library(taxotools)

source("~/GeoBatDiv/R/02b_targetSpecies.R")

# List known fossil species -----------------------------------------------

fossilOnlySpecies <- c(
  "Archaeonycteris relicta",
  "Cardioderma leakeyi",
  "Desmodus draculae",
  "Desmodus stocki",
  "Diphylla centralis",
  "Icaronycteris index",
  "Miomyotis floridanus",
  "Myzopoda africana",
  "Nycticeinops serengetiensis",
  "Palaeochiropteryx tupaiodon",
  "Phyllonycteris major",
  "Pseudorhinolophus schlosseri",
  "Pteronotus pristinus",
  "Scotoecus olduvensis"
  )

# Assemble ITIS synonyms list --------------------------------------------------
rerunSynonyms <- FALSE
if(rerunSynonyms) {

  # Apply for all bat species of the world.
  # Batch it up to avoid errors.
  x <- c(seq(1,nrow(batNames), by = 100), nrow(batNames))

  tmpdir <- file.path(tempdir(), "batsynonyms")
  if(!dir.exists(tmpdir)) dir.create(tmpdir)

  for(i in 2:length(x)) {
    target <- x[i-1]:x[i]
    print(paste0("Working on ", min(target), "-",max(target) ))
    targetNames <- batNames$canonical[target]
    batsynonyms <- list_itis_syn(targetNames)
    saveRDS(batsynonyms, file = file.path(tmpdir, paste0(min(target), "_",max(target), ".rds")))
  }

  batsynonyms <- list.files(tmpdir, pattern = ".rds", full.names = T) %>%
    lapply(readRDS) %>%
    bind_rows

  fwrite(batsynonyms, file = file.path(wd$bin, "bat_synonyms.csv"))
}

# Apply ITIS synonyms ----------------------------------------------------------

# Load synonyms
batSynonyms <- fread( file.path(wd$bin, "bat_synonyms.csv") ) %>%
  # 'canonical' is the name used by batnames
  # 'synonym' is the name used by GBIF.
  # Rename for clarity and easy joining.
  dplyr::rename(name_batnames = canonical, name_occs = synonym)

# Load occurrence data.
occs <- fread( file.path(wd$bin, "occurrences_tidy.csv") ) %>%
  dplyr::filter(
    # Remove occs w/o occurrennces.
    !is.na(decimalLongitude),
    !is.na(decimalLatitude),
    # Remove unknown species.
    species != "",
    # Remove fossil specimens and species.
    basisOfRecord != "FOSSIL_SPECIMEN",
    !species %in% fossilOnlySpecies
    )

# Species to look up
toLookUp <- unique(occs$species)
length(toLookUp)

# How many need a synonym?
needsSynonym <- toLookUp[!toLookUp %in% batNames$canonical]
length(needsSynonym)

# How many of those have synonyms?
hasSynonym <- needsSynonym[needsSynonym %in% batSynonyms$name_occs]
needsMore  <- needsSynonym[!needsSynonym %in% batSynonyms$name_occs]
length(needsMore)
# How many data are affected?
sum(occs$species %in% needsMore)

# List the most common examples.
dplyr::filter(occs, species %in% needsMore) %>% dplyr::group_by(species) %>% dplyr::summarise(n=n()) %>% arrange(desc(n))


# Get wiki synonyms -------------------------------------------------------

# Search for remaining species to match. Retain those with synonyms in batNames.
if(rerunSynonyms) {

  mylist <- lapply(needsMore, function(query){
    try({response <- list_wiki_syn(query)})
    if(exists("response")) return(response)
  })

  wiki <- mylist %>%
    purrr::discard(is.null) %>%
    dplyr::bind_rows()

  wikiSyns <- full_join(
    dplyr::select(batNames, canonical),
    dplyr::select(wiki, WikiName, Name),
    by = c("canonical" = "WikiName")
    ) %>%
    dplyr::filter(!is.na(Name), canonical != Name) %>%
    distinct() %>%
    dplyr::rename(name_batnames = canonical, name_occs = Name)
  fwrite(wikiSyns, file = file.path(wd$bin, "wikiSyns.csv"))
} else {
  if(!exists("wikiSyns")) wikiSyns <- fread( file.path(wd$bin, "wikiSyns.csv") )
}


# Apply manual corrections ------------------------------------------------
# What's left?
needsMore2 <- needsMore[!needsMore %in% wikiSyns$name_occs]

manualSyns <- data.frame(
  # Chaerephon -> Mops
  name_occs     = needsMore2[grep("Chaerephon", needsMore2)],
  name_batnames = gsub("Chaerephon", "Mops", needsMore2[grep("Chaerephon", needsMore2)])
  )
manualSyns <-  rbind(manualSyns, data.frame(
  # Histiotus -> Eptesicus
  name_occs     = needsMore2[grep("Histiotus", needsMore2)],
  name_batnames = gsub("Histiotus", "Eptesicus", needsMore2[grep("Histiotus", needsMore2)])
))

manualSyns <- rbind(manualSyns, data.frame(
  name_occs = c(
    "Anoura aequatoris",            # synonym of Anoura caudifer
    "Anoura peruana",               # synonym of Anoura geoffroyi
    "Rhinolophus achilles",         # subspecies of Rhinolophus philippinensis
    "Rhinolophus alticolus",        # subspecies of Rhinolophus simulator
    "Rhinolophus axillaris",        # synonym of Rhinolophus landeri
    "Rhinolophus gorongosae",       # synonym of Rhinolophus simulator
    "Rhinolophus ferrum-eqtjinum",  # Typo
    "Rhinolophus nanus",            # synonym of Rhinolophus keyensis
    "Myotis keenii",                # subspecies of Myotis evotis
    "Myotis peninsularis",          # subspecies of Myotis velifer
    "Myotis browni",                # subspecies of Myotis muricola
    "Trachops fuliginosus",         # synonym of Trachops cirrhosus
    "Tonatia saurophila",           # subspecies of Tonatia bakeri
    "Lasiurus blossebillii",        # Typo
    "Thainycteris torquatus",       # Presence of multiple spellings?
    "Anthops omatus",               # Another spelling issue?
    "Roussettus aegyptiacus",       # Spelling -- too many s's
    "Rhinolophus robertsi",         # subspecies of Rhinolophus philippinensis
    "Rhinolophus perniger",         # subspecies of Rhinolophus luctus
    "Rhinolophus monticolus",       # synonym of Rhinolophus chutamasae
    "Rhinolophus kahuzi",           # synonym of Rhinolophus ruwenzorii
    "Rhinolophus bembanicus",       # synonym of Rhinolophus simulator
    "Rhinolophus brockmani",        # brockmani - synonym of Rhinolophus blasii
    "Rhinolophus dobsoni",          # dobsoni - synonym of Rhinolophus landeri
    "Rhinolophus indorouxii",       # indorouxii - synonym of Rhinolophus rouxii
    "Balantiopteryx infosca",       # Spelling
    "Chironax tumulus" ,            # subspecies of Chironax melanocephalus
    "Kerivoula malpasi",            # synonym of Kerivoula hardwickii
    "Kerivoula crypta",             # synonym of Kerivoula hardwickii
    "Afronycteris nanus",           # apparent spelling issue, nanus v. nana (same name, date on type specimen)
    "Miniopterus blepotis",         # subspecies of Miniopterus fuliginosus
    "Harpyiocephalus esperoptenus", # Subspecies of Hesperoptenus doriae
    "Hipposideros pratti" ,         # synonym of Hipposideros swinhoei
    "Miniopterus arenarius",        # subspecies of Miniopterus natalensis
    "Miniopterus eschscholtzii",    # subspecies of Miniopterus fuliginosus
    "Rhinolophus geoffroyii",       # Weird spelling prevented linking w/ synonym of R. clivosus
    "Rhinolophus etjryale",         # wiki synonym of Rhinolophus euryale
    "Rhogeessa gracilis",           # synonym of Baeodon gracilis
    "†Pipistrellus murrayi"         # Weird extra character

    ),
  name_batnames = c(
    "Anoura caudifer",
    "Anoura peruana",
    "Rhinolophus philippinensis",
    "Rhinolophus simulator",
    "Rhinolophus landeri",
    "Rhinolophus simulator",
    "Rhinolophus ferrumequinum",
    "Rhinolophus keyensis",
    "Myotis evotis",
    "Myotis velifer",
    "Myotis muricola",
    "Trachops cirrhosus",
    "Tonatia bakeri",
    "Lasiurus blossevillii",
    "Thainycteris torquatus",
    "Anthops ornatus",
    "Roussettus aegyptiacus" ,
    "Rhinolophus philippinensis",
    "Rhinolophus luctus",
    "Rhinolophus chutamasae",
    "Rhinolophus ruwenzorii",
    "Rhinolophus simulator",
    "Rhinolophus blasii",
    "Rhinolophus landeri",
    "Rhinolophus rouxii",
    "Balantiopteryx infusca",
    "Chironax melanocephalus",
    "Kerivoula hardwickii",
    "Kerivoula hardwickii",
    "Afronycteris nana",
    "Miniopterus fuliginosus",
    "Hesperoptenus doriae",
    "Hipposideros swinhoei",
    "Miniopterus natalensis",
    "Miniopterus fuliginosus",
    "Rhinolophus clivosus",
    "Rhinolophus euryale",
    "Baeodon gracilis",
    "Pipistrellus murrayi"
  )
))

(needsMore3 <- needsMore2[!needsMore2 %in% manualSyns$name_occs] %>% sort() )

# Notes: no synonyms available for:
# "Nycteris madagascariensis"
# "Pipistrellus montanus"
# "Xeronycteris vierai" 5 records



# Assemble w/ synonyms ----------------------------------------------------

allSynonyms <- bind_rows(batSynonyms, wikiSyns, manualSyns) %>% distinct()

newTargetCols <- c(
  "id", "institutionCode", "collectionCode", "basisOfRecord", "occurrenceID",
  "order", "family", "genus", "species", "specificEpithet", "scientificName",
  "identifiedBy", "eventDate", "year", "month", "day",  "verbatimEventDate",
  "recordedBy", "lifeStage", "decimalLongitude", "decimalLatitude", "locality",
  "county", "coordinateUncertaintyInMeters", "source",
  "canonical"
)

occs_syn <- occs %>%
  dplyr::left_join(., allSynonyms, by = c("species" = "name_occs")) %>%
  dplyr::left_join(., batNames,    by = c("name_batnames" = "canonical")) %>%
  dplyr::mutate(
    canonical = case_when(
      is.na(name_batnames)  ~ species,
      !is.na(name_batnames) ~ name_batnames
    )
  ) %>%
  dplyr::select( any_of(newTargetCols) )

# Manual fix of weird character.
occs_syn$canonical[occs_syn$canonical == "†Pipistrellus murrayi"] <- "Pipistrellus murrayi"

# Write.
fwrite(occs_syn, file = file.path(wd$bin, "occs_syn.csv"))
