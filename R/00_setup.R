
library(magrittr)
library(dplyr)
library(terra)
library(ggplot2)
library(data.table)
library(sf)

my_dir_path <- getwd()
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$fun     <- file.path( my_dir_path, "functions" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$coords  <- file.path( wd$bin, "coords" )
wd$coords2 <- file.path( wd$bin, "coords2" )
wd$preds   <- file.path( my_dir_path, "preds" )
wd$data    <- file.path( my_dir_path, "data" )
wd$occs    <- file.path( wd$data, "occurrences" )
wd$cleaned <- file.path( wd$data, "occurrences_clean" )
wd$names   <- file.path( wd$data, "batNames")
wd$figs    <- file.path( my_dir_path, "figs" )
wd$out     <- file.path( my_dir_path, "out" )
invisible({ lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) ) })

proj4_eqe <- "+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs"
proj4_wgs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
