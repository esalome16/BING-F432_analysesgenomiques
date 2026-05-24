#génomique

# ajouter une colonne dns CSV pour faire correspondre fas et CSV
# Installer les packages si besoin
if (!require(seqinr)) install.packages("seqinr")

library(seqinr)


setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque/doc pour beauti ")#pour donner la direction

# Charger le FASTA, fait donc plus besoin de le run 
#sequences <- read.fasta("EBOV_genomic_analyses_simulated_dataset_3-1.fas", seqtype = "DNA")

# Vérifier les noms
#names(sequences)  # tu verras : H-4, H-2, H-3, etc.

# Renommer : extraire juste le numéro après H-
#new_names = gsub("H-", "", names(sequences))
#names(sequences) = new_names

# Sauvegarder
#write.fasta(sequences, names = names(sequences),  file.out = "EBOV_OK.fas")

# Vérifier
#names(sequences)

#supp les lignes pas nécessaire dans data
library(readr)
library(dplyr)
library(stringr)
# Lire les fichiers
 csv <- read_csv("EBOV_genomic_analyses_simulated_dataset_3-1.csv")
 fas <- readLines("EBOV_OK.fas")
# Extraire les noms des séquences du FASTA
taxa <- fas[startsWith(fas, ">")] |>
     str_remove("^>") |>
     str_trim()
# Garder uniquement les lignes du CSV présentes dans le FASTA
csv_filtre <- csv |>
     mutate(sequence_name = as.character(sequence_name)) |>
      filter(sequence_name %in% taxa)
# Exporter
 write_csv(csv_filtre, "EBOV_OK.csv")
 
length(taxa)
nrow(csv_filtre)
setdiff(taxa, csv_filtre$sequence_name)
setdiff(csv_filtre$sequence_name, taxa)
 
 # juste date: supprimer lon et lat et convertir en TSV

 data <- read.csv("EBOV_OK.csv")
 # Supprimer les colonnes longitude et latitude
data <- data[, !(names(data) %in% c("longitude", "latitude"))]
# Sauvegarder le nouveau fichier CSV
write_tsv(data, "EBOV_OK_date.tsv")

 
 # same pour juste lon et lat 
traits <- read_csv("EBOV_OK.csv")
# Créer un fichier avec uniquement taxon, longitude, latitude
geo_traits <- traits %>%
   select(taxon = sequence_name, longitude, latitude)
 # Exporter en TSV pour BEAUti
  write_tsv(geo_traits, "EBOV_OK_loc.tsv")
  