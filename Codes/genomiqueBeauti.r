#génomique:préparation des données pour BEAUti

setwd("/Users/salomeembise/Desktop/épidémio/TP/Binome_03/génomique/ProjetBinome3-Arbovirus-Risque/doc pour beauti ")
#supp H-4
library(seqinr)
sequences=read.fasta("EBOV_genomic_analyses_simulated_dataset_3-1.fas", seqtype = "DNA")
new_names=gsub("H-", "", names(sequences))
names(sequences)=new_names
#sauvegarde
write.fasta(sequences, names = names(sequences),  file.out = "EBOV_OK.fas")
#verif
names(sequences)

#supp les lignes pas nécessaire dans data
library(readr)
library(dplyr)
library(stringr)
csv=read_csv("EBOV_genomic_analyses_simulated_dataset_3-1.csv")
fas=readLines("EBOV_OK.fas")
#extraire les noms des séquences du fasta
taxa <- fas[startsWith(fas, ">")] |>
     str_remove("^>") |>
     str_trim()
#garder lignes du CSV pstes dans fasta
csv_ok <- csv |>
     mutate(sequence_name = as.character(sequence_name)) |>
      filter(sequence_name %in% taxa)
#sauvegarde
write_csv(csv_ok, "EBOV_OK.csv")
 
# juste date: supprimer lon et lat et convertir en TSV
data=read.csv("EBOV_OK.csv")
data=data[, !(names(data) %in% c("longitude", "latitude"))]
write_tsv(data, "EBOV_OK_date.tsv")

# same pour juste lon et lat 
traits=read_csv("EBOV_OK.csv")
geo_traits=traits %>%
   select(taxon = sequence_name, longitude, latitude)
#exporter en tsv pour beauti
write_tsv(geo_traits, "EBOV_OK_loc.tsv")
  