## Load packages
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "here",
  "tidyverse",
  "lubridate",
  "dplyr",
  "here"
)

pacman::p_load(pkgs, character.only = T)


## if multiple in put sheets, load and merge sheets
##harry
studies_h = read.csv("../Data/raw_trapping_data/studies_h.csv")
studies_h = studies_h[0:8]
host_h = read.csv("./Data/raw_trapping_data/host_h.csv")
host_h = host_h[0:13]
pathogen_h = read.csv("./Data/raw_trapping_data/pathogen_h.csv")
pathogen_h = pathogen_h[0:12]

##david
studies_d = read.csv("./Data/raw_trapping_data/studies_d.csv")
studies_d = studies_d[0:8]
host_d = read.csv("./Data/raw_trapping_data/host_d.csv")
host_d = host_d[0:13]
pathogen_d = read.csv("./Data/raw_trapping_data/pathogen_d.csv")
pathogen_d = pathogen_d[0:12]

##ana
studies_a = read.csv("./Data/raw_trapping_data/studies_a.csv")
studies_a = studies_a[0:8]
host_a = read.csv("./Data/raw_trapping_data/host_a.csv")
host_a = host_a[0:13]
pathogen_a = read.csv("./Data/raw_trapping_data/pathogen_a.csv")
pathogen_a = pathogen_a[0:12]


studies <-  rbind(studies_h, studies_d, studies_a)
host <-  rbind(host_h, host_d, host_a)
pathogen <-  rbind(pathogen_h, pathogen_d, pathogen_a)



## Filter for David's data


## Filter for site resolution

#studies <- subset(studies, data_resolution == "site-session")

## Filter for site resolution <1kmsq
studies <- subset(studies, kmsq_site == "y")

# keep only host and pathogen data with matching study ID from this subset

host <- host[host$study_id %in% studies$study_id, ]
pathogen <- pathogen[pathogen$study_id %in% studies$study_id, ]

# how many unique sites per species do we have?



#saveRDS(pathogen, file="./data/pathogen_assay_record.rds")

# Keep only 1 pathogen row per rodent
## check for multiple pathogen entries for each rodent_id and take max positives

pathogen <- pathogen %>%
  group_by(associated_rodent_record_id) %>%
  slice(which.max(positive)) %>%
  ungroup()

## Impute pathogen negatives
## get rodent record IDs with no corresponding pathogen entry
host_impute <- host[!(host$rodent_record_id %in% pathogen$associated_rodent_record_id), ]

## Get relevant columns from this dataframe
new_pathogen <- data.frame(host_impute$rodent_record_id,
                           host_impute$study_id,
                           host_impute$scientificName,
                           host_impute$individualCount)

# rename columns to match pathogen sheet
new_pathogen <- new_pathogen %>% 
  rename(
    associated_rodent_record_id = host_impute.rodent_record_id,
    study_id = host_impute.study_id,
    associatedTaxa = host_impute.scientificName,
    tested = host_impute.individualCount
  )

new_pathogen <- merge(new_pathogen, pathogen[, c("study_id", "family", "scientificName", "assay")], by = "study_id", all.x = TRUE)

new_pathogen <- unique(new_pathogen) ### THINK ABOUT NEGATIVE ARENAVIRUS/HANTAVIRUS TESTS WHERE BOTH ARE TESTED BUT ONLY 1 REPORTED

# add pathogen ID, positives and negatives and other (blank) columns
df_length <- length(new_pathogen[,1])

new_pathogen$positive <- rep(0, df_length)
new_pathogen$negative <- new_pathogen$tested
new_pathogen$number_inconclusive <- NA
new_pathogen$note <- NA
new_pathogen$pathogen_record_id <- paste("i", 1:df_length, sep = "_")

# merge pathogen and imputed pathogen dataframes

pathogen <- rbind(pathogen, new_pathogen)



#########

pathogen <- pathogen %>% 
  rename(
    rodent_record_id = associated_rodent_record_id
  )

host <- merge(host,pathogen,by="rodent_record_id")

## convert coordinates to numeric

host$decimalLatitude <- as.numeric(host$decimalLatitude)
host$decimalLongitude <- as.numeric(host$decimalLongitude)

### correct column names and remove excess columns

host <- subset(host, select = -c(verbatimLocality,
                                 pathogen_record_id,
                                 study_id.y,
                                 associatedTaxa,
                                 number_inconclusive,
                                 note,
                                 trapEffort,
                                 trapEffortResolution))

host <- host %>% 
  rename(
    host_name = scientificName.x,
    pathogen_name = scientificName.y,
    study_id = study_id.x
  )

# Harmonise virus names in the dataset
host_path_wide <- host %>% 
  mutate(
    pathogen_abbrev = case_when(
      pathogen_name == "Sin Nombre Virus" ~ "SNV",
      pathogen_name == "Sin Nombre Orthohantavirus" ~ "SNV", 
      pathogen_name == "Sin Nombre virus" ~ "SNV", 
      pathogen_name == "Catarina Virus" ~ "CTNV",
      pathogen_name == "Oliveros Virus" ~ "OLV",
      pathogen_name == "Latino mammarenavirus" ~ "LATV",
      pathogen_name == "Orthohantavirus dobravaense" ~ "DOBV",
      pathogen_name == "Dobrava-Belgrade orthohantavirus" ~ "DOBV",
      pathogen_name == "Mammarenavirus bearense" ~ "BCNV",
      pathogen_name == "Puumala orthohantavirus" ~ "PUUV",
      pathogen_name == "Puumala virus" ~ "PUUV",
      pathogen_name == "Seoul virus" ~ "SEOV",
      pathogen_name == "Guanarito virus" ~ "GTOV",
      pathogen_name == "Lassa mammarenavirus" ~ "LASV",
      pathogen_name == "Pirital virus" ~ "PIRV",
      pathogen_name == "Rat hepatitits E virus" ~ "RHEV",
      pathogen_name == "Hantaan virus" ~ "HTNV",
      pathogen_name == "Hantaviridae" ~ "HTNV",
      pathogen_name == "mammarenavirus" ~ "",
      pathogen_name == "Leptospriosa interrogans" ~ "",
      # Keep all other names the same
      TRUE ~ pathogen_name
    ),
  )

## Harmonise host species names in the dataset

host_path_wide <- host_path_wide %>% 
  mutate(
    host_name = case_when(
      host_name == "Akodon azare" ~ "Akodon azarae",
      host_name == "Clethryonomys glareolus" ~ "Myodes glareolus", 
      host_name == "Spermophilus mexicanus" ~ "Ictidomys mexicanus",
      host_name == "P. californicus" ~ "Peromyscus californicus",
      host_name == "Brush mice" ~ "Peromyscus boylii",
      host_name == "dusky-footed" ~ "Neotoma fuscipes",
      host_name == "R. norvegicus" ~ "Rattus norvegicus",
      host_name == "bank vole" ~ "Myodes glareolus",
      host_name == "Oligoryzomys formesi" ~ "Oligoryzomys fornesi",
      host_name == "Crocidura cf. yaldeni" ~ "Crocidura yaldeni",
      host_name == "Crocidura cf. macmillani" ~ "Crocidura macmillani",
      # Keep all other names the same
      TRUE ~ host_name
    ),
  )

## translate species name into genus/species names

host_path_wide <- host_path_wide %>%
  separate(host_name, into = c("genus", "species", "other"), sep = " ", remove = FALSE)

# save to rds for phylogeny

write_rds(host_path_wide, file="./data/trapping_data_NC.rds")


### how many communities do we have?

dim(unique(host[c('start_date','locality')]))
