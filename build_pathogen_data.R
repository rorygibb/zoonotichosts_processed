

# =============== Reads host-pathogen data from CLOVER and VIRION and turns into a zoonotic host edgelist ================

library(RCurl); library(vroom); library(dplyr); library(magrittr)


# CLOVER files are stored at: https://github.com/viralemergence/clover/blob/main/clover/clover_1.0_allpathogens/

# VIRION files are stored at: https://github.com/viralemergence/virion/tree/main/Virion



# --------- 1. CLOVER --------------

# read in bacteria, fungi and helminth files
clov_loc = "C:/Users/roryj/Documents/PhD/202011_clover/repos/clover/clover/clover_1.0_allpathogens/"
clover = do.call(
  rbind.data.frame,
  lapply(
    list.files(clov_loc, pattern="csv", full.names = TRUE)[ 1:3 ], read.csv
  )
)

# remove viruses as these come from VIRION and filter to key fields
clover = clover %>%
  dplyr::filter(PathogenType != "virus") %>%
  dplyr::filter(HostNCBIResolved == TRUE & PathogenNCBIResolved == TRUE) %>% # keep tax resolved only
  dplyr::select(
    Host, HostClass, HostOrder, HostFamily, HostGenus,
    Pathogen, PathogenType, PathogenClass, PathogenOrder, PathogenFamily, PathogenGenus,
    DetectionMethod
  ) %>%
  unique() %>%
  dplyr::filter(!is.na(Host) & !is.na(Pathogen))


# -------- 2. VIRION -------------

# read in VIRION datafile
virion = vroom::vroom("C:/Users/roryj/Documents/PhD/202104_virion/virion/Virion/Virion.csv.gz")

virion = virion %>%
  dplyr::rename(
    "Pathogen"=Virus, "PathogenClass"=VirusClass, 
    "PathogenOrder"=VirusOrder, "PathogenFamily"=VirusFamily, "PathogenGenus"=VirusGenus, "PathogenNCBIResolved" = VirusNCBIResolved
  ) %>%
  dplyr::mutate(PathogenType = "virus", 
                PredictFlag = stringr::str_detect(Pathogen, "predict\\_")) %>%
  dplyr::filter(ICTVRatified | PredictFlag) %>% # keep only properly taxonomically resolved viruses
  dplyr::filter(HostNCBIResolved == TRUE) %>% # keep tax resolved hosts only
  dplyr::select(
    Host, HostClass, HostOrder, HostFamily, HostGenus,
    Pathogen, PathogenType, PathogenClass, PathogenOrder, PathogenFamily, PathogenGenus,
    DetectionMethod
  ) %>%
  unique() %>%
  dplyr::filter(!is.na(Host) & !is.na(Pathogen))


# ------- Combine into edgelist (host * pathogen * detection method) ----------

# combine
hp = rbind(clover, virion)

# identify "human pathogen" (detected in humans via PCR or isolation)
zoono1 = hp %>% 
  dplyr::filter(Host == "homo sapiens" & DetectionMethod %in% c("Isolation/Observation", "PCR/Sequencing"))

# identify
hp$InfectsHumans = hp$Pathogen %in% zoono1$Pathogen

# exclude humans from list
hp = hp %>%
  dplyr::filter(Host != "homo sapiens")




# ------- Subset to zoonotic hosts list -----------

# broad detection criterion: serology, PCR or isolation
zh1 = hp %>%
  dplyr::filter(DetectionMethod != "Not specified") %>%
  dplyr::select(-DetectionMethod) %>%
  unique() %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(HostClass = head(HostClass, 1),
                   HostOrder = head(HostOrder, 1), 
                   HostFamily = head(HostFamily, 1), 
                   HostGenus = head(HostGenus, 1),
                   ZoonoticHost = any(InfectsHumans == TRUE),
                   PathogenRichness = n_distinct(Pathogen),
                   ZoonoticPathogenRichness = n_distinct(Pathogen[ InfectsHumans == TRUE ]),
                   Pathogens = paste(Pathogen, collapse=", "),
                   ZoonoticPathogens = paste(Pathogen[ InfectsHumans == TRUE ], collapse=", ")) %>%
  dplyr::mutate(HostCriteria = "Broad (serology, sequencing or isolation)")

# narrower detection criterion 
zh2 = hp %>%
  dplyr::filter(DetectionMethod != c("Not specified", "Antibodies")) %>%
  dplyr::select(-DetectionMethod) %>%
  unique() %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(HostClass = head(HostClass, 1),
                   HostOrder = head(HostOrder, 1), 
                   HostFamily = head(HostFamily, 1), 
                   HostGenus = head(HostGenus, 1),
                   ZoonoticHost = any(InfectsHumans == TRUE),
                   PathogenRichness = n_distinct(Pathogen),
                   ZoonoticPathogenRichness = n_distinct(Pathogen[ InfectsHumans == TRUE ]),
                   Pathogens = paste(Pathogen, collapse=", "),
                   ZoonoticPathogens = paste(Pathogen[ InfectsHumans == TRUE ], collapse=", ")) %>%
  dplyr::mutate(HostCriteria = "Stricter (sequencing or isolation)")


# combine
hosts = rbind(zh1, zh2)

# save file
write.csv(hosts, "./output/zoonotichosts_edgelist_VIRIONplusCLOVER.csv", row.names=FALSE)

# metadata
# ZoonoticHost = is species identified as zoonotic host under given criterion
# PathogenRichness = overall pathogen richness (num pathogens)
# ZoonoticPathogenRichness = zoonotic pathogen richness (num pathogens)
# Pathogens and ZoonoticPathogens = names of pathogens
# HostCriteria = describes detection criteria, either broad or stricter, filter to one of these before further processing
