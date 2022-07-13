

# ================================================================================================

# Scrape PubMed for citation counts for all hosts
# n.b. do not include synonyms for now as taxize db lookup runs v slowly

# ================================================================================================

# root dir and dependencies
# dependencies and basedir
setwd("C:/Users/roryj/Documents/PhD/supervision/2022_CiaraDuggan_UCL/zoonotichosts_processed/")
pacman::p_load("dplyr", "magrittr", "stringr", "ggplot2", "inlabru", "INLA", "lemon", "mgcv", "vroom", "RISmed", "taxize", "Hmisc")

# read hosts list
spp = read.csv("./output/zoonotichosts_edgelist_VIRIONplusCLOVER.csv") %>%
  dplyr::filter(HostClass %in% c("mammalia", "aves"))



# # ------------------ 1. Access synonyms for all host species using taxize -----------------
# 
# # run for all species to access synonyms
# findSyns4 = function(spp){
#   
#   # query taxize using itis
#   syns = taxize::synonyms(spp, db="itis", accepted=TRUE)
#   
#   if(nrow(syns[[1]]) == 0){
#     res = data.frame(
#       Host = spp,
#       Synonym = spp
#     )
#   } else{
#     res = data.frame(
#       Host = spp,
#       Synonym = c(spp, tolower(syns[[1]]$syn_name))
#     )
#   }
# 
#   return(res)
#   Sys.sleep(0.25)
# }
# 
# # run and save for each
# for(h in virion$Host[ 1:nrow(virion) ]){
#   e = simpleError("lookup error")
#   sx = tryCatch(findSyns4(h), error = function(e) e)
#   if(class(sx)[1] == "simpleError"){ sx = data.frame(Host = h, Synonym = "lookup fail")}
#   write.csv(sx, paste("./output/host_synonyms/syns_", h, ".csv", sep=""))
# }
# 
# # combine and save
# do.call(rbind.data.frame, lapply(list.files("./output/host_synonyms/perhost/", full.names=TRUE), read.csv)) %>%
#   dplyr::select(-X) %>%
#   write.csv("./output/host_synonyms/hostsyns_VIRION.csv", row.names=FALSE)




# ================= run pubmed scrape for annual publication counts ==============

# # unique host species separated by host synonyms
# spp = read.csv("./data/clover/Clover_v1.0_NBCIreconciled_20201218.csv", stringsAsFactors = FALSE) %>%
#   dplyr::select(Host, HostClass, HostOrder, HostFamily, HostSynonyms) %>%
#   dplyr::filter(!duplicated(Host)) %>%
#   dplyr::arrange(desc(HostClass), HostOrder, HostFamily, Host) %>%
#   tidyr::separate_rows(HostSynonyms, sep=",")

# # read in species and syns
# syns = read.csv("./output/host_synonyms/hostsyns_VIRION.csv")
# syns$Synonym[ syns$Synonym == "lookup fail" ] = syns$Host[ syns$Synonym == "lookup fail" ]
# syns = syns %>% dplyr::rename("HostSynonyms"=Synonym)
# spp = syns



# ======================= function to extract cites for each species ====================

# function to scrape pubmed results
getPubMedCites = function(species_binomial, taxid_search = FALSE, virus_related = FALSE){
  
  # print species name
  print(sprintf("Processing: %s", species_binomial))
  
  # create search term with all synonyms (if including)
  # datx = spp[ spp$Host == species_binomial, ]
  # sppx = Hmisc::capitalize(str_trim(unique(c(datx$Host[1], datx$HostSynonyms[ datx$HostSynonyms != ""]))))
  # sppx = strsplit(sppx, " ")
  # sppx = lapply(sppx, "[", 1:2)
  # #sppx = sppx[ which(unlist(lapply(sppx, length)) < 3) ]
  # sppx = unique(sppx)
  
  # search term with only main binomial
  datx = spp[ spp$Host == species_binomial, ]
  sppx = strsplit(unique(datx$Host), " ")
  
  if (!taxid_search) {
    # Search by name(s)
    createSearchTerm = function(x){ paste("(", x[1], " [TIAB] AND ", x[2], " [TIAB])", sep="") }
    search_term = paste(unlist(lapply(sppx, createSearchTerm)), collapse=" OR ")
    search_db = "pubmed"
    if(virus_related){
      search_term = paste(search_term, "AND (virus [TIAB] OR viral [TIAB])")
    }
      
  } else {
    # Search by taxid
    taxonomy_id = taxize::get_uid(species_binomial, rank_query = "species", 
                                  ask = FALSE, messages = FALSE)
    
    if (attr(taxonomy_id, "multiple_matches"))
      warning("Multiple names matched for ", species_binomial, ", returning NA")
    
    if (attr(taxonomy_id, "pattern_match")) {
      # Don't want partial matches
      warning("Partial match for", species_binomial, ", returning NA")
      taxonomy_id = NA_character_
    }
    
    search_term = paste0("txid", taxonomy_id,"[Organism:exp]") # Papers mentioning this species or any subspecies
    search_db = "PMC"  # This query only works on Pubmed Central
  }
  
  # run searches: firstly to get number of publications, then set this as retmax to ensure all are returned 
  e = simpleError("test error")
  search1 = tryCatch(RISmed::EUtilsSummary(search_term, type="esearch", db=search_db, datetype='pdat', mindate=1930, maxdate=2020), error=function(e) e)
  if(class(search1)[1] == "simpleError"){ return(data.frame(Year=NA, NumPubs=NA, Host = species_binomial, Note="Lookup error")) }
  numpubs_total = RISmed::QueryCount(search1)
  
  # result df
  resx = data.frame(Host = species_binomial, NumPubs = numpubs_total)
  
  # sleep for 0.25 seconds to prevent over-requesting
  Sys.sleep(0.75)
  return(resx)
}



# ============== run scrape and append records to csv ============

# create filenames
output_loc = "./output/"
save_file = paste(output_loc, "hostcitationcounts_pubmed_append.csv", sep="")

# for all host species get time series of publication effort (axis axis == impossible to lookup)
species_vec = unique(spp$Host)
species_vec = species_vec[ !species_vec %in% c("axis axis", "indicator indicator", "canis lupus familiaris")]

# append each new query to csv
for(i in 1:length(species_vec)){

  # run query (5 attempts)
  cat(paste(i, "...", sep=""))
  e = simpleError("test error")
  for(attempt in 1:5){
    
    resx = tryCatch(getPubMedYears(species_vec[i]), error=function(e) e)
    if(class(resx)[1] != "simpleError"){ 
      break
    } else{
      Sys.sleep(0.5)
      resx = tryCatch(getPubMedYears(species_vec[i]), error=function(e) e)
      }
    }
  
  # return lookup error if still throwing errors
  if(class(resx)[1] == "simpleError"){ 
    resx = data.frame(Year=NA, NumPubs=NA, Host = species_vec[i], Note="Lookup error")
  }
  
  # initialise file on first iteration, and then append
  if(i == 1){
    write.csv(resx, save_file, row.names=FALSE)
  } else{
    write.table(resx, save_file, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE, quote=TRUE) # append
  }
}


# -------------- read in scrape and ensure all species are included ------------

# effort cites
effort = read.csv("./output/hostcitationcounts_pubmed.csv")

# any NAs
if(any(is.na(effort$NumPubs))){ print("NAs found") } else{ print("Data complete") }
if(all(species_vec %in% effort$Host)){ print("All hosts included") } else{ print("Hosts missing") }
