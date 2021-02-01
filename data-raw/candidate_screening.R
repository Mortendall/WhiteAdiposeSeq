#load in reactomedata
 reactomedata <- list(TimeBComp = NA,
                      PR_effect = NA,
                      CR_effect = NA,
                      PR_A_vs_CR_A = NA,
                      Interaction = NA)
  for (i in 1:5){
    reactomedata[[i]] <- openxlsx::read.xlsx(here("data/210201Reactome_data_annotated.xlsx"), sheet = i)
  }

 #load in cpm data
 cpm_data <- openxlsx::read.xlsx(here("data/CPM_matrix.xlsx"))


candidates_test <- unlist(strsplit(reactomedata[[1]]$genesInTerm[1], ", "))
candidates_test <- as_tibble(candidates_test)

candidate_data <- cpm_data %>%
  filter(SYMBOL %in% candidates_test$value)


