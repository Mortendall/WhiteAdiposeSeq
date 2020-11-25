
#load count matrix and setup data

count_matrix <- count_matrix_assembly("Human_RNAseq_counts_CR.xlsx")

setup <- load_metadata("setup_data_CR.xlsx")


#Quality_control_plots(count_matrix, setup)

#both PR1 samples have a very low number of counts. Same for CR5L and CR10L
setup <- setup %>%
  dplyr::filter(!ID == "PR1" & !ID =="PR1L" & !ID == "CR5L" & !ID == "CR10L")
design <- Generate_design_matrix(setup)

count_matrix <- count_matrix %>%
  dplyr::select(-PR1, -PR1L, -CR5L, -CR10L)


#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeAComp = CR_A - PR_A,
  PR_effect = PR_B - PR_A,
  CR_effect = CR_B - CR_A,
  PR_B_vs_CR_B = PR_B - CR_B,
  Interaction = (PR_B - PR_A) - (CR_B-CR_A),
  levels = design)

all(setup$ID == colnames(count_matrix))


dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)
write.xlsx(dgeResults, file = here("data/edgeR_PR_CR.xlsx"), asTable = TRUE)

#goResult <- goAnalysis(dgeResults)

#printGOterms(goResult)


