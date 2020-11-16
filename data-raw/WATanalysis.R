#load count matrix and setup data

count_matrix <- count_matrix_assembly("Human_RNAseq_counts.xlsx")

setup <- load_metadata("setup_data.xlsx")


#Quality_control_plots(count_matrix, setup)

#S42 placebo A look like it is seperaring from the remaining groups. Consider removing it
setup <- setup %>%
  dplyr::filter(!ID == "S42")
design <- Generate_design_matrix(setup)

count_matrix <- count_matrix %>%
  dplyr::select(-S42)

#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeAComp = NR_A - Placebo_A,
  NR_effect = NR_B - NR_A,
  Placebo_effect = Placebo_B-Placebo_A,
  NR_vs_placebo_B = NR_B-Placebo_B,
  Interaction = (NR_B - NR_A) - (Placebo_B-Placebo_A),
  levels = design)

all(setup$ID == colnames(count_matrix))

dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)
write.xlsx(dgeResults, file = here("data/edgeR_afterQC.xlsx"), asTable = TRUE)

#goResult <- goAnalysis(dgeResults)

printGOterms(goResult)



