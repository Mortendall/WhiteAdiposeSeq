#load count matrix and setup data

count_matrix <- count_matrix_assembly("Human_RNAseq_counts.xlsx")

setup <- load_metadata("setup_data.xlsx")
design <- Generate_design_matrix(setup)

Quality_control_plots(count_matrix, setup)

#S10 and S42 placebo A look like they are seperaring from the remaining groups. Consider removing them
setup <- setup %>%
  dplyr::filter(!ID == "S42"& !ID == "S10")

count_matrix <- count_matrix %>%
  dplyr::select(-S10,-S42)

#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeAComp = NR_A - Placebo_A,
  NR_effect = NR_B - NR_A,
  Placebo_effect = Placebo_B-Placebo_A,
  NR_vs_placebo_B = NR_B-Placebo_B,
  Interaction = (NR_B - NR_A) - (Placebo_B-Placebo_A),
  levels = design)

dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)
write.xlsx(dgeResults, file = here("data/edgeR_afterQC.xlsx"), asTable = TRUE)

#goResult <- goAnalysis(dgeResults)

printGOterms(goResult)



