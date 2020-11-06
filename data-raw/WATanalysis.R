#load count matrix and setup data

count_matrix <- count_matrix_assembly("Human_RNAseq_counts.xlsx")

setup <- load_metadata("setup_data.xlsx")
design <- Generate_design_matrix(setup)

#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeAComp = NR_A - Placebo_A,
  NR_effect = NR_B - NR_A,
  Placebo_effect = Placebo_B-Placebo_A,
  NR_vs_placebo_B = NR_B-Placebo_B,
  Interaction = (NR_B - NR_A) - (Placebo_B-Placebo_A),
  levels = design)

dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)
write.xlsx(dgeResults, file = here("data/edgeR.xlsx"), asTable = TRUE)

goResult <- goAnalysis(dgeResults)

printGOterms(goResult)



