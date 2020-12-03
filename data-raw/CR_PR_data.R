
#NB! Recently learned that the B samples are "Before" and "A" are after. Hence, this code should be modified to reflect this

#load count matrix and setup data

count_matrix <- count_matrix_assembly("Human_RNAseq_counts_CR.xlsx")

setup <- load_metadata("setup_data_CR.xlsx")


Quality_control_plots(count_matrix, setup)

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

#####data analysis on re-align#####
count_matrix_raw <- count_matrix_assembly_realign("realign_counts")

setup_raw <- load_metadata_realign("025_Metadata_Human_realign")

setup <- setup_raw %>%
  dplyr::filter(Condition2 == "Calorie restriction" | Condition2 == "Protein restriction")

count_matrix <- count_matrix_raw %>%
  dplyr::select(setup$Sample_ID)

all(setup$Sample_ID==colnames(count_matrix))

setup <- setup %>%
  dplyr::mutate(Condition1 = case_when(Condition1 == "1 - Before" ~ "B",
                                       Condition1 == "2 - After" ~ "A"),
                Condition2 = case_when(Condition2 == "Calorie restriction"  ~ "CR",
                                       Condition2 == "Protein restriction" ~ "PR"))
setup <- setup %>%
  dplyr::mutate(Group = paste(Condition2, Condition1, sep = "_"))

#Quality_control_plots(count_matrix, setup)

#QC shows sample 76, 77, 78 and 96 looks strange. These are the same samples that we previously found to look odd(PR1 and CR5L, CR10L)


setup <- setup %>%
  dplyr::filter(!Sample_ID == "025_76" & !Sample_ID =="025_77" & !Sample_ID == "025_78" & !Sample_ID == "025_96")
design <- Generate_design_matrix(setup)

count_matrix <- count_matrix %>%
  dplyr::select(-"025_76", -"025_77", -"025_78", -"025_96")
all(colnames(count_matrix)==setup$Sample_ID)

#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeBComp = CR_B - PR_B,
  PR_effect = PR_A - PR_B,
  CR_effect = CR_A - CR_B,
  PR_A_vs_CR_A = PR_A - CR_A,
  Interaction = (PR_A - PR_B) - (CR_A-CR_B),
  levels = design)

rownames(count_matrix) <- rownames(count_matrix_raw)

dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)

annotated <- dgeResults

  for (i in 1:length(annotated)){
  key <- bitr(annotated[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  annotated[[i]] <- left_join(annotated[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
      }

write.xlsx(annotated, file = here("data/edgeR_PR_CR_realign.xlsx"), asTable = TRUE)
