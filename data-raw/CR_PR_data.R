
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
#Further talks suggest that 79 and 80 could be excluded as well

setup <- setup %>%
  dplyr::filter(
    !Sample_ID == "025_76" &
      !Sample_ID == "025_77" &
      !Sample_ID == "025_78" &
      !Sample_ID == "025_96" &
      !Sample_ID == "025_79" &
      !Sample_ID == "025_80"
  )
design <- Generate_design_matrix(setup)

count_matrix <- count_matrix %>%
  dplyr::select(-"025_76", -"025_77", -"025_78", -"025_96", -"025_79", -"025_80")
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

#Quality_control_plots(count_matrix, setup)

all(colnames(count_matrix)==setup$Sample_ID)


dgeResults <- RNAseq_processing(count_matrix, setup, design, ctrsts)

annotated <- dgeResults

  for (i in 1:length(annotated)){
  key <- bitr(annotated[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  annotated[[i]] <- left_join(annotated[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
      }

#write.xlsx(annotated, file = here("data/edgeR_PR_CR_realign.xlsx"), asTable = TRUE)

#a bit troubling distribution of P-values suggested to Lars that an analysis through Limma could perhaps cause a better result

group <- as.matrix(setup[4])
RNAseq_limma <- edgeR::DGEList(counts = count_matrix, group = group)
keep <- edgeR::filterByExpr(RNAseq_limma, design = design)
RNAseq_limma <- RNAseq_limma[keep, , keep.lib.sizes = F]

count_matrix_limma <- voomWithQualityWeights(RNAseq_limma$counts, design = design, plot = T)

fit <- lmFit(count_matrix_limma, design = design, method = "robust")
fit2 <- contrasts.fit(fit, ctrsts)
fit2 <- eBayes(fit2, trend = T, robust = T)
resultTable <- list(TimeBComp = topTable(fit2, coef = "TimeBComp", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_effect = topTable(fit2, coef = "PR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    CR_effect = topTable(fit2, coef = "CR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_A_vs_CR_A = topTable(fit2, coef = "PR_A_vs_CR_A", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE))


resultTable_export <- resultTable


for (i in 1:length(resultTable_export)){
  key <- bitr(resultTable_export[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  resultTable_export[[i]] <- left_join(resultTable_export[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
}

#write.xlsx(resultTable_export, file = here("data/edgeR_PR_CR_limma.xlsx"), asTable = TRUE)

#####modelling for info on participants#####

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

ID_key_file <- fs::dir_ls(here::here("data-raw/"),
                                  regexp = "ID_Person",
ID_key <- openxlsx::read.xlsx(xlsxFile = ID_key_file) %>%
                                  recurse = T)
  dplyr::select(-Sample)

setup_person <- left_join(setup, ID_key, by = c("Sample_ID" = "ID"))

#Quality_control_plots(count_matrix, setup)

#QC shows sample 76, 77, 78 and 96 looks strange. These are the same samples that we previously found to look odd(PR1 and CR5L, CR10L)
#Further talks suggest that 79 and 80 could be excluded as well

setup_person <- setup_person %>%
  dplyr::filter(
    !Sample_ID == "025_76" &
      !Sample_ID == "025_77" &
      !Sample_ID == "025_78" &
      !Sample_ID == "025_96" &
      !Sample_ID == "025_79" &
      !Sample_ID == "025_80"
  )


colnames(setup_person)[5]<-"ID"

design <- Generate_design_matrix_with_patient(setup_person)

count_matrix <- count_matrix %>%
  dplyr::select(-"025_76", -"025_77", -"025_78", -"025_96", -"025_79", -"025_80")
all(colnames(count_matrix)==setup_person$Sample_ID)

#select what groups from design yuo wish to compare
ctrsts <- makeContrasts(
  TimeBComp = CR_B - PR_B,
  PR_effect = PR_A - PR_B,
  CR_effect = CR_A - CR_B,
  PR_A_vs_CR_A = PR_A - CR_A,
  Interaction = (PR_A - PR_B) - (CR_A-CR_B),
  levels = design)

rownames(count_matrix) <- rownames(count_matrix_raw)

#Quality_control_plots(count_matrix, setup)

all(colnames(count_matrix)==setup_person$Sample_ID)


#Attempt to analyze data with voom with qualityweights + duplicateCorrelation as pr https://support.bioconductor.org/p/59700/

group <- as.matrix(setup[4])
RNAseq_limma <- edgeR::DGEList(counts = count_matrix)
RNAseq_limma<- calcNormFactors(RNAseq_limma)
design <-  stats::model.matrix( ~0+Group, setup_person)
colnames(design) <-
  stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")

keep <- edgeR::filterByExpr(RNAseq_limma)
RNAseq_limma <- RNAseq_limma[keep, , keep.lib.sizes = F]

count_matrix_limma <- voomWithQualityWeights(RNAseq_limma, design = design, plot = T)
corfit <- duplicateCorrelation(count_matrix_limma, design,block = setup_person$ID)
count_matrix_limma <- voomWithQualityWeights(RNAseq_limma, design = design, plot = T, block = setup_person$ID, correlation = corfit$consensus.correlation)
fit <- lmFit(count_matrix_limma, design = design, block = setup_person$ID, correlation = corfit$consensus.correlation)



fit2 <- contrasts.fit(fit, ctrsts)
fit2 <- eBayes(fit2)
resultTable <- list(TimeBComp = topTable(fit2, coef = "TimeBComp", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_effect = topTable(fit2, coef = "PR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    CR_effect = topTable(fit2, coef = "CR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_A_vs_CR_A = topTable(fit2, coef = "PR_A_vs_CR_A", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE))


resultTable_export <- resultTable


for (i in 1:length(resultTable_export)){
  key <- bitr(resultTable_export[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  resultTable_export[[i]] <- left_join(resultTable_export[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
}

#write.xlsx(resultTable_export, file = here("data/edgeR_PR_CR_doubleVoom.xlsx"), asTable = TRUE)
#####Clinical data for PR/CR experiment#####
clinical_data <- read.xlsx(here::here("data-raw/collected_clinical_data_trimmed.xlsx"))

roworder <- colnames(count_matrix_test)
#modified data to de-seelct T2DM
clinical_data<- clinical_data %>%
  dplyr::mutate(Sample.ID = factor(Sample.ID, levels = roworder)) %>%
  dplyr::select(-T2DM.diagnostic, -Age)
clinical_data <- clinical_data[match(roworder, clinical_data$Sample.ID),]

clinical_matrix <- clinical_data %>%
  dplyr::select(-Sample.ID, -Gender)

rownames(clinical_matrix)<-clinical_data$Sample.ID
clinical_matrix<- as.matrix(clinical_matrix)



sample_info <- load_metadata("025_Metadata_Human_realign")
sample_info <- sample_info %>%
  dplyr::filter(
    !Sample_ID == "025_76" &
      !Sample_ID == "025_77" &
      !Sample_ID == "025_78" &
      !Sample_ID == "025_96" &
      !Sample_ID == "025_79" &
      !Sample_ID == "025_80"
  )



sample_info <- sample_info %>%
  dplyr::filter(Sample_ID %in% setup_person$Sample_ID)


sample_info <- left_join(sample_info, setup_person, by = "Sample_ID")


test <- filterByExpr(count_matrix, design = design)
count_matrix_test <- count_matrix[test, , keep.lib.sizes = F, keep.rownames = T]

all(sample_info$Sample.ID==colnames(count_matrix_test))
all(colnames(count_matrix_test)==rownames(clinical_matrix))

y <- DGEList(counts = count_matrix_test,group = sample_info$Group)
cpm_matrix <- cpm(y, log = T)
rownames(cpm_matrix)<-rownames(count_matrix_test)

cpm_matrix_wo_batch <-
  removeBatchEffect(
    cpm_matrix,
    batch = sample_info$ID,
    covariates = clinical_matrix
  )

plotMDS(cpm_matrix_wo_batch, labels = NULL,
        dim.plot = c(1, 2),
        pch = 16,
        main = "MDS with batch correction",
        col = as.numeric(as.factor(sample_info$Group))
)

mdsData <- plotMDS(cpm_matrix_wo_batch,
        ndim = 3, plot = F)
mdsData <-
  mdsData$cmdscale.out %>% data.table(keep.rownames = TRUE) %>%
  mutate(ID = rownames(y$samples)) %>%
  dplyr::select(-rn) %>%
  mutate(Group = sample_info$Group)

setnames(mdsData, c("V1","V2", "V3"),c("dim1","dim2", "dim3"))
# plot1 <- ggplot(mdsData, aes(x = dim1, y = dim2, colour = Group)) +
#   geom_point(size = 5) +
#   #geom_label(show.legend = FALSE, size = 5) +
#   theme_bw()
#
#
# plot2<- ggplot(mdsData, aes(x = dim1, y = dim3, colour = Group)) +
#   geom_point(size = 5) +
#   #geom_label(show.legend = FALSE, size = 5) +
#   theme_bw()
# plot3 <- ggplot(mdsData, aes(x = dim2, y = dim3, colour = Group)) +
#   geom_point(size = 5) +
#   #geom_label(show.legend = FALSE, size = 5) +
#   theme_bw()

# pdf(file.path(here("data/figures/QCplots"), "MDSplots_after_model_adjustment.pdf"), width = 4, height = 4)
# oldpar <- par()$mfrow
# par(mfrow = c(1, 1))
# plot(plot1)
# plot(plot2)
# plot(plot3)
# par(mfrow = oldpar)
# dev.off()


View(mdsData)

?plotMDS

plotMDS(cpm_matrix, labels = NULL,
        dim.plot = c(1, 2),
        pch = 16,
        main = "MDS with batch correction gender adjustment",
        col = as.numeric(as.factor(clinical_data$Gender))
)
#MDS plots reveal that gender is a major cause of clustering
group <- as.matrix(setup_person[4])
y <- edgeR::calcNormFactors(y)
setup_person <- setup_person %>%
  dplyr::mutate(Gender = clinical_data$Gender)
design <- stats::model.matrix( ~0+Group+Gender, setup_person)
colnames(design) <-
  stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")

colnames(design) <-
  stringr::str_remove_all(colnames(design), "\\(|\\)|Gender|:")
y <- edgeR::estimateDisp(y,design)

efit <- edgeR::glmQLFit(y, design)
ctrsts <- makeContrasts(
  TimeBComp = CR_B - PR_B,
  PR_effect = PR_A - PR_B,
  CR_effect = CR_A - CR_B,
  PR_A_vs_CR_A = PR_A - CR_A,
  Interaction = (PR_A - PR_B) - (CR_A-CR_B),
  levels = design)
dgeResults <- apply(ctrsts, 2, . %>%
                      edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                      edgeR::topTags(n = Inf, p.value = 1) %>%
                      magrittr::extract2("table") %>%
                      data.table::as.data.table(keep.rownames = TRUE))
#this gives more or less the same result. Try with voom instead


count_matrix_limma <- voomWithQualityWeights(y$counts, design = design, plot = T)
corfit <- duplicateCorrelation(count_matrix_limma, design,block = setup_person$ID)
count_matrix_limma <- voomWithQualityWeights(RNAseq_limma, design = design, plot = T, block = setup_person$ID, correlation = corfit$consensus.correlation)
fit <- lmFit(count_matrix_limma, design = design, block = setup_person$ID, correlation = corfit$consensus.correlation)

fit2 <- contrasts.fit(fit, ctrsts)
fit2 <- eBayes(fit2)
resultTable <- list(TimeBComp = topTable(fit2, coef = "TimeBComp", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_effect = topTable(fit2, coef = "PR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    CR_effect = topTable(fit2, coef = "CR_effect", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    PR_A_vs_CR_A = topTable(fit2, coef = "PR_A_vs_CR_A", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE),
                    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1) %>% data.table(keep.rownames = TRUE))


resultTable_export <- resultTable


for (i in 1:length(resultTable_export)){
  key <- bitr(resultTable_export[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  resultTable_export[[i]] <- left_join(resultTable_export[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
}

limma::plotMA(fit2)
