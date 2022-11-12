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
                            regexp = "ID_Person",recurse = T)
  ID_key <- openxlsx::read.xlsx(xlsxFile = ID_key_file) %>%
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

  #design <- Generate_design_matrix_with_patient(setup_person)

  count_matrix <- count_matrix %>%
    dplyr::select(-"025_76", -"025_77", -"025_78", -"025_96", -"025_79", -"025_80")
  all(colnames(count_matrix)==setup_person$Sample_ID)



  rownames(count_matrix) <- rownames(count_matrix_raw)

  #Quality_control_plots(count_matrix, setup)

  all(colnames(count_matrix)==setup_person$Sample_ID)


  #Attempt to analyze data with voom with qualityweights + duplicateCorrelation as pr https://support.bioconductor.org/p/59700/

  #####Clinical data for PR/CR experiment#####
  clinical_data <- read.xlsx(here::here("data-raw/collected_clinical_data_trimmed.xlsx"))

  roworder <- colnames(count_matrix)

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

setup_person <- setup_person %>%
  dplyr::mutate(Gender = clinical_data$Gender)
design <- stats::model.matrix( ~0+Group+Gender, setup_person)
colnames(design) <-
  stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")

colnames(design) <-
  stringr::str_remove_all(colnames(design), "\\(|\\)|Gender|:")


  rownames(count_matrix) <- rownames(count_matrix_raw)
  y <- DGEList(counts = count_matrix,group = sample_info$Group)

  test <- filterByExpr(y, design = design)
  y<-y[test, , keep.lib.sizes = F]

count_matrix_test <- y$counts



all(sample_info$Sample.ID==colnames(count_matrix_test))
all(colnames(count_matrix_test)==rownames(clinical_matrix))

group <- as.matrix(setup_person[4])
y <- edgeR::calcNormFactors(y)
cpm_matrix <- cpm(y, log = T)

#y <- edgeR::estimateDisp(y,design)

#efit <- edgeR::glmQLFit(y, design)
ctrsts <- makeContrasts(
  TimeBComp = CR_B - PR_B,
  PR_effect = PR_A - PR_B,
  CR_effect = CR_A - CR_B,
  PR_A_vs_CR_A = PR_A - CR_A,
  Interaction = (PR_A - PR_B) - (CR_A-CR_B),
  levels = design)
#dgeResults <- apply(ctrsts, 2, . %>%
# edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
# edgeR::topTags(n = Inf, p.value = 1) %>%
# magrittr::extract2("table") %>%
# data.table::as.data.table(keep.rownames = TRUE))
#this gives more or less the same result. Try with voom instead


count_matrix_limma <-
  voomWithQualityWeights(y$counts, design = design, plot = T)
corfit <-
  duplicateCorrelation(count_matrix_limma, design, block = setup_person$ID)
count_matrix_limma <-
  voomWithQualityWeights(
    y$counts,
    design = design,
    plot = F,
    block = setup_person$ID,
    correlation = corfit$consensus.correlation
  )

fit <-
  lmFit(
    count_matrix_limma,
    design = design,
    block = setup_person$ID,
    correlation = corfit$consensus.correlation
  )

fit2 <- contrasts.fit(fit, ctrsts)
fit2 <- eBayes(fit2)
resultTable <- list(TimeBComp = topTable(fit2, coef = "TimeBComp", number = Inf, p.value = 1,sort.by = "p") %>% data.table(keep.rownames = TRUE),
                    PR_effect = topTable(fit2, coef = "PR_effect", number = Inf, p.value = 1,sort.by = "p") %>% data.table(keep.rownames = TRUE),
                    CR_effect = topTable(fit2, coef = "CR_effect", number = Inf, p.value = 1,sort.by = "p") %>% data.table(keep.rownames = TRUE),
                    PR_A_vs_CR_A = topTable(fit2, coef = "PR_A_vs_CR_A", number = Inf, p.value = 1,sort.by = "p") %>% data.table(keep.rownames = TRUE),
                    Interaction = topTable(fit2, coef = "Interaction", number = Inf, p.value = 1, sort.by = "p") %>% data.table(keep.rownames = TRUE))


resultTable_export <- resultTable



for (i in 1:length(resultTable_export)){
  key <- bitr(resultTable_export[[i]]$rn, fromType='ENSEMBL', toType='SYMBOL', OrgDb = "org.Hs.eg.db", drop = F)
  resultTable_export[[i]] <- left_join(resultTable_export[[i]], key, by = c("rn"="ENSEMBL"), keep = F)
}


#write.xlsx(resultTable_export, file = here::here("data/edgeR_PR_CR_201812_sort.xlsx"), asTable = TRUE)

#####Data analysis with reactome####
#rLst <- fread("https://reactome.org/download/current/Ensembl2Reactome.txt", header = FALSE)
rLst <- openxlsx::read.xlsx("C:/Users/tvb217/Documents/R/tmp/Reactomedatabase.xlsx")
rLst <- as.data.table(rLst, keep.rownames = T)
rLst <- rLst[V6 == "Homo sapiens"]
reactomeName <- rLst[, .(ID = unique(V2), TERM = unique(V4))]

reactomeList <- tapply(rLst$V1, rLst$V2, list)
reactomeList <- Filter(. %>% length %>% is_greater_than(4), reactomeList) # Remove small categories

reactomeList <- Filter(. %>% length %>% is_less_than(501), reactomeList) # Remove small categories

camera_test <- apply(ctrsts, 2, camera, index = reactomeList, y = count_matrix_limma, design = design)
camera_test <- lapply(camera_test, data.table, keep.rownames = T)
camera_test <- lapply(camera_test, setnames, old = "rn", new = "ID")
camera_test <- lapply(camera_test, extract, !is.na(PValue))
camera_test <- lapply(camera_test, extract, reactomeName, on = "ID", nomatch = FALSE)
camera_test <- lapply(camera_test, extract, order(PValue, decreasing = FALSE))

#### Gene Ontology (only BP)
keysGO <- keys(GO.db)
termGO <- AnnotationDbi::select(GO.db, keys=keysGO, columns=c("TERM", "ONTOLOGY")) %>% data.table
termGO <- termGO[ONTOLOGY == "BP"]

termGO[, ONTOLOGY:=NULL]
setnames(termGO, "GOID", "ID")

cyt.go.genes <- as.list(org.Hs.egGO2ALLEGS)
cyt.go.genes <- cyt.go.genes[names(cyt.go.genes) %in% termGO$ID]
cyt.go.genes <- Filter(. %>% length %>% is_greater_than(4), cyt.go.genes) # Remove small categories

cyt.go.genes <- Filter(. %>% length %>% is_less_than(501), cyt.go.genes) # Remove large categories

entrez_matrix <- count_matrix_limma
conv <- bitr(rownames(entrez_matrix$E), fromType='ENSEMBL', toType='ENTREZID', OrgDb = "org.Hs.eg.db") %>%
  data.table(key = "ENSEMBL")
rownames(entrez_matrix$E) <- conv[rownames(entrez_matrix$E), ENTREZID, mult = "first"]

GO_test <- apply(ctrsts, 2, camera, index = cyt.go.genes, y = entrez_matrix, design = design)


GO_test <- lapply(GO_test, data.table, keep.rownames = T)
GO_test <- lapply(GO_test, setnames, old = "rn", new = "ID")
GO_test <- lapply(GO_test, extract, !is.na(PValue))
GO_test <- lapply(GO_test, extract, termGO, on = "ID", nomatch = FALSE)
GO_test <- lapply(GO_test, extract, order(PValue, decreasing = FALSE))


#write.xlsx(camera_test, here::here("data/210201Reactome_data.xlsx"), asTable = TRUE)
#write.xlsx(GO_test, here::here("data/210201GO_data.xlsx"), asTable = TRUE)

####GO-term extraction#####
#code from Lars to extract genes
annotateWithGenes <- function(tests, termList, fromType){
  conv <- bitr(unique(unlist(termList)), fromType = fromType, toType = 'SYMBOL', OrgDb = "org.Hs.eg.db") %>%
    data.table(key = fromType)
  pasteSymbols <- function(x) conv[termList[[x]], paste0(SYMBOL, collapse = ", ")]
  termSymbols <- lapply(names(termList), pasteSymbols) %>% data.table
  termSymbols[, id:=names(termList)]
  setnames(termSymbols, c("symbolList", "id"))
  setkey(termSymbols, id)

  tests <- lapply(tests, copy)
  for (i in names(tests)){
    tests[[i]][, genesInTerm:=termSymbols[ID, symbolList]]
  }
  tests
}

cameraGoAnnotated <- annotateWithGenes(GO_test, cyt.go.genes, fromType = 'ENTREZID')
cameraReactomeAnnotated <- annotateWithGenes(camera_test, reactomeList, fromType = "ENSEMBL")

# write.xlsx(cameraReactomeAnnotated, here::here("data/210201Reactome_data_annotated.xlsx"), asTable = TRUE)
# write.xlsx(cameraGoAnnotated , here::here("data/210201GO_data_annotated.xlsx"), asTable = TRUE)

#####prepare volcano plots with FDR threshold####
#PR volcano plot
resultTable_export <- list(TimeBComp = NA,
                           PR_effect = NA,
                           CR_effect = NA,
                           PR_A_vs_CR_A = NA,
                           Interaction = NA)
for (i in 1:length(resultTable_export)){
  resultTable_export[[i]]<-openxlsx::read.xlsx(here::here("data/edgeR_PR_CR_201812_sort.xlsx"),sheet = i)
}
resultTable_export[[2]]<-resultTable_export[[2]] %>%
  dplyr::mutate(sig_group = case_when(
    logFC>1 & -log10(P.Value)>1.30 ~"High",
    logFC< -1 & -log10(P.Value)>1.30 ~"Low",
    -log10(P.Value)<1.30 ~"nonSig",
    TRUE ~"Mid"
  ),
  SYMBOL = case_when(
    is.na(SYMBOL)~rn,
    TRUE ~SYMBOL
  ))
colors <- c("High"= "Blue", "Low" = "red", "Mid" = "brown", "nonSig"="Grey")
PR_volc <- ggplot2::ggplot(resultTable_export[[2]], aes(x = logFC, y = -log10(P.Value), color = sig_group))+
  ggplot2::geom_point(size = 1)+
  ggplot2::ggtitle("Protein restriction treatment")+
  ggplot2::theme(plot.title = element_text(hjust = 0.5,
                                           size = 14),
                 legend.position = "none"
                 )+
  ggplot2::geom_text(aes(y = 5.8, label = "FDR Threshold", x = 0.5), color = "black")+
  ggplot2::geom_hline(yintercept =  5.560481, linetype = "dashed")+
  ggplot2::xlim(c(-2.5,2.5))+
  ggplot2::ylim(0,6.5)+
  ggplot2::scale_color_manual(values =  colors)+
  ggrepel::geom_text_repel(
    data = subset(resultTable_export[[2]], logFC < -1.8 &
                    -log10(P.Value) > 1.3),
    aes(logFC,-log10(P.Value), label = SYMBOL),
    color = "black",
    fontface = "bold",
    max.overlaps = 20,
    arrow = arrow(length = unit(0.02, "npc")),
    box.padding = 1
  )+
  ggrepel::geom_text_repel(
    data = subset(resultTable_export[[2]], logFC > 1.8 &
                    -log10(P.Value) > 1.3),
    aes(logFC,-log10(P.Value), label = SYMBOL),
    color = "black",
    fontface = "bold",
    max.overlaps = 20,
    arrow = arrow(length = unit(0.02, "npc")),
    box.padding = 1
  )+
  ggrepel::geom_text_repel(
    data = subset(resultTable_export[[2]],
                    -log10(P.Value) > 6),
    aes(logFC,-log10(P.Value), label = SYMBOL),
    color = "black",
    fontface = "bold"
  )

resultTable_export[[3]]<-resultTable_export[[3]] %>%
  dplyr::mutate(sig_group = case_when(
    logFC>1 & -log10(P.Value)>1.30 ~"High",
    logFC< -1 & -log10(P.Value)>1.30 ~"Low",
    -log10(P.Value)<1.30 ~"nonSig",
    TRUE ~"Mid"
  ),
  SYMBOL = case_when(
    is.na(SYMBOL)~rn,
    TRUE ~SYMBOL
  )
  )


CR_volc <- ggplot2::ggplot(resultTable_export[[3]], aes(x = logFC, y = -log10(P.Value), color = sig_group))+
  ggplot2::geom_point(size = 1)+
  ggplot2::ggtitle("Caloric restriction treatment")+
  ggplot2::theme(plot.title = element_text(hjust = 0.5,
                                           size = 14),
                 legend.position = "none")+
  ggplot2::geom_text(aes(y = 5.8, label = "FDR Threshold", x = 0), color = "black")+
  ggplot2::geom_hline(yintercept =  5.560481, linetype = "dashed")+
  ggplot2::xlim(c(-2.5,2.5))+
  ggplot2::ylim(0,6.5)+
  ggplot2::scale_color_manual(values =  colors)+
  ggrepel::geom_text_repel(
    data = subset(resultTable_export[[3]], logFC < -1.5 &
                    -log10(P.Value) > 1.3),
    aes(logFC,-log10(P.Value), label = SYMBOL),
    color = "black",
    fontface = "bold",
    max.overlaps = 20,
    arrow = arrow(length = unit(0.02, "npc")),
    box.padding = 1

  )+
  ggrepel::geom_text_repel(
    data = subset(resultTable_export[[3]], logFC > 1.8 &
                    -log10(P.Value) > 1.3),
    aes(logFC,-log10(P.Value), label = SYMBOL),
    color = "black",
    fontface = "bold",
    max.overlaps = 20,
    arrow = arrow(length = unit(0.02, "npc")),
    box.padding = 1

  )

PR_volc+CR_volc

tiff(here::here("data/figures_PR_CR/Volcanotplots.tif"), res = 300, height = 20, width = 30, units = "cm")
PR_volc+CR_volc
dev.off()

#####Create heatmaps#####

#To run this code, first run the data processing in the top of the script.


all(colnames(cpm_matrix)==setup_person$Sample_ID)
#create heatmap for PR
setup_person_heatmap_PR <- setup_person %>%
  dplyr::filter(Group == "PR_A"|Group == "PR_B") %>%
  dplyr::arrange(desc(Condition1),ID)

cpm_matrix_PR <- as.data.frame(cpm_matrix) %>%
  dplyr::select(setup_person_heatmap_PR$Sample_ID)
all(colnames(cpm_matrix_PR)==setup_person_heatmap_PR$Sample_ID)

PR_key <- setup_person_heatmap_PR
rownames(PR_key) <- PR_key$Sample_ID
PR_key$ID <- stringr::str_remove_all(PR_key$ID, "Person_")
PR_key <- PR_key %>%   dplyr::select(ID, Condition1)
PR_key$Condition1<-factor(PR_key$Condition1, c("B","A"))
colnames(PR_key)[2]<-"Treatment"
cpm_matrix_PR<-as.matrix(cpm_matrix_PR)

PR_hm <- pheatmap::pheatmap(cpm_matrix_PR,
                            # treeheight_col = 0,
                            # treeheight_row = 0,
                            scale = "row",
                            legend = T,
                            na_col = "white",
                            Colv = NA,
                            na.rm = T,
                            cluster_cols = F,
                            cluster_rows = F,
                            fontsize_row = 8,
                            fontsize_col = 14,
                            cellwidth = 16,
                            cellheight = 0.025,
                            annotation_col = PR_key,
                            show_rownames = F,
                            show_colnames = F,
                            main = "Protein Restriction"
                            )
#make the plot for CR
setup_person_heatmap_CR <- setup_person %>%
  dplyr::filter(Group == "CR_A"|Group == "CR_B") %>%
  dplyr::arrange(desc(Condition1),ID)

cpm_matrix_CR <- as.data.frame(cpm_matrix) %>%
  dplyr::select(setup_person_heatmap_CR$Sample_ID)
all(colnames(cpm_matrix_CR)==setup_person_heatmap_CR$Sample_ID)
colnames()
CR_key <- setup_person_heatmap_CR
rownames(CR_key) <- CR_key$Sample_ID
#correct typo in CR_key
CR_key$ID <- case_when(CR_key$ID == "Percon_8C"~"Person_8C",
                       TRUE ~ as.character(CR_key$ID))
CR_key$ID <- stringr::str_remove_all(CR_key$ID, "Person_")
CR_key <- CR_key %>%   dplyr::select(ID, Condition1)
CR_key$Condition1<-factor(CR_key$Condition1, c("B","A"))
colnames(CR_key)[2]<-"Treatment"
cpm_matrix_CR<-as.matrix(cpm_matrix_CR)

CR_hm <- pheatmap::pheatmap(cpm_matrix_CR,
                            # treeheight_col = 0,
                            # treeheight_row = 0,
                            scale = "row",
                            legend = T,
                            na_col = "white",
                            Colv = NA,
                            na.rm = T,
                            cluster_cols = F,
                            cluster_rows = T,
                            fontsize_row = 8,
                            fontsize_col = 14,
                            cellwidth = 16,
                            cellheight = 0.025,
                            annotation_col = CR_key,
                            show_rownames = F,
                            show_colnames = F,
                            main = "Caloric Restriction",
)
PR_hm <- ggplotify::as.ggplot(PR_hm)
CR_hm <- ggplotify::as.ggplot(CR_hm)



tiff(here::here("data/figures_PR_CR/Heatmaps.tif"), width = 30, height = 30, res = 200, units = "cm", )
PR_hm+CR_hm
dev.off()

#run cor and as.dist goes into hclust

cor()

#####Create MDS plots####
#run the code above to generate the necessary designs and data files
mdsData <- plotMDS(y, ndim = 3, plot = FALSE)
group <- as.matrix(setup_person)



mdsData <-
  mdsData$eigen.vectors %>% as.data.table() %>%
  dplyr::mutate(ID = rownames(y$samples)) %>%
  dplyr::mutate(Group = setup_person$Group,
                Sex = setup_person$Gender,
                Person = setup_person$ID) %>%
  dplyr::select(ID, Group,Sex, Person, V1, V2, V3)


setnames(mdsData,
         c("V1", "V2", "V3", "ID", "Group", "Sex", "Person"),
         c("dim1", "dim2", "dim3", "ID", "Group", "Sex", "Person"))

pBase <-
  ggplot(mdsData, aes(x = dim1, y = dim2, colour = Person)) +
  geom_point(size = 10) +
  #geom_label(show.legend = FALSE, size = 5) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5))+
  ggtitle("MDS Plot")
pBase

