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

 #load in sample info
 sample_info <- openxlsx::read.xlsx(here("data/sample_info.xlsx"))

 #extract genes in term
candidates_test <- unlist(strsplit(reactomedata[[1]]$genesInTerm[1], ", "))
candidates_test <- as_tibble(candidates_test)

#select candidates from extracted
candidate_data <- cpm_data %>%
  filter(SYMBOL %in% candidates_test$value)
candidate_data <- candidate_data %>% distinct(SYMBOL, .keep_all = T)

rownames(candidate_data)<-candidate_data$SYMBOL
candidate_data <- candidate_data %>%
  dplyr::select(-c(rn, SYMBOL))

#select ID and group


roworder_test <- c("CR_B", "CR_A", "PR_B", "PR_A")
 candidate_groups <- sample_info %>%
   dplyr::select(Sample.Name, Group, Sample_ID) %>%
   dplyr::arrange(Group = factor(Group, levels = roworder_test))

 colorder <- candidate_groups$Sample_ID




candidate_data <- candidate_data %>% relocate(colorder)

colnames(candidate_data)==candidate_groups$Sample_ID
candidate_data <- as.matrix(candidate_data)
candidate_groups <- candidate_groups %>%
  dplyr::select(-Sample_ID)

rownames(candidate_groups) <- candidate_groups$Sample.Name
colnames(candidate_data)<-candidate_groups$Sample.Name
candidate_groups <- candidate_groups %>%
  dplyr::select(-Sample.Name)
heatmap <- pheatmap(candidate_data,
                    treeheight_col = 0,
                    treeheight_row = 0,
                    scale = "row",
                    legend = T,
                    na_col = "white",
                    Colv = NA,
                    na.rm = T,
                    cluster_cols = F,
                    fontsize_row = 5,
                    fontsize_col = 7,
                    cellwidth = 7,
                    cellheight = 5,
                    main = "Reactome Genes",
                    annotation_col = candidate_groups
)

