library(fs)
library(vroom)
library(Biobase)
library(here)
library(GEOquery)
library(tidyverse)
library(fs)
library(vroom)
library(Biobase)
library(here)
library(magrittr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)
library(openxlsx)
library(pheatmap)

#' count_matrix_loader
#'
#' @param file_type your raw data file type (presently optimized for txt)
#'
#' @return a count matrix containing raw count values for the experiment with group annotations.

count_matrix_assembly <- function(file_type){
  count_file <- fs::dir_ls(here("data-raw/"),
                         regexp = file_type,
                         recurse = TRUE)
  count_matrix <- openxlsx::read.xlsx(xlsxFile = count_file,sheet = 1)
  rownames(count_matrix) <- count_matrix$Gene.name
  count_matrix <- count_matrix %>%
    dplyr::select(-Gene.name)

  return(count_matrix)
}

#' Load metadata and sort them according to count matrix input
#'
#' @param file_name the name of teh metadata file (default "metadata.csv")
#'
#' @return metadata file sorted according to count matrix order

load_metadata <- function(file_name) {
  data_file <- fs::dir_ls(here::here("data-raw/"),
                          regexp = file_name,
                          recurse = T)
  metadata <- openxlsx::read.xlsx(xlsxFile = data_file)
  return(metadata)
}

#' Generate design matrix
#'
#' @param metadata a metadata object generated through the load_metadata function
#'
#' @return a design matrix file


Generate_design_matrix <- function(metadata){
  design <- stats::model.matrix( ~0+Group, setup)
    colnames(design) <-
    stringr::str_remove_all(colnames(design), "\\(|\\)|Group|:")
  return(design)
}

RNAseq_processing <- function(count_matrix, metadata, design, ctrsts) {
  group <- as.matrix(metadata[4])
  RNAseq <- edgeR::DGEList(counts = count_matrix, group = group)
  keep <- edgeR::filterByExpr(RNAseq)
  RNAseq <- RNAseq[keep, , keep.lib.sizes = F]
  RNAseq <- edgeR::calcNormFactors(RNAseq)
  RNAseq <- edgeR::estimateDisp(RNAseq,design, robust = T)
  et <- edgeR::exactTest(RNAseq)
  efit <- edgeR::glmQLFit(RNAseq, design, robust = T)
  dgeResults <- apply(ctrsts, 2, . %>%
                        edgeR::glmQLFTest(glmfit = efit, contrast = .) %>%
                        edgeR::topTags(n = Inf, p.value = 1) %>%
                        magrittr::extract2("table") %>%
                        data.table::as.data.table(keep.rownames = TRUE))
  return(dgeResults)
}


goAnalysis <- function(result_list){
  bg <- result_list[[1]]
  bg_list <- clusterProfiler::bitr(
    bg$rn,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db",
    drop = T
  )

  goResult_list <- vector(mode = "list", length = length(result_list))
  for(i in 1:length(result_list)){
    sig_list<- result_list[[i]] %>%
      dplyr::filter(FDR<0.05)

    eg <- clusterProfiler::bitr(
      sig_list$rn,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = "org.Hs.eg.db",
      drop = T
    )
    goResults <- clusterProfiler::enrichGO(gene = eg$ENTREZID,
                                           universe = bg_list$ENTREZID,
                                           OrgDb = org.Hs.eg.db,
                                           ont = "BP")
    goResult_list[[i]]<- goResults
  }
  for (i in 1:length(goResult_list)){
    names(goResult_list)[i]<-names(result_list)[i]
  }
  return(goResult_list)

}

#' File exporter - exports GOresults as an excel sheet, and prints dotplot and cnet plots
#'
#' @param goList a list object containing one or more enrichResults
#'
#' @return

printGOterms <- function(goList){
  goSheets<- vector(mode = "list", length = length(goList))
  for (i in 1:length(goSheets)){
    goSheets[[i]] <- goList[[i]]@result
    names(goSheets)[i]<-names(goList)[i]
  }
  openxlsx::write.xlsx(goSheets, file = here("data/NASH_NAFLD_GOterms.xlsx"), asTable = TRUE)
  dir.create(here("data/figures"), showWarnings = F)
  for (i in 1:length(goList)){
    dotplot <- enrichplot::dotplot(goList[[i]], title = names(goList)[i],size = 1)
    ggplot2::ggsave(dotplot, filename = paste(here("data/figures"),"/dotplot_",names(goList[i]),".png", sep = ""),width = 12, height = 12, units = "cm", scale = 2.5)
    goList_anno <- clusterProfiler::setReadable(goList[[i]], OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    cnetplot <- enrichplot::cnetplot(goList_anno, title = names(goList)[i], size = 1)
    ggplot2::ggsave(cnetplot, filename = paste(here("data/figures"),"/cnetplot_",names(goList[i]),".png", sep = ""),scale = 2.5)
  }
}


