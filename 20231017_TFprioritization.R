library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(stringr)

######## Ensembl ##########
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
ensembl_ms <- useEnsembl(biomart = "genes") 
ensembl_ms <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl_ms)

######## Functions #######################
make_count_mat <- function(TF_intersect_mat) {
  TF_count_mat <- data.frame(matrix(ncol = length(unique(TF_intersect_mat$RE_name)), 
                                    nrow = length(unique(TF_intersect_mat$TF_name))))
  row.names(TF_count_mat) <- c(unique(TF_intersect_mat$TF_name))
  colnames(TF_count_mat) <- unique(TF_intersect_mat$RE_name)
  
  for (x in unique(TF_intersect_mat$RE_name)) {
    df <- subset(TF_intersect_mat, RE_name == x)
    for (y in unique(TF_intersect_mat$TF_name)) {
      TF_count_mat[y, x] = sum(df$TF_name == y)
    }
  }
  return(TF_count_mat)
}

find_enrichment <- function(TF_count_mat, low_RE_names, hi_RE_names) {
  TF_count_mat$high_total <- rowSums(TF_count_mat[, hi_RE_names])
  TF_count_mat$low_total <- rowSums(TF_count_mat[, low_RE_names])
  TF_count_mat$hgnc_symbol <- rownames(TF_count_mat)
  TF_count_mat <- TF_count_mat[(TF_count_mat$high_total > TF_count_mat$low_total) & 
                                                   (TF_count_mat$high_total > 1),]
  return(TF_count_mat)
}

# TF_count_mat: requires that row names are gene names in all caps; no duplicates
#     the gene names from ReMap ChIP data often have something like "PBX1-2-3" or just "PBX"
#     meant to illustrate that any PBX TF could bind there, but to correlate expression, need individual entries for each gene
#     so, duplicate the PBX entry to appear as three rows for PBX1, PBX2, and PBX3
# gene_list: list of all gene names in appropriate organism
# returns TF_count_mat
# TO DO: make run automatically from within make_count_mat function
replace_nonspecific_ChIP_TFs <- function(TF_count_mat, gene_list) {
  genes_to_remove <- c()
  for (gene in rownames(TF_count_mat)) {
    if (!(gene %in% gene_list)) {
      loc <- unlist(gregexpr("-", gene))
      # genes are separated by dashes
      if (loc[1] > 0) { 
        base <- substr(gene, 1, (loc[1] - 2))
        new_gene <- substr(gene, 1, loc[1] - 1)
        if (!(new_gene %in% rownames(TF_count_mat))){
          TF_count_mat[new_gene, ] <- TF_count_mat[gene, ]
        }
        for (i in c(1:length(loc))) {
          new_gene <- paste0(base, substr(gene, loc[i] + 1, loc[i] + 1))
          if (!(new_gene %in% rownames(TF_count_mat))){
            TF_count_mat[new_gene, ] <- TF_count_mat[gene, ]
          }
        }
      }
      # only gene base name is provided
      else {
        I_loc <- grep(gene, gene_list)
        for (i in I_loc) {
          new_gene <- gene_list[i]
          if (!(new_gene %in% rownames(TF_count_mat))){
            TF_count_mat[new_gene, ] <- TF_count_mat[gene, ]
          }
        }
      }
      genes_to_remove <- c(genes_to_remove, gene)
    }
  }
  TF_count_mat <- TF_count_mat[!(rownames(TF_count_mat) %in% genes_to_remove), ]
  return(TF_count_mat)
}

# expression_mat: requires rownames are gene symbols
# cols_for_cor: columns to use to correlate gene expression (i.e. specific cell types)
correlation <- function(TFs, gene_interest, expression_mat, cols_for_cor) {
  corr <- data.frame(hgnc_symbol = TFs, correlation = rep(0, length(TFs)))
  rownames(corr) <- corr$hgnc_symbol
  for (TF in TFs){
    a <- expression_mat[TF,][,cols_for_cor]
    b <- expression_mat[gene_interest, ][,cols_for_cor]
    corr[TF, ]$correlation <- cor(t(a), t(b))
  }
  return(corr)
}

rm_zero <- function(y, x){
  x <- x[which(y != 0)]
  x <- x + 0.1
  return(x)
}

#cell types must be names of cells in Immgen dataset
make_dotplot <- function(TF_count, TF_subset, ccRE_names, cell_types, Immgen) {
  celltype_TF_count <- TF_count[TF_subset,]
  
  df <- as.data.frame(matrix(0, nrow = nrow(celltype_TF_count), ncol = 2 * length(cell_types)))
  colnames(df) <- c(paste0(cell_types, ".expr"), paste0(cell_types, ".sum"))
  celltype_TF_count <- cbind(celltype_TF_count, df)
  
  for(x in seq(1:nrow(celltype_TF_count))){
    gene <- celltype_TF_count[x,]$hgnc_symbol
    for(cell in cell_types){
      celltype_TF_count[gene, ][, paste0(cell, ".expr")] <- Immgen[gene, ][, cell]
      if(Immgen[gene,][, cell] > 100){
        celltype_TF_count[gene, ][, paste0(cell, ".sum")] <- 1.0
      }
    }
  }
  
  cell_list <- c()
  for(cell in cell_types){
    cell_list <- c(cell_list, rep(cell, length(ccRE_names)))
  }
  
  dotplot <- data.frame(RE = rep(ccRE_names, length(cell_types)), 
                        sum = rep(0, length(ccRE_names) * length(cell_types)), 
                        expr = rep(0, length(ccRE_names) * length(cell_types)),
                        cell = cell_list)
  dotplot$RE <- factor(dotplot$RE, levels = ccRE_names)
  
  for(x in seq(1:length(ccRE_names))){
    for(cell in cell_types){
      dotplot[dotplot$RE == ccRE_names[x] & dotplot$cell == cell,]$sum <- sum(celltype_TF_count[, x] * celltype_TF_count[, paste0(cell, ".sum")])
      a <- rm_zero(celltype_TF_count[, x], celltype_TF_count[, paste0(cell, ".expr")])
      if(length(a) > 0){
        dotplot[dotplot$RE == ccRE_names[x] & dotplot$cell == cell,]$expr <- exp(mean(log(a)))
      }
    }
  }
  
  dotplot$expr_color <- cut(dotplot$expr, breaks = seq(from = -1, to = max(dotplot$expr + 10), length.out = 101), labels = seq(1:100))
  
  return(dotplot)
}

############### Set up expression datasets from public data #########################

############# Immgen splenic NK Cell IFN-g 2 h treatment #######################
#####Series GSE75306
#[MoGene-1_0-st] Affymetrix Mouse Gene 1.0 ST Array [transcript (gene) version]
NK.SP.2h.IFN.rep1 <- read.delim("Immgen_NK_IFNg_2h_normcounts/NK.IFNg.2hr.Sp.rep1.txt", col.names = c("ID_REF", "NK.IFNg.1"))
NK.SP.2h.IFN.rep2 <- read.delim("Immgen_NK_IFNg_2h_normcounts/NK.IFNg.2hr.Sp.rep2.txt", col.names = c("ID_REF", "NK.IFNg.2"))
NK.SP.rep1 <- read.delim("Immgen_NK_IFNg_2h_normcounts/NK.SP.SPF.rep1.txt", col.names = c("ID_REF", "NK.untx.1"))
NK.SP.rep2 <- read.delim("Immgen_NK_IFNg_2h_normcounts/NK.SP.SPF.rep2.txt", col.names = c("ID_REF", "NK.untx.2"))
NK.SP.rep3 <- read.delim("Immgen_NK_IFNg_2h_normcounts/NK.SP.SPF.rep3.txt", col.names = c("ID_REF", "NK.untx.3"))

NK.Immgen <- left_join(NK.SP.2h.IFN.rep1, NK.SP.2h.IFN.rep2, by = "ID_REF")
NK.Immgen <- left_join(NK.Immgen, NK.SP.rep1, by = "ID_REF")
NK.Immgen <- left_join(NK.Immgen, NK.SP.rep2, by = "ID_REF")
NK.Immgen <- left_join(NK.Immgen, NK.SP.rep3, by = "ID_REF")

NK.Immgen$FC <- log((rowMeans(NK.Immgen[, c("NK.IFNg.1", "NK.IFNg.2")]) / rowMeans(NK.Immgen[, c("NK.untx.1", "NK.untx.2", "NK.untx.3")])), 2)

affy_IDs <- biomaRt::select(ensembl_ms, keys = NK.Immgen$ID_REF, columns = c("affy_mogene_1_0_st_v1", "mgi_symbol", "hgnc_symbol"), keytype = "affy_mogene_1_0_st_v1")
NK.Immgen <- left_join(NK.Immgen, affy_IDs, join_by("ID_REF" == "affy_mogene_1_0_st_v1"))
NK.Immgen$hgnc_symbol <- toupper(NK.Immgen$mgi_symbol)

# BMDMs treated with IFNg
GSE199128_deseq2_IFNg_vs_UT <- read.csv2("GSE199128_deseq2_IFNg_vs_UT.csv")
GSE199128_deseq2_IFNg_vs_UT$hgnc_symbol <- toupper(GSE199128_deseq2_IFNg_vs_UT$X)
rownames(GSE199128_deseq2_IFNg_vs_UT) <- GSE199128_deseq2_IFNg_vs_UT$hgnc_symbol

# MDMs treated with LPS (and other TLR agonists)
MDM_RNAseq <- readxl::read_excel("GSE147311_RNA-seq_Mono_and_MDMs_time_course.xlsx")
MDM_RNAseq[, c(2:18)] <- MDM_RNAseq[, c(2:18)] + 0.1
MDM_RNAseq$Mono.TLR4.18hr.FC <- log2(MDM_RNAseq$`RNA-seq_D6_mono_TLR4_18hr` / MDM_RNAseq$`RNA-seq_D6_mono_unstimulated`)
MDM_RNAseq$MDM.TLR4.18hr.FC <- log2(MDM_RNAseq$`RNA-seq_D6_MDMs_TLR4_18hr` / MDM_RNAseq$`RNA-seq_D6_MDMs_unstimulated`)
row.names(MDM_RNAseq) <- MDM_RNAseq$`Feature ID`

############# Human Protein Atlas and Immgen Human #####################
GSE227743_Normalized_Gene_count_table <- read.csv("GSE227743_Normalized_Gene_count_table.csv")
Immgen_human_PBMCs <- GSE227743_Normalized_Gene_count_table
rownames(Immgen_human_PBMCs) <- Immgen_human_PBMCs$gene_symbol

Immgen_human_PBMCs$NK.imm <- rowMeans(Immgen_human_PBMCs[, c(10:11)])
Immgen_human_PBMCs$NK.mat <- rowMeans(Immgen_human_PBMCs[, c(12:13)])
Immgen_human_PBMCs$T.4Nve <- rowMeans(Immgen_human_PBMCs[, c(16:17)])
Immgen_human_PBMCs$T.4EffMem <- rowMeans(Immgen_human_PBMCs[, c(18:19)])
Immgen_human_PBMCs$Treg <- rowMeans(Immgen_human_PBMCs[, c(32:33)])
Immgen_human_PBMCs$Mo.16 <- rowMeans(Immgen_human_PBMCs[, c(34:35)])
Immgen_human_PBMCs$Mo.14 <- rowMeans(Immgen_human_PBMCs[, c(36:37)])

############# Mouse Immgen ########################

Immgen <- read.csv("GSE109125_Normalized_Gene_count_table.csv")
Immgen$Gene.upper <- toupper(Immgen$gene_symbol)
# remap symbols are in uppercase, even for mouse genes
rownames(Immgen) <- Immgen$Gene.upper

Immgen$MF.Alv.Lu <- rowMeans(Immgen[, c(69:70)])
Immgen$NK.imm <- rowMeans(Immgen[, c(81:82)])
Immgen$NK.mat <- rowMeans(Immgen[, c(89:90)])
Immgen$T.4.Nve.Sp <- rowMeans(Immgen[, c(115:116)])
Immgen$T.4.Th <- rowMeans(Immgen[, c(119:120)])
Immgen$MF.RP.Sp <- rowMeans(Immgen[, c(163:164)])
Immgen$Ly6C.neg <- rowMeans(Immgen[, c(169:170)])
Immgen$Ly6C.pos <- rowMeans(Immgen[, c(170:171)])
Immgen$Treg.4.25hi.Sp <- rowMeans(Immgen[, c(155:156)])

  
# ################## from bedtools intersect, pull out TFs #####################
# # from nonredundant peak set
# # created via bedtools intersect with ReMap dataset and tested human luciferase regions
TF_intersect <- read.csv("ReMap_intersect_nonredundant_hg38_luciferaseregions.txt",
                         sep = "\t", header = FALSE)
get_name_tissue <- data.frame(do.call('rbind', strsplit(as.character(TF_intersect$V11), ":", fixed = TRUE)))
TF_intersect <- cbind(TF_intersect, get_name_tissue$X1)
TF_intersect <- cbind(TF_intersect, get_name_tissue$X2)
TF_intersect <- TF_intersect[, c(1, 2, 3, 4, 8, 9, 10, 17, 18)]
colnames(TF_intersect) <- c("RE_chr", "RE_start", "RE_stop", "RE_name", "TF_chr",
                            "TF_start", "TF_stop", "TF_name", "TF_tissue")

# subset TF_intersect by LAIR1 and LAIR2 REs, since need to select for TF expression in different tissue types
LAIR1_REs <- c("RE50_53",   "RE54_55",   "RE56_61",   "RE62_63",   "RE64_69",   
               "RE70_71",   "RE72_73",   "RE74_75",   "RE74_83",  
               "RE84_87",   "RE86_87",   "RE88_91")
LAIR2_REs <- c("RE92_93",   "RE94_95",   "RE94_97",   "RE96_97",   "RE98_101",  "RE102_103",
               "RE104_105", "RE106_107", "RE106_113", "RE110_111", "RE114_121", "RE116_121",
               "RE118_121", "RE122_125", "RE130_133", "RE134_137", "RE140_143", "RE142_143")

# high vs low REs
LAIR1_high_REs <- c("RE62_63", "RE72_73", "RE74_75",  "RE74_83", "RE84_87", "RE86_87") 
LAIR2_high_REs <- c("RE94_95",   "RE94_97",   "RE96_97",  "RE102_103", "RE104_105", "RE134_137", "RE142_143")
LAIR1_low_REs <- LAIR1_REs[!(LAIR1_REs %in% LAIR1_high_REs)]
LAIR2_low_REs <- LAIR2_REs[!(LAIR2_REs %in% LAIR2_high_REs)]
TF_intersect_LAIR1 <- TF_intersect[TF_intersect$RE_name %in% LAIR1_REs,]
TF_intersect_LAIR2 <- TF_intersect[TF_intersect$RE_name %in% LAIR2_REs,]

############# create TF counts for LAIR1 and LAIR2
## LAIR1
TF_count_LAIR1 <- make_count_mat(TF_intersect_LAIR1)
TF_count_LAIR1 <- replace_nonspecific_ChIP_TFs(TF_count_LAIR1, GSE227743_Normalized_Gene_count_table$gene_symbol)
TF_count_LAIR1 <- find_enrichment(TF_count_LAIR1, LAIR1_low_REs, LAIR1_high_REs)
LAIR1_correlation <- correlation(TF_count_LAIR1$hgnc_symbol, "LAIR1", Immgen_human_PBMCs, c(10:13, 16:19, 32:37))
TF_count_LAIR1 <- left_join(TF_count_LAIR1, LAIR1_correlation, by = "hgnc_symbol")
TF_count_LAIR1 <- left_join(TF_count_LAIR1, MDM_RNAseq, join_by("hgnc_symbol" == "Feature ID"))

# LAIR1: correlated with LAIR1, expressed in MDM and monocytes, increases with LPS stimulation in MDM and monocytes
TF_LAIR1_prioritized <- dplyr::filter(TF_count_LAIR1, (correlation > 0.01)
                                        & (MDM.TLR4.18hr.FC > 0)
                                        & (Mono.TLR4.18hr.FC > 0)
                                        & (`RNA-seq_D6_MDMs_unstimulated` > 10)
                                        & (`RNA-seq_D6_mono_unstimulated` > 10)
                                      )
# remove TMEDx, BRDx, SMCx, MEDx, transcriptional coactivators, members of cohesin complex, etc
TF_LAIR1_prioritized <- as.data.frame(TF_LAIR1_prioritized[!grepl('TMED|BRD|SMC|MED|CDK|CTCF|P300|HDAC|CREB|BHLH', TF_LAIR1_prioritized$hgnc_symbol),])
TF_count_LAIR1_subset <- head(TF_LAIR1_prioritized[order(TF_LAIR1_prioritized$Mono.TLR4.18hr.FC, decreasing = TRUE),], n = 5)$hgnc_symbol

################# LAIR1 plot ########################################

LAIR1_formal_REs <- c("h1",   "h2",   "h3",   "h4",   "h5",  "h6",
                      "h7", "h8", "h9", "h10", "h11", "h12")
cell_types <- c("NK.mat", "T.4EffMem", "Treg", "Mo.14")

expr_data <- make_dotplot(TF_count_LAIR1, TF_count_LAIR1_subset, LAIR1_formal_REs, cell_types, Immgen_human_PBMCs)

pal <- colorRampPalette(brewer.pal(n = 9, name = "BuPu"))

ggplot(expr_data, aes(x=RE, y=cell, color = expr_color)) +
  geom_point(aes(size = sum)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(range = c(0, 20)) +
  scale_color_manual(values = pal(19))
ggsave("hLAIR1_expression.pdf")

############################## LAIR2 ##############################################
# set up count matrix
TF_count_LAIR2 <- make_count_mat(TF_intersect_LAIR2)
TF_count_LAIR2 <- replace_nonspecific_ChIP_TFs(TF_count_LAIR2, GSE227743_Normalized_Gene_count_table$gene_symbol)
TF_count_LAIR2 <- find_enrichment(TF_count_LAIR2, LAIR2_low_REs, LAIR2_high_REs)
LAIR2_correlation <- correlation(TF_count_LAIR2$hgnc_symbol, "LAIR2", Immgen_human_PBMCs, c(10:13, 16:19, 32:37))

# define gene sets
# LAIR2: correlated with LAIR2, expressed in T cells, increases after IFNg stimulation in NK cells
LAIR2_correlated <- LAIR2_correlation[(LAIR2_correlation$correlation > 0.01),]$hgnc_symbol
T_cell_expressed <- Immgen_human_PBMCs[Immgen_human_PBMCs$T.4EffMem.CD3.4.RA.62L..Bl.1 > 10 & Immgen_human_PBMCs$T.4EffMem.CD3.4.RA.62L..Bl.2 > 10,]$gene_symbol
TF_count_LAIR2 <- left_join(TF_count_LAIR2, NK.Immgen, by = "hgnc_symbol")

TF_LAIR2_prioritized <- dplyr::filter(TF_count_LAIR2, (hgnc_symbol %in% LAIR2_correlated) & 
                                            (hgnc_symbol %in% T_cell_expressed))

TF_LAIR2_prioritized <- as.data.frame(TF_LAIR2_prioritized[!grepl('TMED|BRD|SMC|MED|CDK|CTCF|P300|HDAC|RAD|CREB|BHLH', TF_LAIR2_prioritized$hgnc_symbol),])
TF_count_LAIR2_subset <- head(TF_LAIR2_prioritized[order(TF_LAIR2_prioritized$FC, decreasing = TRUE),], n = 5)$hgnc_symbol

########### LAIR2 plot #################################
LAIR2_formal_REs <- c("h13",   "h14",   "h15",   "h16",   "h17",  "h18",
                      "h19", "h20", "h21", "h22", "h23", "h24",
                      "h25", "h26", "h27", "h28", "h29", "h30")
cell_types <- c("NK.mat", "T.4EffMem", "Treg", "Mo.14")

expr_data <- make_dotplot(TF_count_LAIR2, TF_count_LAIR2_subset, LAIR2_formal_REs, cell_types, Immgen_human_PBMCs)

pal <- colorRampPalette(brewer.pal(n = 9, name = "BuPu"))

ggplot(expr_data, aes(x=RE, y=cell, color = expr_color)) +
  geom_point(aes(size = sum)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(range = c(0, 20)) +
  scale_color_manual(values = pal(19))
ggsave("hLAIR2_expression.pdf")

################### mouse Lair1 #######################
TF_intersect <- read.csv("ReMap_intersect_nonredundant_mm10_luciferaseregions.txt", 
                         sep = "\t", header = FALSE)
get_name_tissue <- data.frame(do.call('rbind', strsplit(as.character(TF_intersect$V13), ":", fixed = TRUE)))
TF_intersect <- cbind(TF_intersect, get_name_tissue$X1)
TF_intersect <- cbind(TF_intersect, get_name_tissue$X2)
TF_intersect <- TF_intersect[, c(1, 2, 3, 4, 10, 11, 12, 19, 20)]
colnames(TF_intersect) <- c("RE_chr", "RE_start", "RE_stop", "RE_name", "TF_chr", 
                            "TF_start", "TF_stop", "TF_name", "TF_tissue")

# define high and low REs
mLAIR1_high_REs <- c("RE_5-6", "RE_9-10_?", "RE_11-12_?")
mLAIR1_low_REs <- c("RE_1-2", "RE_3-4", "RE_7-8", "RE_13-14_?")

# set up count matrix
TF_count_mLAIR1 <- make_count_mat(TF_intersect)
TF_count_mLAIR1 <- replace_nonspecific_ChIP_TFs(TF_count_mLAIR1, Immgen$Gene.upper)
TF_count_mLAIR1 <- find_enrichment(TF_count_mLAIR1, mLAIR1_low_REs, mLAIR1_high_REs)

# define gene sets
# AVMs, naive NK cells, mature NK cells, T.4.Nve.Sp, T.4.Th, MF.RP.Sp, Ly6C low and high monocytes
Immgen_cols <- c(69:70, 81:82, 89:90, 115:116, 119:120, 163:164, 169:172)
mLAIR1_correlation <- correlation(TF_count_mLAIR1$hgnc_symbol, "LAIR1", Immgen, Immgen_cols)
mLAIR1_correlated <- mLAIR1_correlation[(mLAIR1_correlation$correlation > 0.01),]$hgnc_symbol
TF_count_mLAIR1 <- left_join(TF_count_mLAIR1, GSE199128_deseq2_IFNg_vs_UT, by = "hgnc_symbol")
MF_expressed <- Immgen[Immgen$MF.RP.Sp.1 > 10 & Immgen$MF.RP.Sp.2 > 10, ]$Gene.upper
# subset count matrix
TF_mLAIR1_prioritized <- TF_count_mLAIR1[(TF_count_mLAIR1$hgnc_symbol %in% mLAIR1_correlated) &
                                            (TF_count_mLAIR1$hgnc_symbol %in% MF_expressed)
                                          , ]
TF_LAIR2_prioritized <- as.data.frame(TF_mLAIR1_prioritized[!grepl('TMED|BRD|SMC|MED|CDK|CTCF|P300|HDAC|RAD|CREB|BHLH', TF_mLAIR1_prioritized$hgnc_symbol),])
TF_count_mLAIR1_subset <- head(TF_mLAIR1_prioritized[order(TF_mLAIR1_prioritized$log2FoldChange, decreasing = TRUE),], n = 5)$hgnc_symbol

###############################

Lair1_formal_REs <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7")

cell_types <- c("Treg.4.25hi.Sp", "NK.mat", "MF.RP.Sp", "Ly6C.pos")
rownames(TF_count_mLAIR1) <- TF_count_mLAIR1$hgnc_symbol
expr_data <- make_dotplot(TF_count_mLAIR1, TF_count_mLAIR1_subset, Lair1_formal_REs, cell_types, Immgen)

pal <- colorRampPalette(brewer.pal(n = 9, name = "BuPu"))

ggplot(expr_data, aes(x=RE, y=cell, color = expr_color)) +
  geom_point(aes(size = sum)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_size(range = c(0, 20)) +
  scale_color_manual(values = pal(17))
ggsave("mLair1_expression.pdf")

#################### generate info for UCSC tracks ################################
# create tracks for UCSC browser from non-redundant peak set 
# only open needed datasets

ReMap_all <- read.csv("remap2022_nr_macs2_hg38_v1_0.chr19.bed", 
                      sep = "\t", header = FALSE)
get_name_tissue <- data.frame(do.call('rbind', strsplit(as.character(ReMap_all$V4), ":", fixed = TRUE)))
colnames(ReMap_all) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
ReMap_all$name <- get_name_tissue$X1

for (x in c(TF_count_LAIR1_subset, TF_count_LAIR2_subset)) {
   write.table(ReMap_all[ReMap_all$name == x, ], file = paste0("Tracks_human/", x, ".bed"),
               row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

# mouse
ReMap_all <- read.csv("remap2022_nr_macs2_mm10_v1_0.chr7.bed", 
                      sep = "\t", header = FALSE)
get_name_tissue <- data.frame(do.call('rbind', strsplit(as.character(ReMap_all$V4), ":", fixed = TRUE)))
colnames(ReMap_all) <- c("chr", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")
ReMap_all$name <- get_name_tissue$X1

for (x in TF_count_mLAIR1_subset) {
  write.table(ReMap_all[ReMap_all$name == x, ], file = paste0("Tracks_mouse/", x, ".bed"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}




