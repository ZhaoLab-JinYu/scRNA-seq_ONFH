#January 5
#Single-cell RNA-seq workflow

##################### Load librarie ####################

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(SingleR)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)
library(NMF)
library(ComplexHeatmap)
library(CellChat)
rm(list = ls())
setwd("D:/research/")
getwd()
memory.limit(58000)

################################ 1.Reading Rawdata #############################

for (file in c('onfh1','onfh2','onfh3')){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)
}

merged_seurat <- merge(test_oa,test_op,add.cell.id = c('onfh1','onfh2','onfh3'))

merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

metadata <- merged_seurat@meta.data

metadata$cells <- rownames(metadata)

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "onfh1"))] <- 'onfh1'
metadata$sample[which(str_detect(metadata$cells, "onfh2"))] <- 'onfh2'
metadata$sample[which(str_detect(metadata$cells, "onfh3"))] <- 'onfh3'

merged_seurat@meta.data <- metadata

################################### 2.QC ######################################

filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.25))

counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

################################ 3.Normalize ##########################

seurat_phase <- NormalizeData(filtered_seurat)

############################### 4.Cell cycle (optional) ####################

# Load cell cycle markers(Download from github,markers of human cell G2/M and S phage)
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/cell_cycle_scoring.md
load("data/cycle.rda")

seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

seurat_phase <- ScaleData(seurat_phase)

seurat_phase <- RunPCA(seurat_phase)

DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

############################ 5.SCTransForm ############################

options(future.globals.maxSize = 4000 * 1024^6)

split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- list(group1,group2)

split_seurat <- split_seurat[c("onfh1",'onfh2','onfh3')]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

############################ 6.Integrate ###########################

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

################################ 7.Reducing dimensions ########################

seurat_integrated <- RunPCA(seurat_integrated, 
                            verbose = FALSE)

DimHeatmap(seurat_integrated, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)

ElbowPlot(object = seurat_integrated, 
          ndims = 40)

pc <- 40

seurat_integrated <- RunTSNE(seurat_integrated, 
                             reduction = "pca",
                             dims = 1:pc)

seurat_integrated <- RunUMAP(seurat_integrated, 
                             reduction = "pca", 
                             dims = 1:pc)

DimPlot(seurat_integrated, 
        reduction = "tsne",
        group.by = "ident",
        label=T,
        pt.size = 0.5)

DimPlot(seurat_integrated, 
        reduction = "umap",
        label=T,
        pt.size = 0.3)

############################ 7.Clustering ############################


seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:pc)

#Determine the clusters for various resolutions(optional)
# seurat_integrated <- FindClusters(object = seurat_integrated,
#                                  resolution = c(seq(.2,1.6,.2)))
# clustree test
#clustree(seurat_integrated)

seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.5)

seurat_integrated <- RunUMAP(seurat_integrated, 
                             reduction = "pca", 
                             dims = 1:40)

seurat_integrated <- RunTSNE(seurat_integrated,
                             reduction = "pca",
                             dims = 1:40)

DimPlot(seurat_integrated,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

######################## 8.SingleR (optional) #########################

# Load reference data(downloaded from SingleR github)
hpca.se <- readRDS(file = "ref/hpca.se.rds")
bp.se <- readRDS(file = "ref/bp.se.rds")

singleR_counts <- seurat_integrated@assays$RNA@counts

anno.total.bycluster <- SingleR(test = singleR_counts,
                                ref = list(HP = hpca.se , BP = bp.se),
                                labels = list(hpca.se$label.main , bp.se$label.main),
                                method = "cluster",
                                cluster = seurat_integrated$seurat_clusters)

plotScoreHeatmap(anno.total.bycluster)

anno.total.bycluster$cluster <- rownames(anno.total.bycluster)
anno <- anno.total.bycluster %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids <- anno$labels
names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
seurat_integrated <- AddMetaData(object = seurat_integrated, 
                                 metadata = seurat_integrated@active.ident,
                                 col.name = "celltype")

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

DimPlot(seurat_integrated,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)

####################### 9.Findmarkers ###################

DefaultAssay(seurat_integrated) <- "RNA"

seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  

compare_markers <- FindMarkers(object = seurat_integrated,
                               ident.1 = 2)

# annotation.csv was downloaded from github
# https://raw.githubusercontent.com/hbctraining/scRNA-seq/master/data/annotation.csv
annotations <- read.csv("data/annotation.csv")
markers <- markers %>% 
  left_join(gene,y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

write.csv(markers,file = "results/onfh_markers.csv")

top10 <- markers %>% 
  mutate(avg_log2FC) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, 
        wt = avg_log2FC)

seurat_integrated <- ScaleData(seurat_integrated)

DoHeatmap(object = seurat_integrated,
          features = top10$gene,
          label = FALSE) +
  scale_fill_gradientn(colors = c("#1976d2","white", "#c62828"))

############################ 10.Monocle ##############################

library(monocle)

Mono_tj<- seurat_integrated

Mono_matrix<-as(as.matrix(Mono_tj@assays$RNA@counts), 'sparseMatrix')

feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))

rownames(feature_ann)<-rownames(Mono_matrix)

Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)

sample_ann<-Mono_tj@meta.data

Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)

Mono.cds<-newCellDataSet(Mono_matrix,phenoData = Mono_pd,featureData = Mono_fd,expressionFamily=negbinomial.size())

head(pData(Mono.cds))

head(fData(Mono.cds))

Mono.cds <- estimateSizeFactors(Mono.cds)

Mono.cds <- estimateDispersions(Mono.cds)

disp_table <- dispersionTable(Mono.cds)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

markers <- read.csv(file = "results/onfh_markers.csv")

ordering_genes <- markers$gene

Mono.cds <- setOrderingFilter(Mono.cds,ordering_genes)

Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 2,
  method = 'DDRTree')

Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))

plot_cell_trajectory(Mono.cds,color_by = 'groups') + facet_wrap(~groups,nrow =1 )
plot_cell_trajectory(Mono.cds,color_by = 'Pseudotime')
plot_cell_trajectory(Mono.cds,color_by = 'State')
plot_cell_trajectory(Mono.cds,color_by = 'seurat_clusters')
plot_cell_trajectory(Mono.cds,color_by = 'Pseudotime',markers = 'ADIPOQ')

BEAM_res <- BEAM(Mono.cds, branch_point = 1, cores = 4)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]

BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(Mono.cds[row.names(subset(BEAM_res100,
                                                      qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 5,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)


############################ 11.CellChat #########################

library(CellChat)
library(ggplot2)
library(ggalluvial)
options(stringsAsFactors = FALSE)

seurat_integrated <- readRDS(file = "results/onfh+hoa+fnf_ec+msc_res0.1.rds")
data.input  <- seurat_integrated@assays$RNA@data
identity <- data.frame(group =seurat_integrated@active.ident,
                       row.names = names(seurat_integrated@active.ident))

unique(identity$group)

cellchat <- createCellChat(object = data.input)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human 
# Secreted Signaling or ECM-Receptor or Cell-Cell Contact
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_signalingRole(cellchat,
                                      slot.name = "netP")
cellchat@netP$pathways

vertex.receiver <- seq(1,15)
pathname <- cellchat@netP$pathways
pathways.show <- 'CHEMERIN'

netVisual_aggregate(cellchat, 
                    signaling = 'VISFATIN',
                    vertex.receiver = vertex.receiver, 
                    vertex.size = groupSize)

netVisual_aggregate(cellchat, 
                    signaling = 'VISFATIN', 
                    layout = "circle", 
                    vertex.size = groupSize)

netAnalysis_contribution(cellchat,
                         signaling = pathways.show)

cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")

netVisual_signalingRole(cellchat, signaling = pathways.show)

rankNet(cellchat,
        mode = "comparison", 
        sources.use = c('A+EC'), 
        targets.use = c('ADIPO','OSTEO','CHONDRO','MSC'), 
        comparison = c('ONFH_III','ONFH_IV','FNF','HOA'),
        stacked = F, 
        do.stat = TRUE)

plotGeneExpression(cellchat, features = "AEBP1")

########################### 12.Spearman #######################

library('openxlsx')
library('corrplot')
library('PerformanceAnalytics') 
seurat_integrated <- readRDS('results/MSC_ONFH_PC40_RES1.02.rds')
seurat_integrated <- subset(seurat_integrated, idents = list('26'))
DefaultAssay(seurat_integrated) <- "RNA"
markers <- read.csv('results/onfh_markers.csv')
markers <- markers[(markers$cluster == 26),]
markers <- markers[markers$avg_logFC>=1.5,]
mycors.p <- c('ITGA5','ITGB1','INSR','CD44','GLG1','SPP1','SP7','COL1A1','COL2A1','COL1A2','RUNX2','BGLAP','IBSP','PPARG','PDGFRB','PDGFRA')
cors_gene <- data.frame(FetchData(object = seurat_integrated, vars = mycors.p))
cor <- cor(cors_gene,method = 'spearman')   
cor <- data.frame(cor)
cor <- t(cor)
corrplot(cor,
         method="color",
         type = "upper",
         addCoef.col="grey",
         diag = F)

############################  13.Go ###########################

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)

seurat_integrated <- readRDS('results/onfh+hoa_MSC.rds')
seurat_integrated <- subset(seurat_integrated,idents = '1')
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

deg <- FindMarkers(seurat_integrated, 
                   ident.1 = 'ONFH',
                   ident.2 = 'HOA',
                   group.by = 'group')

deg <- read.csv(file = 'results/onfh_markers.csv')
deg <- data.frame(gene = rownames(deg), deg)
deg1 <- deg
logFC_t <- 0.5
P.Value_t <- 0.05
k1 <- (deg1$p_val_adj < P.Value_t)&(deg1$avg_logFC < -logFC_t)
k2 <- (deg1$p_val_adj < P.Value_t)&(deg1$avg_logFC > logFC_t)
change <- ifelse(k1,"down",ifelse(k2,"up","stable"))
deg1$change <- change

s2e <- bitr(deg1$X, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db) 

deg1 <- inner_join(deg1,s2e,by=c("X"="SYMBOL"))

gene_up <- deg1[deg1$change == 'up','ENTREZID'] 
gene_down <- deg1[deg1$change == 'down','ENTREZID'] 
gene_diff <- c(gene_up,gene_down)
gene_all <- deg1[,'ENTREZID']

enrich<-enricher(gene=gene_up,)

cnetplot(enrich,
         vertex.label.cex=0.1,
         colorEdge = TRUE)

ego_CC <- enrichGO(gene = gene_up,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene_up,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene_up,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

##########################  14.Plots #######################

DefaultAssay(seurat_integrated) <- "RNA"

#DOTPLOT 
DotPlot(seurat_integrated,
        idents = '1',
        features = 'PPARG',
        cols = c('#ffd54f','#d32f2f'),
        dot.min = 0,
        group.by = 'groups')

# Bar plot 
ggplot(ggplot_data,aes(x = seurat_clusters, fill = Phase)) + 
  geom_bar(position = 'dodge') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label=..count..),stat = 'count')

# BOX PLOT 
genelist <- c('APOE','PPARG','LEPR','CEBPA','CEBPB','TFAP2A','SLC2A4','ADIPOQ','LEP','PNPLA2','LPL','PLIN1','FABP4')
ggplot_gene <- data.frame(FetchData(object = seurat_integrated, vars = genelist,slot = 'data'))
ggplot_gene <- data.frame(ggplot_gene,clusters = seurat_integrated@meta.data$seurat_clusters)
ggplot_gene <- data.frame(ggplot_gene,groups = seurat_integrated@meta.data$groups)
ggplot_gene <- data.frame(ggplot_gene,group = seurat_integrated@meta.data$group)

my_comparisons <- list( c("HOA", "ONFH_III"), c("HOA", "ONFH_IV"), c("ONFH_III", "ONFH_IV") )

ggboxplot(ggplot_gene,
          x = 'groups',
          y = genelist,
          combine = TRUE,) + 
stat_compare_means(comparisons = my_comparisons, label =  "p.signif")#+ # Add pairwise comparisons p-value

#Feature plot 
FeaturePlot(object = seurat_integrated, 
            features = 'COL6A2',
            reduction = "umap",
            sort.cell = TRUE,
            #split.by = 'group',
            blend = FALSE,
            max.cutoff = "q90",
            label = FALSE,
            repel = TRUE,
            cols = c("light grey","#d50000"),
            by.col = TRUE)

# Double gene Feature plot
FeaturePlot(object = seurat_integrated, 
            features = c("LEPR","ADIPOQ"),
            reduction = "umap",
            sort.cell = TRUE,
            blend = TRUE,
            min.cutoff = "q10",
            max.cutoff = "q90",
            label = TRUE,
            repel = TRUE,
            cols = c("#388e3c","#d32f2f"),
            by.col = TRUE)


#VLNPLOT
features <- c('BGLAP','RUNX2','IBSP','PPARG','ADIPOQ','APOE','NOTCH3','MCAM','RGS5','COL2A1','ACAN','SOX9','COL1A1','PRG4','COL3A1')
VlnPlot(object = seurat_integrated, 
        pt.size = 0,
        features = features)


#Volcano plot
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
compare_markers <- FindMarkers(object = seurat_integrated,
                               ident.1 = 'ONFH_III',
                               ident.2 = 'FNF',
                               group.by = 'groups')
EnhancedVolcano(compare_markers,
                lab = rownames(compare_markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(-1.5,2),
                FCcutoff = 0.5,
                ylim = c(0,60),
                colAlpha = 1)
