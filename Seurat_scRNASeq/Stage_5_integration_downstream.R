# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <in_dir> <out_prefix> <cluster res> <SampleSheet_file> <Output_dir>\n"

args = commandArgs(trailingOnly = TRUE)
if(length(args)!=5) {
  stop("\n    Wrong parameters. \n    ", usage)
}


# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RCurl))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(AnnotationHub))


#in_dir = "/scratch1/01775/saathe/DG_china_cohort/"
#SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
#out_prefix = "Pat_0013"
#res = "0.4"
#res_string <- paste0("integrated_snn_res.",res)

in_dir = args[1]
#SamplesSheet_file = args[2]
#SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = args[2]
res = args[3]
SamplesSheet_file = args[4]
out_dir = args[5]

res_string <- paste0("integrated_snn_res.",res)

in_dir
out_prefix
SamplesSheet_file
res

setwd(in_dir)
#in_dir <- "D:/Data/HY_Yapeng/CUS47_10671/"
#out_dir <- paste0(in_dir,"/Seurat_output/data/results/")
dir.create(out_dir, recursive = TRUE)
dir.create(paste0(out_dir,"QC/"))

SamplesSheet <- read.csv(SamplesSheet_file,header=TRUE)

files <- list.files(path = paste0(in_dir,"/Seurat_input_data/"), pattern = "*.h5", full.names = T)

paste0("Input Directory: ",in_dir)
paste0("Out prefix: ",out_prefix)
paste0("SamplesSheet file: ",SamplesSheet_file)
paste0("Output Directory: ",out_dir)
paste0("Output QC Directory: ",out_dir,"QC/")

paste0("SamplesSheet input: \n")
SamplesSheet
paste0("Input Hdf5r files: \n")
files


# Save integrated seurat object
seurat_integrated <- readRDS(paste0(out_dir,"/",out_prefix,"_integrated_seurat.rds"))

seurat_integrated <- FindClusters(seurat_integrated, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4))

# Assign identity of clusters
Idents(object = seurat_integrated) <- res_string

#seurat_integrated$Group <- NA
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "M11914_5d"))] <- "M11914_5d"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "T9236_2wk"))] <- "T9236_2wk"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "T9233_4wk"))] <- "T9233_4wk"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "AL_Control"))] <- "AL_Control"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "Triple_TDT_mouse"))] <- "Triple_TDT_mouse"

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
        dplyr::count(ident, sample) %>%
        tidyr::spread(ident, n)

# View table
#View(n_cells)

# create color palette:
library(RColorBrewer)
#coul <- brewer.pal(3, "Pastel2") 
coul <- c("#F8766D", "#00BFC4")
# Transform this data in %
data_percentage <- apply(data.matrix(n_cells[,-1]), 2, function(x){x*100/sum(x,na.rm=T)})
rownames(data_percentage) <- n_cells[,1];  colnames(data_percentage) <- colnames(n_cells[,-1])

# Make a stacked barplot--> it will be in %!
png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_BarPlot_percent_bycondition.png"), width = 12, height = 12, units = 'in', res = 300)
par(xpd=T, mar=par()$mar+c(0,0,0,6))
barplot(data_percentage, col=coul , border="white", xlab="condition", ylab= "percentage of cells"); legend("topright", legend=rownames(data_percentage), inset=c(-0.02,0), col=coul, pch = 19, bty='n')
dev.off()


# UMAP of cells in each cluster by sample
png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_integrated_snn_res.",res,"_bysample.png"), width = 18, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
dev.off()
		
# Explore whether clusters segregate by cell cycle phase
png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_integrated_snn_res.",res,"_bycellcycphase.png"), width = 14, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
dev.off()
		
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_integrated_snn_res.",res,"_FeaturePlot.png"), width = 12, height = 12, units = 'in', res = 300)                           
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
					 
# Extract the UMAP coordinates for the first 10 cells
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
# Plotting a UMAP plot for each of the PCs
png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_umap_pcs.png"), width = 14, height = 12, units = 'in', res = 300)                           
map(paste0("PC_", 1:16), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)
dev.off()
		
# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

#Exploring known cell type markers
png(paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_dimplot-markers.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()
		
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

#CD14+ monocyte markers
#FeaturePlot(seurat_integrated, 
#            reduction = "umap", 
#            features = c("TP53", "LYZ"), 
#            sort.cell = TRUE,
#            min.cutoff = 'q10', 
#            label = TRUE)
			
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_integrated, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)           
						  
DefaultAssay(seurat_integrated) <- "RNA"
tmp <- paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_samples_dge_markers.txt")
write.table(markers,tmp,sep="\t",row.names=TRUE)

##################################################################

sel.clust = res_string
markers[order(markers$p_val_adj),] -> cluster.markers
head(cluster.markers)

cluster.markers %>%
  group_by(cluster) %>%
    slice(1:25)  -> topgene.per.cluster
	#%>%
    #  pull(gene)

#We can now select the top 25 up regulated genes for plotting.



#par(1, 5, mar = c(4, 6, 3, 1))
tmp <- paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_samples_dge_barplot.pdf")
pdf(tmp, width = 16, height = 14)
par(1, 5, mar = c(4, 6, 3, 1))
par(mfrow=c(3,6),mgp=c(2, 0.2, 0))
#plot_grid(ncol = 3, for (i in unique(top25$cluster)) {
#    barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T, 
#        las = 1, main = paste0(i, " vs. rest"), xlim = c(0,as.numeric(max(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F)))), border = "white", yaxs = "i")
#    abline(v = c(0, 0.25), lty = c(1, 2))
#})
#top25$avg_log2FC[is.infinite(top25$avg_log2FC)] <- 600  
#lapply(unique(top25$cluster), function(i) {barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T,las = 1, xlab = "log2(Fold Change)", main = paste0(i, " vs. rest"), xlim = c(0,as.numeric(max(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F)))), border = "white", yaxs = "i"); abline(v = c(0, 0.25), lty = c(1, 2))})
lapply(unique(topgene.per.cluster$cluster), function(i) {barplot(sort(setNames(topgene.per.cluster$avg_log2FC, topgene.per.cluster$gene)[topgene.per.cluster$cluster == i], F), horiz = T,las = 1, xlab = "log2(Fold Change)", main = paste0(i, " vs. rest"), xlim = c(0,as.numeric(max(sort(setNames(topgene.per.cluster$avg_log2FC, topgene.per.cluster$gene)[topgene.per.cluster$cluster == i], F)))), border = "white", yaxs = "i"); abline(v = c(0, 0.25), lty = c(1, 2))})
dev.off()

#We can visualize them as a heatmap. Here we are selecting the top 5.



cluster.markers %>%
  group_by(cluster) %>%
    slice(1:5) -> top5
	
cluster.markers %>%
  group_by(cluster) %>%
    slice(1:3) -> top3
	
cluster.markers %>%
  group_by(cluster) %>%
    slice(1:2) -> top2
	
#create a scale.data slot for the selected genes
seurat_integrated <- ScaleData(seurat_integrated, features = as.character(unique(top5$gene)), assay = "RNA")
tmp <- paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_samples_dge_heatmap.png")
png(tmp, width = 12, height = 12, units = "in", res=300)
DoHeatmap(seurat_integrated, features = as.character(unique(top3$gene)), group.by = sel.clust, 
    assay = "RNA")
dev.off()

#Another way is by representing the overal group expression and detection rates in a dot-plot.

tmp <- paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_samples_dge_dotplot.pdf")
pdf(tmp, width = 14, height = 16)
DotPlot(seurat_integrated, features = as.character(unique(top2$gene)), group.by = sel.clust, 
    assay = "RNA") + coord_flip()
dev.off()

#We can also plot a violin plot for each gene.


# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1

tmp <- paste0(out_dir,"/Integrated_data_figures/",out_prefix,"_samples_dge_vlnplot.pdf")
pdf(tmp, width = 26, height = 26)
p = VlnPlot(object = seurat_integrated, pt.size = 0, same.y.lims = FALSE, features = as.character(unique(top2$gene)), ncol = 6, group.by = sel.clust, assay = "RNA", combine = FALSE)
p1 = lapply(X = p, FUN = function(x) x + theme(axis.text.x = element_text(size = 16,  family="Helvetica"), axis.text.y = element_text(size = 20,  family="Helvetica"), axis.title = element_text(size = 0,  family="Helvetica"), legend.text = element_text(size = 0,  family="Helvetica")) + theme(legend.position = "none"))
p2 = CombinePlots(plots = p1)
p2
dev.off()

##################################################################

# Create function to get conserved markers for any given cluster

Per_cluster_conserved <- function(cluster){
        FindConservedMarkers(seurat_integrated,
                             ident.1 = cluster,
                             grouping.var = "sample",
                             only.pos = TRUE) %>%
                rownames_to_column(var = "gene") %>%
                cbind(cluster_id = cluster, .)
}

#map_dfr(inputs_to_function, name_of_function)

# Iterate function across desired clusters
#cons <- map_dfr(levels(seurat_integrated@active.ident), Per_cluster_conserved)
#
#Conserved_genes_per_cluster <- table(cons$cluster_id)
#tmp <- paste0(out_dir,"Integrated_data_figures/",out_prefix,"_conserved_markers.txt")
#write.table(cons,tmp,sep="\t",row.names=TRUE)


# Extract top 10 markers per cluster
#top25 <- cons %>% 
#  mutate(avg_fc = (KO_avg_log2FC + WT_avg_log2FC) /2) %>% 
#  group_by(cluster_id) %>% 
#  top_n(n = 25, 
#        wt = avg_fc)

# Visualize top 10 markers per cluster
#View(top25)  

#tmp <- paste0(out_dir,"Integrated_data_figures/",out_prefix,"_samples_conservedgenes_barplot.pdf")
#pdf(tmp, width = 14, height = 14)
#par(1, 5, mar = c(4, 6, 3, 1))
#par(mfrow=c(3,6))
#lapply(unique(top25$cluster_id), 
#function(i) {
#barplot(sort(setNames(((top25$KO_avg_log2FC + top25$WT_avg_log2FC) /2), top25$gene)[top25$cluster_id == i], F), horiz = T,las = 1, xlab = "log2(Fold Change)", main = paste0(i, " conserved in WT and KO"), xlim = c(0,as.numeric(max(sort(setNames(((top25$KO_avg_log2FC + top25$WT_avg_log2FC) /2), top25$gene)[top25$cluster_id == i], F)))), border = "white", yaxs = "i"); abline(v = c(0, 0.25), lty = c(1, 2))
#})
#dev.off()


# Save integrated seurat object
saveRDS(seurat_integrated, paste0(out_dir,"/",out_prefix,"_final_integrated_seurat.rds"))

# Plot interesting marker gene expression for cluster 20
#FeaturePlot(object = seurat_integrated, 
#                        features = c("TPSAB1", "TPSB2", "FCER1A", "GATA1", "GATA2"),
#                         sort.cell = TRUE,
#                         min.cutoff = 'q10', 
#                         label = TRUE,
#			 repel = TRUE)
			 
# Vln plot - cluster 20
#VlnPlot(object = seurat_integrated, 
#        features = c("TPSAB1", "TPSB2", "FCER1A", "GATA1", "GATA2"))

# Rename all identities
#seurat_integrated <- RenameIdents(object = seurat_integrated, 
#                               "0" = "Naive or memory CD4+ T cells",
#                               "1" = "CD14+ monocytes",
#                               "2" = "Naive or memory CD4+ T cells",
#                               "3" = "CD14+ monocytes",
#                               "4" = "CD4+ T cells",
#                               "5" = "CD8+ T cells",
#                               "6" = "B cells",
#                               "7" = "Stressed cells / Activated T cells",
#                               "8" = "NK cells",
#                               "9" = "FCGR3A+ monocytes",
#                               "10" = "CD4+ T cells",
#                               "11" = "B cells",
#                               "12" = "NK cells",
#                               "13" = "CD8+ T cells",
#                               "14" = "CD14+ monocytes",
#                               "15" = "Conventional dendritic cells",
#			       "16" = "Megakaryocytes",
#			       "17" = "B cells", 
#			       "18" = "CD4+ T cells", 
#			       "19" = "Plasmacytoid dendritic cells", 
#			       "20" = "Mast cells")


# Plot the UMAP
#DimPlot(object = seurat_integrated, 
#        reduction = "umap", 
#        label = TRUE,
#        label.size = 3,
#        repel = TRUE)

# Remove the stressed or dying cells
#seurat_subset_labeled <- subset(seurat_integrated,
#                               idents = "Stressed cells / Activated T cells", invert = TRUE)

# Re-visualize the clusters
#DimPlot(object = seurat_subset_labeled, 
#        reduction = "umap", 
#        label = TRUE,
#        label.size = 3,
#	repel = TRUE)


# Save final R object
#write_rds(seurat_integrated,
#          path = "results/seurat_labelled.rds")  
