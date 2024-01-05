# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <in_dir> <out_prefix> <cluster_res> <organism>\n"

args = commandArgs(trailingOnly = TRUE)
if(length(args)!=4) {
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
SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = args[2]
res = args[3]
#res = "0.8"
res_string <- paste0("SCT_snn_res.",res)
#res_string <- paste0("integrated_snn_res.",res)
org = args[4]

in_dir
out_prefix
SamplesSheet_file
res

setwd(in_dir)
#in_dir <- "D:/Data/HY_Yapeng/CUS47_10671/"
out_dir <- paste0(in_dir,"/Seurat_output/data/results/")
dir.create(out_dir, recursive = TRUE)
dir.create(paste0(out_dir,"QC/"))

#SamplesSheet <- read.csv(SamplesSheet_file,header=TRUE)

#files <- list.files(path = paste0(in_dir,"/Seurat_input_data/"), pattern = "*.h5", full.names = T)

paste0("Input Directory: ",in_dir)
paste0("Out prefix: ",out_prefix)
#paste0("SamplesSheet file: ",SamplesSheet_file)
paste0("Output Directory: ",out_dir)
paste0("Output QC Directory: ",out_dir,"QC/")

#paste0("SamplesSheet input: \n")
#SamplesSheet
#paste0("Input Hdf5r files: \n")
#files


# Read split individual seurat object
seurat_indiv <- readRDS(paste0(out_dir,"/",out_prefix,"_seurat.rds"))
# Single-cell RNA-seq analysis - clustering analysis

if (org == "mouse") {
	m.cc.genes <- readRDS("/home1/08270/yanliu74/Seurat_scRNASeq/mouse_cell_cycle_genes.rds")
	s_genes <- m.cc.genes$s.genes
	g2m_genes <- m.cc.genes$g2m.genes
} else if (org == "human") {
	cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
	cell_cycle_genes <- read.csv(text = cc_file)

	# Connect to AnnotationHub
	ah <- AnnotationHub()

	# Access the Ensembl database for organism
	ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

	# Acquire the latest annotation files
	id <- ahDb %>%
        	mcols() %>%
        	rownames() %>%
        	tail(n = 1)

	# Download the appropriate Ensembldb database
	edb <- ah[[id]]

	# Extract gene-level information from database
	annotations <- genes(edb, 
                     return.type = "data.frame")

	# Select annotations of interest
	annotations <- annotations %>%
        		dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

	# Get gene names for Ensembl IDs for each gene
	cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

	# Acquire the S phase genes
	s_genes <- cell_cycle_markers %>%
        	dplyr::filter(phase == "S") %>%
        	pull("gene_name")
        
	# Acquire the G2M phase genes        
	g2m_genes <- cell_cycle_markers %>%
        	dplyr::filter(phase == "G2/M") %>%
        	pull("gene_name")
} else {
print("Organism other than human or mouse")
}

seurat_indiv <- NormalizeData(seurat_indiv, verbose = TRUE)
seurat_indiv <- CellCycleScoring(seurat_indiv, g2m.features=g2m_genes, s.features=s_genes)
seurat_indiv <- SCTransform(seurat_indiv, vars.to.regress = c("mitoRatio"))

# Run PCA
seurat_indiv <- RunPCA(object = seurat_indiv)

dir.create(paste0(out_dir,"Individual_data_figures/"))
dir.create(paste0(out_dir,"Individual_data_figures/",out_prefix))

# Plot PCA
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_PCAPlot.png"), width = 12, height = 12, units = 'in', res = 300)    
PCAPlot(seurat_indiv)
dev.off()		
		
# Run UMAP
seurat_indiv <- RunUMAP(seurat_indiv, 
                             dims = 1:40,
			     reduction = "pca")

# Plot UMAP  
#png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_DimPlot_bycondition.png"), width = 12, height = 12, units = 'in', res = 300)                           
#DimPlot(seurat_indiv)
#dev.off()
		
		
# Explore heatmap of PCs
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_DimHeatmap_PCA.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimHeatmap(seurat_indiv, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()
		   
# Printing out the most variable genes driving PCs
print(x = seurat_indiv[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
	  
# Plot the elbow plot
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_ElbowPlot.png"), width = 12, height = 12, units = 'in', res = 300)                           
ElbowPlot(object = seurat_indiv, 
          ndims = 40)
dev.off()
		  
# Determine the K-nearest neighbor graph
seurat_indiv <- FindNeighbors(object = seurat_indiv, 
                                dims = 1:40)
                                
# Determine the clusters for various resolutions                                
seurat_indiv <- FindClusters(object = seurat_indiv,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

suppressPackageStartupMessages(library(clustree))
tmp <- paste(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_integrated_cl_clustree.png",sep="")
png(tmp, width = 12, height = 8, units = "in", res=300)
clustree(seurat_indiv@meta.data, prefix = "SCT_snn_res.")
dev.off()
							   
# Explore resolutions
#seurat_indiv@meta.data %>% 
#        View()
		
# Assign identity of clusters
Idents(object = seurat_indiv) <- "SCT_snn_res.0.8"

# Plot the UMAP
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_SCT_snn_res.0.8.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_indiv,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()
		
# Assign identity of clusters
Idents(object = seurat_indiv) <- "SCT_snn_res.0.4"

# Plot the UMAP
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_SCT_snn_res.0.4.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_indiv,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()		
			
# Assign identity of clusters
Idents(object = seurat_indiv) <- res_string

# UMAP of cells in each cluster by sample
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_snn_res.",res,".png"), width = 18, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_indiv, 
        label = TRUE)  + NoLegend()
dev.off()
		
# Explore whether clusters segregate by cell cycle phase
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_snn_res.",res,"_bycellcycphase.png"), width = 14, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_indiv,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
dev.off()
		
# Determine metrics to plot present in seurat_indiv@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")


png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_snn_res.",res,"_FeaturePlot.png"), width = 12, height = 12, units = 'in', res = 300)                           
FeaturePlot(seurat_indiv, 
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
pc_data <- FetchData(seurat_indiv, 
                     vars = columns)
					 
# Extract the UMAP coordinates for the first 10 cells
seurat_indiv@reductions$umap@cell.embeddings[1:10, 1:2]

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_indiv, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
# Plotting a UMAP plot for each of the PCs
png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_umap_pcs.png"), width = 14, height = 12, units = 'in', res = 300)                           
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
print(seurat_indiv[["pca"]], dims = 1:5, nfeatures = 5)

#Exploring known cell type markers
#png(paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_dimplot-markers.png"), width = 12, height = 12, units = 'in', res = 300)                           
#DimPlot(object = seurat_indiv, 
#        reduction = "umap", 
#        label = TRUE) + NoLegend()
		
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_indiv) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_indiv <- NormalizeData(seurat_indiv, verbose = FALSE)

#CD14+ monocyte markers
#FeaturePlot(seurat_indiv, 
#            reduction = "umap", 
#            features = c("TP53", "LYZ"), 
#            sort.cell = TRUE,
#            min.cutoff = 'q10', 
#            label = TRUE)
			
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_indiv, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)           
						  
DefaultAssay(seurat_indiv) <- "RNA"
tmp <- paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_samples_dge_markers.txt")
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
tmp <- paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_samples_dge_barplot.pdf")
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
seurat_indiv <- ScaleData(seurat_indiv, features = as.character(unique(top5$gene)), assay = "RNA")
tmp <- paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_samples_dge_heatmap.png")
png(tmp, width = 12, height = 12, units = "in", res=300)
DoHeatmap(seurat_indiv, features = as.character(unique(top3$gene)), group.by = sel.clust, 
    assay = "RNA")
dev.off()

#Another way is by representing the overal group expression and detection rates in a dot-plot.

tmp <- paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_samples_dge_dotplot.pdf")
pdf(tmp, width = 14, height = 16)
DotPlot(seurat_indiv, features = as.character(unique(top2$gene)), group.by = sel.clust, 
    assay = "RNA") + coord_flip()
dev.off()

#We can also plot a violin plot for each gene.


# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1

tmp <- paste0(out_dir,"Individual_data_figures/",out_prefix,"/",out_prefix,"_dge_vlnplot.pdf")
pdf(tmp, width = 26, height = 26)
p = VlnPlot(object = seurat_indiv, pt.size = 0, same.y.lims = FALSE, features = as.character(unique(top2$gene)), ncol = 6, group.by = sel.clust, assay = "RNA", combine = FALSE)
p1 = lapply(X = p, FUN = function(x) x + theme(axis.text.x = element_text(size = 16,  family="Helvetica"), axis.text.y = element_text(size = 20,  family="Helvetica"), axis.title = element_text(size = 0,  family="Helvetica"), legend.text = element_text(size = 0,  family="Helvetica")) + theme(legend.position = "none"))
p2 = CombinePlots(plots = p1)
p2
dev.off()

saveRDS(seurat_indiv,paste0(out_dir,"/",out_prefix,"final_seurat.rds"))



