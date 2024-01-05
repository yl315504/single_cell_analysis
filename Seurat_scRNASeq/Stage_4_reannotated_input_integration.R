
# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <in_dir> <out_prefix> <Samplesheet_metascape>\n"

args = commandArgs(trailingOnly = TRUE)
if(length(args)!=3) {
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

in_dir = args[1]
SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = args[2]
org = args[3]

setwd(in_dir)
#in_dir <- "D:/Data/HY_Yapeng/CUS47_10671/"
out_dir <- paste0(in_dir,"/Seurat_output/data/results/")


#org="mouse"

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

# The first sample
E15_5_creneg <- readRDS("/scratch/08270/yanliu74/KD/2021_09_28_10X32_11234/Seurat_output/data/results/E15_5_crenegfinal_reclustered_annotated.rds",ref=NULL)

# The second sample
E15_5_mutant <- readRDS("/scratch/08270/yanliu74/KD/2021_09_28_10X32_11234/Seurat_output/data/results/E15_5_mutantfinal_reclustered_annotated.rds",ref=NULL)

split_seurat <- list(E15_5_creneg,E15_5_mutant)

for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
    split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
    }

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000) 
											
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
								   
# Find best buddies - can take a while to run
k.filter <- min(200, min(sapply(split_seurat, ncol)))

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", k.filter = k.filter,
                                        anchor.features = integ_features)
										
memory.limit(56000)									
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
		

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

dir.create(paste0(out_dir,"Integrated_data_figures_after_reannot/"))
# Plot PCA
png(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_PCAPlot.png"), width = 12, height = 12, units = 'in', res = 300)    
PCAPlot(seurat_integrated, split.by = "sample")
dev.off()		
		
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:40, reduction = "pca")

# Plot UMAP
png(paste0(out_dir,"Integrated_data_figures_after_reannot/", out_prefix, "_DimPlot_bycondition.png"), width = 12, height = 12, units = 'in', res = 300)
DimPlot(seurat_integrated)
dev.off()

pdf(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_DimPlot_bysample.pdf"), width = 20, height = 16)
DimPlot(seurat_integrated, split.by = "sample")
dev.off()		
		
		
# Explore heatmap of PCs
png(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_DimHeatmap_PCA.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()
		   
# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
	  
# Plot the elbow plot
png(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_ElbowPlot.png"), width = 12, height = 12, units = 'in', res = 300)                           
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()
		  
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:40)
                                
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

suppressPackageStartupMessages(library(clustree))
tmp <- paste(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_integrated_cl_clustree.png",sep="")
png(tmp, width = 12, height = 8, units = "in", res=300)
clustree(seurat_integrated@meta.data, prefix = "integrated_snn_res.")
dev.off()
							   
# Explore resolutions
#seurat_integrated@meta.data %>% 
#        View()
		
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
png(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_integrated_snn_res.0.8.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()
		
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
png(paste0(out_dir,"Integrated_data_figures_after_reannot/",out_prefix,"_integrated_snn_res.0.4.png"), width = 12, height = 12, units = 'in', res = 300)                           
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()		
			
		
		
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
		
# Save integrated seurat object
saveRDS(seurat_integrated, paste0(out_dir,"/",out_prefix,"_integrated_reclustered_annotated.rds"))
