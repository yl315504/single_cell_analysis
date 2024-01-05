# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <COUNTS_FILE> <SAMPLE_FILE> <OUT_DIR> <EDGER_BCV> <CONTRAST_FILE> <TPM_file> <AnnotationFile|none>\n"

args = commandArgs(trailingOnly = TRUE)
if(length(args)!=2) {
  cat(usage1)
  stop("\n Wrong parameters. See usage above.\n")
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


#in_dir = "/scratch/01775/saathe/Ralf_Kittler/Bachoo_Mouse_samples/"
#SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
#out_prefix = "Mouse_5sam"

in_dir = args[1]
SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = args[2]


setwd(in_dir)
#in_dir <- "D:/Data/HY_Yapeng/CUS47_10671/"
out_dir <- paste0(in_dir,"/Seurat_output/data/results/")
dir.create(out_dir, recursive = TRUE)
dir.create(paste0(out_dir,"QC/"))

SamplesSheet <- read.csv(SamplesSheet_file,header=TRUE)

#files <- list.files(path = paste0(in_dir,"/Seurat_input_data/"), pattern = "*.h5", full.names = T)
files <- list.dirs(path = paste0(in_dir,"/Seurat_input_data/"), full.names = TRUE, recursive = FALSE)

paste0("Input Directory: ",in_dir)
paste0("Out prefix: ",out_prefix)
paste0("SamplesSheet file: ",SamplesSheet_file)
paste0("Output Directory: ",out_dir)
paste0("Output QC Directory: ",out_dir,"QC/")

paste0("SamplesSheet input: \n")
SamplesSheet
paste0("Input Hdf5r files: \n")
files

# Read merged before filter seurat object
merged_seurat <- readRDS(file=paste0(out_dir,"/merged_filtered_seurat.rds"),ref=NULL)

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))

#Specific for Tyler Miller Audrey Chang kidney cells # Filter out low quality reads using selected thresholds - these will change with experiment 
#filtered_seurat <- subset(x = merged_seurat, 
#                         subset= (nUMI >= 200) & 
#                           (nGene >= 100) & 
#                           (log10GenesPerUMI > 0.80) & 
#                           (mitoRatio < 0.20))

						   
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Create .RData object to load at any time
saveRDS(filtered_seurat, file=paste0(out_dir,"seurat_filtered.rds"))


######################################################
# Create metadata dataframe
metadata <- filtered_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
#metadata <- metadata %>%
#        dplyr::rename(seq_folder = orig.ident,
#                      nUMI = nCount_RNA,
#                      nGene = nFeature_RNA)
					  
# Create sample column
metadata$sample <- NA
# Creat group column; This is for Chunyang only.
metadata$group <- NA

for(x in 1:length(files)){
metadata$sample[which(str_detect(metadata$cells, SamplesSheet$Type[x]))] <- SamplesSheet$Type[x]
# This line is for Chunyang only!!!!
metadata$group[which(str_detect(metadata$cells, SamplesSheet$Group[x]))] <- SamplesSheet$Group[x]
}

# Add metadata back to Seurat object
filtered_seurat@meta.data <- metadata

saveRDS(filtered_seurat, file=paste0(out_dir,"seurat_filtered.rds"))
# run garbage collect to free up memory
gc()

#paste0(in_dir,"data/")
#out_dir <- paste0(in_dir,"data/results/")
#dir.create(out_dir)
dir.create(paste0(out_dir,"QC/After_filter/"))


# Visualize the number of cell counts per sample
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_Cellnums_Sample.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
	ggtitle("Number of cell counts per sample") + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(size=22,  family="Helvetica", angle = 45, vjust = 1, hjust=1)) +
	theme(text = element_text(size=22, family="Helvetica", colour="black")) +
	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
dev.off()

# Visualize the number UMIs/transcripts per cell
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_UMIsvstranscripts_Sample.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
	ggtitle("Number UMIs/transcripts per cell") + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
	theme(text = element_text(size=22,  family="Helvetica", colour="black")) +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
dev.off()

# Visualize the distribution of genes detected per cell via histogram
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_hist_genespercell.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
	ggtitle("Histogram of genes per cell") + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
	theme(text = element_text(size=22, family="Helvetica", colour="black")) +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_boxp_genespercell.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
	ggtitle("Boxplot of genes per cell") + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
	theme(text = element_text(size=22,  family="Helvetica", colour="black")) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")
dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_genesvsUMIs.png"), width = 18, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
	ggtitle("Correlation between genes and UMIs") + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
	theme(text = element_text(size=22,  family="Helvetica", colour="black")) +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)
dev.off()


# Visualize the distribution of mitochondrial gene expression detected per cell
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_MTgenesvscells.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
	ggtitle("Mitochondrial gene expression per cell") + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
	theme(text = element_text(size=22,  family="Helvetica", colour="black")) +
  	geom_vline(xintercept = 0.2)
dev.off()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
png(paste0(out_dir,"QC/After_filter/",out_prefix,"_complexity.png"), width = 14, height = 12, units = 'in', res = 300)    
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
	ggtitle("Overall complexity (Good > 0.8)") + 
  	geom_density(alpha = 0.2) +
  	theme_classic() +
	theme(text = element_text(size=22,  family="Helvetica", colour="black")) +
  	geom_vline(xintercept = 0.8)
dev.off()
	

######################################################
