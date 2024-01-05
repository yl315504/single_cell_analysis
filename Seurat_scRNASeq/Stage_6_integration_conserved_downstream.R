# Rscript edgeR.R $COUNTS_FILE $SAMPLE_FILE $OUT_DIR $EDGER_BCV $CONTRAST_FILE
rm(list=ls())
options(echo=TRUE)
# usage:
usage = "Usage: Rscript script.R  <in_dir> <out_prefix> <cluster res>\n"

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


#in_dir = "/scratch1/01775/saathe/DG_china_cohort/"
#SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
#out_prefix = "Pat_0012"
#res = "0.4"
#res_string <- paste0("integrated_snn_res.",res)

in_dir = args[1]
#SamplesSheet_file = args[2]
SamplesSheet_file = paste0(in_dir,"/Seurat_input_data/SamplesSheet.csv")
out_prefix = args[2]
res = args[3]
res_string <- paste0("integrated_snn_res.",res)

in_dir
out_prefix
SamplesSheet_file
res
res_string

setwd(in_dir)
#in_dir <- "D:/Data/HY_Yapeng/CUS47_10671/"
out_dir <- paste0(in_dir,"/Seurat_output/data/results/")
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

# Assign identity of clusters
Idents(object = seurat_integrated) <- res_string

#seurat_integrated$Group <- NA
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "${sample1}"))] <- "${sample1}"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "${sample2}"))] <- "${sample2}"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "R159C"))] <- "R159C"


#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "CASE"))] <- "CASE"
#seurat_integrated$Group[which(str_detect(seurat_integrated$orig.ident, "CTRL"))] <- "CTRL"

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

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

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
 
# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

#Exploring known cell type markers
#png(paste0(out_dir,"Integrated_data_figures/",out_prefix,"_dimplot-markers.png"), width = 12, height = 12, units = 'in', res = 300)                           
#DimPlot(object = seurat_integrated, 
#        reduction = "umap", 
#        label = TRUE) + NoLegend()
		
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
##################################################################

sel.clust = res_string
##################################################################

# Create function to get conserved markers for any given cluster
#Please delete min.cells.group = 0, when done

Per_cluster_conserved <- function(cluster){
        FindConservedMarkers(seurat_integrated,
                             ident.1 = cluster,
                             grouping.var = "sample", min.cells.group = 0,
                             only.pos = TRUE) %>%
                rownames_to_column(var = "gene") %>%
                cbind(cluster_id = cluster, .)
}

#map_dfr(inputs_to_function, name_of_function)

# Iterate function across desired clusters
cons <- map_dfr(levels(seurat_integrated@active.ident), Per_cluster_conserved)

Conserved_genes_per_cluster <- table(cons$cluster_id)
tmp <- paste0(out_dir,"Integrated_data_figures/",out_prefix,"_conserved_markers.txt")
write.table(cons,tmp,sep="\t",row.names=TRUE)


# Extract top 10 markers per cluster
top25 <- cons %>% 
  mutate(avg_fc = (Y284L_avg_log2FC+Y284R_avg_log2FC)/2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 25, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
#View(top25)  

tmp <- paste0(out_dir,"Integrated_data_figures/",out_prefix,"_samples_conservedgenes_barplot.pdf")
pdf(tmp, width = 14, height = 14)
par(1, 5, mar = c(4, 6, 3, 1))
par(mfrow=c(3,6))
lapply(unique(top25$cluster_id), 
function(i) {
barplot(sort(setNames(((top25$Y284L_avg_log2FC + top25$Y284R_avg_log2FC)/2), top25$gene)[top25$cluster_id == i], F), horiz = T,las = 1, xlab = "log2(Fold Change)", main = paste0(i, " conserved:Y284L+Y284R"), xlim = c(0,as.numeric(max(sort(setNames(((top25$Y284L_avg_log2FC + top25$Y284R_avg_log2FC)/2), top25$gene)[top25$cluster_id == i], F)))), border = "white", yaxs = "i"); abline(v = c(0, 0.25), lty = c(1, 2))
})
dev.off()

top5 <- cons %>% 
  mutate(avg_fc = (Y284L_avg_log2FC + Y284R_avg_log2FC)/2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 5, 
        wt = avg_fc)


#cols <- viridis(100)[c(1, 50, 100)]

#tmp <- paste0(out_dir,"Integrated_data_figures/",out_prefix,"_samples_conservedgenes_heatmap.pdf")
#pdf(tmp, width = 14, height = 14)

##DoHeatmap(seurat_integrated, genes.use = top5$gene, slim.col.label = TRUE,
#          remove.key = TRUE, col.low = cols[1], col.mid = cols[2],
#          col.high = cols[3])
#DoHeatmap(seurat_integrated, features = top5$gene, group.by = sel.clust, assay = "RNA",col.high = cols[3])

#dev.off()

# Save integrated seurat object
saveRDS(seurat_integrated, paste0(out_dir,"/",out_prefix,"_final_integrated_seurat.rds"))
