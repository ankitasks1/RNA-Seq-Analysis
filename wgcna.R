setwd("/Users/ankitverma/Documents/tutorial/wgcna")
# Uncomment and modify the following to install any missing packages
# install.packages(c("tidyverse", "magrittr", "WGCNA))
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)  
# ==== Load and clean data
data <- data.frame(data.table::fread("GSE61333_ligule_count.txt", header = T))
#> ── Column specification ────────────────────────────────────────────────────────
#> cols(
#>   .default = col_double(),
#>   Count = col_character()
#> )
#> ℹ Use `spec()` for the full column specifications.

data[1:5,1:10]        # Look at first 5 rows and 10 columns
rownames(data) <- data$Count
data <- data[,-1]
dim(data)
tdata <- t(data)
de_input = as.matrix(data)
library(DESeq2)

coldata <- data.frame(Sample = colnames(de_input),
                      Type=factor(sapply(strsplit(colnames(de_input), "\\."), function(x) x[1])))
rownames(coldata) <- coldata$Sample

dds <- DESeqDataSetFromMatrix(round(de_input),
                              coldata,
                              design = ~Type)

dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)
library(genefilter)      # <= why is this here?
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
expr_normalized[1:5,1:10]

# convert longer form for plotting by ggplot2
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

# explore the normalized counts by violin plot
expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )



# WGCNA
input_mat = t(expr_normalized)

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns
# We can see now that the rows = treatments and columns = gene probes. We’re ready to start WGCNA. A correlation network will be a complete network (all genes are connected to all other genes). Ergo we will need to pick a threshhold value (if correlation is below threshold, remove the edge). We assume the true biological network follows a scale-free structure 
#library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)
# steps correlation->adjacency (cor**power)->Topological Overlap Measure -> clusters (with the DynamicTree algorithm)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(12,20,2))
# The soft thresholding, is a value used to power the correlation of the genes to that threshold. The assumption on that by raising the correlation to a power will reduce the noise of the correlations in the adjacency matrix

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

# plot asesses
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


# TOM : This matrix is also a similarity matrix but this one is built on the amount of genes that two genes have in common, therefore the overlap of neighbors between two genes.
# It calculates a correlation matrix from the expression data, calculates a soft threshold and assigns two genes a high topological overlap if they share common neighbourhoods.
# Note: Pick a soft threshold power near the curve of the plot, so here we could pick 7, 8 or 9. We’ll pick 9 but feel free to experiment with other powers to see how it affects your results. Now we can create the network using the blockwiseModules command. The blockwiseModule may take a while to run, since it is constructing the TOM (topological overlap matrix) and several other steps. 
picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

![dendrogram](RNA-Seq-Analysis/Rplot01.png)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]


write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")



# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
# Looking at the heatmap, the green seems positively associated (red shading) with the L groups, turquoise module is positive for the L but negative (blue shading) for the L_L1 groups, and tan module is positive in the S groups but negative elsewhere. There are other patterns, think about what the patterns mean. For this tutorial we will focus on the green, turquoise, and tan module genes.

# Examine Expression Profiles

# pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


# It is possible to plot the modules in green, tan, and turquoise colors but for now we’ll use the auto colors generated by ggplot2. From the expression line graph, the line graph matches the correlation graph where Green module genes are higher in the L groups. Tan genes seem higher in the S groups, and Turquoise genes are low at L_L1 and S_L1 and maybe B_L1 groups.
# Generate and Export Networks
#The network file can be generated for Cytoscape or as an edge/vertices file.
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)
# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")
