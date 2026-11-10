# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichR)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(reshape2)

df <- read.csv("4. Supplementary_Table_all_niche_sig_degs.csv")  


#filter significant genes
sig_genes <- df %>% 
  filter(pval_adj < 0.05, abs(logFC) > 0.5)

#check available databases
dbs <- listEnrichrDbs()
print(head(dbs))

#choose databases 
databases <- c("GO_Biological_Process_2023", 
               "KEGG_2021_Human", 
               "MSigDB_Hallmark_2020")

#get unique niches
niches <- unique(sig_genes$Niche)

#run enrichment for each niche
enrichr_results <- list()
for (niche in niches) {
  genes <- sig_genes %>% 
    filter(Niche == niche) %>% 
    pull(gene)
  
  print(paste("Processing", niche, "with", length(genes), "genes"))
  
  if (length(genes) > 0) {
    enrichr_results[[niche]] <- enrichr(genes, databases)
  }
}

#combine results into a single dataframe 
all_results <- data.frame()
for (niche in names(enrichr_results)) {
  for (db in databases) {
    if (!is.null(enrichr_results[[niche]][[db]]) && 
        nrow(enrichr_results[[niche]][[db]]) > 0) {
      temp <- enrichr_results[[niche]][[db]]
      temp$Niche <- niche
      temp$Database <- db
      all_results <- rbind(all_results, temp)
    }
  }
}

#filter by adjusted p-value
enrichment_filtered <- all_results %>%
  filter(Adjusted.P.value <= 0.05)

print(paste("Total significant pathways found:", nrow(enrichment_filtered)))

#output directory
dir.create("pathway_plots", showWarnings = FALSE)

write.csv(enrichment_filtered, "pathway_plots/all_enrichment_results.csv", row.names = FALSE)

#for visualisation

#data with pathways of interest
pathway_data <- data.frame(
  Niche = c(
    rep("Niche 0", 7),
    rep("Niche 1", 8),
    rep("Niche 2", 6),
    rep("Niche 3", 1),
    rep("Niche 4", 8),
    rep("Niche 5", 9),
    rep("Niche 6", 9)
  ),
  Pathway = c(
    # Niche 0
    "Phagosome", "Allograft Rejection", "Interferon Gamma Response", "Complement", 
    "Lysosome", "Antigen processing and presentation", "Regulation of T cell activation",
    # Niche 1
    "Epithelial Mesenchymal Transition", "Extracellular Matrix Organization", 
    "Elastic Fiber Assembly", "Collagen Fibril Organization", "Extracellular Matrix Assembly",
    "Regulation of ERK1/2 Cascade", "Myogenesis", "Coagulation",
    # Niche 2
    "Positive Regulation of Macrophage Chemotaxis", "Positive Regulation of Macrophage Migration",
    "Regulation of Granulocyte Chemotaxis", "Regulation of Macrophage Chemotaxis",
    "Regulation of Mononuclear Cell Migration", "Positive Regulation of Leukocyte Chemotaxis",
    # Niche 3
    "Myogenesis",
    # Niche 4
    "Myogenesis", "TNF-alpha Signaling via NF-kB", "Striated Muscle Contraction",
    "Muscle Filament Sliding", 
    "Actin-Myosin Filament Sliding", "Myofibril Assembly", "Sarcomere Organization",
    # Niche 5
    "Focal Adhesion", "Vascular Smooth Muscle Contraction", "Apical Junction",
    "Regulation of Cell Migration", "Epithelial Mesenchymal Transition", 
    "Regulation of Actin Cytoskeleton", "Apelin Signaling", "Regulation of Angiogenesis",
    "Positive Regulation of Cell Motility",
    # Niche 6
    "Muscle Contraction", "Myogenesis", "Glycolysis/Gluconeogenesis", "Glucagon Signaling",
    "HIF-1 Signaling", "Hypoxia", "Glycolytic Process", 
    "Carbohydrate Catabolic Process","Calcium signaling pathway"
  ),
  P_value = c(
    # Niche 0
    8.334940e-26
    , 7.463712e-25, 4.404419e-21, 2.814486e-20, 8.441111e-19, 3.661590e-18, 4.230867e-15,
    # Niche 1
    4.422066e-41, 4.854021e-10, 8.133722e-08, 3.423963e-07, 5.828427e-07,1.319354e-06,2.321132e-05, 1.210485e-08,
    # Niche 2
    1.965878e-07, 3.317395e-07, 3.317395e-07, 3.317395e-07, 4.002435e-07, 1.416435e-05,
    # Niche 3
    3.096751e-07,
    # Niche 4
    7.227284e-12, 1.289513e-09, 2.309155e-09, 1.366641e-07, 1.751203e-07, 1.751203e-07, 3.344857e-07,
    # Niche 5
    9.784915e-18, 7.329159e-13, 1.727334e-12, 4.167525e-12, 6.081908e-11, 1.217900e-09, 7.654946e-09, 3.652348e-08, 1.004788e-07,
    # Niche 6
    3.230250e-07, 7.343146e-07, 1.974568e-06, 1.974568e-06, 1.974568e-06, 5.127645e-05, 5.445705e-05,  2.053444e-04,6.292263e-05
  )
)

#calculate -log10(FDR)
pathway_data$neg_log10_FDR <- -log10(pathway_data$P_value)

#create matrix for heatmap
mat <- pathway_data %>%
  dplyr::select(Pathway, Niche, neg_log10_FDR) %>%
  tidyr::pivot_wider(names_from = Niche, values_from = neg_log10_FDR, values_fill = 0) %>%
  as.data.frame()

rownames(mat) <- mat$Pathway
mat <- mat[, -1]

#long format
mat_long <- reshape2::melt(as.matrix(mat))
colnames(mat_long) <- c("Pathway", "Niche", "value")

write.csv(mat_long, file = "niche_functional_finalheatmap.csv")
write.csv(distinct_pathways, file = "niche_functional_distinctpathways.csv")

#heatmap
p <- ggplot(mat_long, aes(x = Niche, y = Pathway, fill = value)) +
  geom_tile(color = "gray", linewidth = 0.8) +  
  scale_fill_gradient(low = "white", high = "darkred", 
                      name = "-log10(FDR)",
                      guide = guide_colorbar(barwidth = 0.8, barheight = 8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               size = 11, color = "black"),  
    axis.text.y = element_text(size = 10, color = "black"), 
    axis.title = element_text(size = 12, color = "black"),  
    plot.title = element_text(hjust = 0.5, margin = margin(b = 20), 
                              size = 14, color = "black"),    
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10, color = "black"),  
    legend.title = element_text(size = 11, color = "black") 
  ) +
  labs(
    x = "Niche",
    y = "Pathway"
  )

tiff("niche_functional_finalheatmap.tiff",height=20,width=20,units='cm',compression = "lzw",res=300)
p
dev.off()


