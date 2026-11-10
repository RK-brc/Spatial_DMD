library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(viridis)
library(ggrepel)
library(scales)
library(RColorBrewer)


# ------Fig.5B Circle plot for niches---------

plot_cellchat_style_circles <- function(df,
                                        validation_filter = NULL,  
                                        min_score = 0.05,          #lrscore
                                        score_column = "lrscore",  
                                        save_file = FALSE,
                                        filename = "cellchat_style_circles") {
  
  library(igraph)
  library(RColorBrewer)
  library(dplyr)
  
  filtered <- df
  
  #validation filter
  if(!is.null(validation_filter) && "validation_status" %in% colnames(df)) {
    filtered <- filtered %>%
      filter(validation_status %in% validation_filter)
  }
  
  #filter by score 
  if(score_column %in% colnames(filtered)) {
    filtered <- filtered %>%
      filter(.data[[score_column]] >= min_score)
  } else {
    warning(paste("Column", score_column, "not found "))
  }
  
  #filter by significance 
  if("is_significant" %in% colnames(filtered)) {
    filtered <- filtered %>%
      filter(is_significant == TRUE)
  }
  
  if(nrow(filtered) == 0) {
    stop("No interactions pass the filters")
  }
  
  #aggregate:number of unique LR pairs per source-target
  count_data <- filtered %>%
    group_by(source, target) %>%
    summarise(count = n(), .groups = "drop") %>%
    rename(source_niche = source, target_niche = target)
  
  #aggregate: Weight (sum of lrscore)
  weight_data <- filtered %>%
    group_by(source, target) %>%
    summarise(weight = sum(.data[[score_column]], na.rm = TRUE), .groups = "drop") %>%
    rename(source_niche = source, target_niche = target)
  
  #get all unique niches
  all_niches <- sort(unique(c(filtered$source, filtered$target)))
  n_niches <- length(all_niches)
  
  #colors
  if(n_niches <= 9) {
    node_colors <- setNames(brewer.pal(max(3, n_niches), "Set1")[1:n_niches], all_niches)
  } else {
    node_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(n_niches), all_niches)
  }
  
  #create circle plot
  plot_circle_network <- function(edge_df, weight_col, title_text, scale_factor = 1) {
    #create graph
    g <- graph_from_data_frame(edge_df[, c("source_niche", "target_niche", weight_col)],
                               directed = TRUE,
                               vertices = data.frame(name = all_niches))
    
    #edge weights
    edge_weights <- edge_df[[weight_col]]
    
    #calculate node sizes based on total outgoing interactions
    node_out_strength <- sapply(all_niches, function(niche) {
      sum(edge_df[[weight_col]][edge_df$source_niche == niche], na.rm = TRUE)
    })
    
    #normalize node sizes
    if(max(node_out_strength) > min(node_out_strength)) {
      node_sizes <- 15 + 35 * (node_out_strength - min(node_out_strength)) /
        (max(node_out_strength) - min(node_out_strength))
    } else {
      node_sizes <- rep(30, length(node_out_strength))
    }
    
    #edge widths
    if(max(edge_weights) > min(edge_weights)) {
      edge_widths <- 1 + 9 * (edge_weights - min(edge_weights)) /
        (max(edge_weights) - min(edge_weights))
    } else {
      edge_widths <- rep(5, length(edge_weights))
    }
    
    #edge transparency based on weight
    if(max(edge_weights) > min(edge_weights)) {
      edge_alpha <- 0.3 + 0.6 * (edge_weights - min(edge_weights)) /
        (max(edge_weights) - min(edge_weights))
    } else {
      edge_alpha <- rep(0.6, length(edge_weights))
    }
    
    #edge colors
    edge_colors <- sapply(1:nrow(edge_df), function(i) {
      adjustcolor(node_colors[edge_df$source_niche[i]], alpha.f = edge_alpha[i])
    })
    
    #circular layout
    layout_circle <- layout_in_circle(g)
    
    #plot
    plot(g,
         layout = layout_circle,
         vertex.size = node_sizes,
         vertex.color = node_colors[V(g)$name],
         vertex.frame.color = "white",
         vertex.label = V(g)$name,
         vertex.label.color = "black",
         vertex.label.cex = 0.9,
         vertex.label.dist = 2.5,
         vertex.label.degree = -pi/2,
         edge.width = edge_widths * scale_factor,
         edge.arrow.size = 0.5,
         edge.arrow.width = 1.5,
         edge.color = edge_colors,
         edge.curved = 0.2,
         main = title_text,
         margin = c(0, 0, 0, 0))
  }
  
  #create plots
  if(save_file) {
    ext <- tolower(tools::file_ext(filename))
    if(ext == "") {
      filename <- paste0(filename, ".pdf")
      ext <- "pdf"
    }
    
    if(ext == "pdf") {
      pdf(filename, width = 16, height = 8)
    } else if(ext == "tiff") {
      tiff(filename, width = 16, height = 8, units = "in", res = 300)
    } else if(ext == "png") {
      png(filename, width = 16, height = 8, units = "in", res = 300)
    }
  }
  
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  
  #Plot 1: Number of interactions
  plot_circle_network(count_data, "count", "Number of interactions", scale_factor = 1)
  
  #Plot 2: Interaction weights/strength
  plot_circle_network(weight_data, "weight", 
                      paste0("Interaction strength (", score_column, ")"), 
                      scale_factor = 0.8)
  
  if(save_file) {
    dev.off()
    cat("Saved to:", filename, "\n")
  }
  
  par(mfrow = c(1, 1))
  
  invisible(list(
    n_interactions = nrow(filtered),
    n_niches = n_niches,
    count_summary = count_data,
    weight_summary = weight_data
  ))
}

#plot using lrscore with default threshold
plot_cellchat_style_circles(
  liana_rank_aggregate_all_samples,
  min_score = 0.05,
  score_column = "lrscore",
  save_file = TRUE,
  filename = "overallniche_communication_overview.pdf"
)


#---- Supplementary Fig.8A Heatmap showing interaction counts between all niche pairs--------


plot_niche_interaction_heatmap <- function(df,
                                           min_score = 0.05,
                                           score_column = "lrscore",
                                           value_type = "count") {
  
  filtered <- df %>%
    filter(.data[[score_column]] >= min_score)
  
  if("is_significant" %in% colnames(filtered)) {
    filtered <- filtered %>% filter(is_significant == TRUE)
  }
  
  #aggregate
  if(value_type == "count") {
    heatmap_data <- filtered %>%
      group_by(source, target) %>%
      summarise(value = n(), .groups = "drop")
    legend_name <- "# Interactions"
  } else {
    heatmap_data <- filtered %>%
      group_by(source, target) %>%
      summarise(value = sum(.data[[score_column]], na.rm = TRUE), .groups = "drop")
    legend_name <- paste0("Total ", score_column)
  }
  
  #create matrix
  all_niches <- sort(unique(c(filtered$source, filtered$target)))
  heatmap_matrix <- heatmap_data %>%
    complete(source = all_niches, target = all_niches, fill = list(value = 0)) %>%
    pivot_wider(names_from = target, values_from = value, values_fill = 0)
  
  #calculate row sums. this is total outgoing
  row_totals <- heatmap_matrix %>%
    select(-source) %>%
    rowSums()
  
  heatmap_matrix$total_outgoing <- row_totals
  
  plot_data <- heatmap_matrix %>%
    select(-total_outgoing) %>%
    pivot_longer(-source, names_to = "target", values_to = "value")
  
  #create heatmap
  p <- ggplot(plot_data, aes(x = target, y = source, fill = value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = ifelse(value > 0, round(value, 1), "")), 
              color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "#E41A1C", name = legend_name) +
    labs(title = "Niche-to-Niche Communication Matrix",
         subtitle = paste0(value_type, " of interactions (min ", score_column, " = ", min_score, ")"),
         x = "Target Niche (Receiver)",
         y = "Source Niche (Sender)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 14),
          panel.grid = element_blank())
  
  #adding marginal totals
  sender_totals <- heatmap_data %>%
    group_by(source) %>%
    summarise(total = sum(value)) %>%
    arrange(desc(total))
  
  print("Top Senders:")
  print(sender_totals)
  
  return(p)
}


p3 <- plot_niche_interaction_heatmap(liana_rank_aggregate_all_samples, 
                                     min_score = 0.05, 
                                     score_column = "lrscore",
                                     value_type = "count")
ggsave("niche_interaction_heatmap.pdf", p3, width = 10, height = 8)


## Fig.5c,d

#input
master_expanded_clean <- read.csv("~/rasya/spatial/Liana_2026March/misty_results/master_expanded_clean.csv", row.names = 1)


#FOR THE BUBBLE PLOT FOR NICHES

plot_niche_bubble <- function(df, source_niche, target_niches = NULL,
                              ligands_of_interest = NULL,
                              validation_filter = c("High confidence", "Medium confidence"),
                              min_score = 0.05,
                              top_n = 30) {
  
  library(dplyr)
  library(ggplot2)
  
  #3 receptor cell type tiers
  expanded_df <- bind_rows(
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype,
        receptor_celltype_pct_selected = receptor_celltype_pct,
        receptor_present_selected = receptor_present_in_target,
        receptor_tier = "primary"
      ),
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype_2nd,
        receptor_celltype_pct_selected = receptor_celltype_2nd_pct,
        receptor_present_selected = receptor_celltype_2nd_present,
        receptor_tier = "secondary"
      ),
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype_3rd,
        receptor_celltype_pct_selected = receptor_celltype_3rd_pct,
        receptor_present_selected = receptor_celltype_3rd_present,
        receptor_tier = "tertiary"
      )
  ) %>%
    filter(!is.na(receptor_celltype_selected),
           receptor_celltype_selected != "")
  
  #filtering by source niche
  filtered <- expanded_df %>%
    filter(source_niche == !!source_niche,
           grepl(paste(validation_filter, collapse = "|"), validation_status_extended),
           mean_score >= min_score)
  
  if(!is.null(target_niches)) {
    filtered <- filtered %>%
      filter(target_niche %in% target_niches)
  }
  
  #ligands_of_interest specified
  if(!is.null(ligands_of_interest)) {
    filtered <- filtered %>%
      filter(ligand %in% ligands_of_interest)
  }
  
  if(nrow(filtered) == 0) {
    stop("No interactions found")
  }
  
  #selecting top N by mean_score
  plot_data <- filtered %>%
    mutate(LR_pair = paste0(ligand, " - ", receptor),
           interaction = paste0(source_niche, " → ", target_niche)) %>%
    arrange(desc(mean_score)) %>%
    head(top_n)
  
  if(nrow(plot_data) == 0) {
    stop("No rows available for plotting")
  }
  
  effect_levels <- unique(as.character(plot_data$best_effect_direction))
  effect_palette <- setNames(
    c("red", "lightblue")[seq_along(effect_levels)],
    effect_levels
  )
  
  #plot
  p <- ggplot(plot_data, aes(x = target_niche, y = reorder(LR_pair, mean_score))) +
    geom_point(aes(size = mean_score, color = best_effect_direction), alpha = 0.8) +
    scale_size_continuous(range = c(2, 12),
                          name = "L-R\nexpression\n(mean_score)") +
    scale_color_manual(values = effect_palette,
                       na.value = "gray80",
                       name = "Cell type\nneighbourhood effect") +
    labs(title = paste0(source_niche, " → Target Niches"),
         x = "",
         y = "Ligand - Receptor") +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.x = element_blank())
  
  return(p)
}

#for niche1
my_ligands <- c("EFEMP2", "LAMA2", "COL6A1", "COL3A1", "DCN", 
                "COL1A2", "C3", "COL1A1","SERPING1")

#for niche0
my_ligands <- c("C1QB", "CALR", "SERPING1", "COL1A1", 
                "C3", "HLA-DRA ", "HLA-DRB1", "FN1", "HLA-DPA1")

p <- plot_niche_bubble(
  master_expanded_clean,
  source_niche = "Niche1",
  #target_niches = c("Niche3", "Niche4"),
  ligands_of_interest = my_ligands,
  validation_filter = c("High confidence", "Medium confidence"),
  min_score = 0.05,
  top_n = 200
)

tiff("niche1_to_all_targets.tiff",height=18,width=18,units='cm',compression = "lzw",res=300)
p
dev.off()


# for supplementary fig. 8b-d

# celltype to target cell types bubble plot

plot_celltype_to_all_bubble <- function(df, 
                                        source_celltype,
                                        target_celltypes = NULL,
                                        lr_pairs_of_interest = NULL,
                                        ligand_pattern = NULL,
                                        receptor_pattern = NULL,
                                        validation_filter = c("High confidence", "Medium confidence"),
                                        min_score = 0.05) {
  
  library(dplyr)
  library(stringr)
  library(ggplot2)
  
  #include all 3 receptor cell type tiers
  expanded_df <- bind_rows(
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype,
        receptor_celltype_pct_selected = receptor_celltype_pct,
        receptor_present_selected = receptor_present_in_target,
        receptor_tier = "primary"
      ),
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype_2nd,
        receptor_celltype_pct_selected = receptor_celltype_2nd_pct,
        receptor_present_selected = receptor_celltype_2nd_present,
        receptor_tier = "secondary"
      ),
    df %>%
      mutate(
        receptor_celltype_selected = receptor_celltype_3rd,
        receptor_celltype_pct_selected = receptor_celltype_3rd_pct,
        receptor_present_selected = receptor_celltype_3rd_present,
        receptor_tier = "tertiary"
      )
  ) %>%
    filter(!is.na(receptor_celltype_selected),
           receptor_celltype_selected != "")
  
  #filter for source cell type sending signals
  filtered <- expanded_df %>%
    filter(ligand_celltype == source_celltype,
           grepl(paste(validation_filter, collapse = "|"), validation_status_extended),
           mean_score >= min_score) %>%
    mutate(LR_pair = paste0(ligand, " - ", receptor))
  
  #filter selected receiving cell types
  if(!is.null(target_celltypes)) {
    filtered <- filtered %>%
      filter(receptor_celltype_selected %in% target_celltypes)
  }
  
  #LR pair filters
  if(!is.null(lr_pairs_of_interest)) {
    filtered <- filtered %>%
      filter(ligand %in% lr_pairs_of_interest | 
               receptor %in% lr_pairs_of_interest |
               LR_pair %in% lr_pairs_of_interest |
               paste0(ligand, "-", receptor) %in% lr_pairs_of_interest)
  }
  
  if(!is.null(ligand_pattern)) {
    filtered <- filtered %>%
      filter(str_detect(ligand, ligand_pattern))
  }
  
  if(!is.null(receptor_pattern)) {
    filtered <- filtered %>%
      filter(str_detect(receptor, receptor_pattern))
  }
  
  if(nrow(filtered) == 0) {
    stop("No interactions found with the specified filters")
  }
  
  #colors
  effect_colors <- c("Promotes accumulation" = "#E41A1C",
                     "Promotes depletion" = "#377EB8",
                     "Mixed/Unclear" = "#FFD700")
  
  filter_info <- c()
  if(!is.null(target_celltypes)) filter_info <- c(filter_info, paste0("targets:", paste(target_celltypes, collapse = ", ")))
  if(!is.null(lr_pairs_of_interest)) filter_info <- c(filter_info, "specific LR pairs")
  if(!is.null(ligand_pattern)) filter_info <- c(filter_info, paste0("ligand:", ligand_pattern))
  if(!is.null(receptor_pattern)) filter_info <- c(filter_info, paste0("receptor:", receptor_pattern))
  
  subtitle_text <- ifelse(length(filter_info) > 0,
                          paste0(nrow(filtered), " interactions (", paste(filter_info, collapse = ", "), ")"),
                          paste0(nrow(filtered), " interactions"))
  
  #plot
  p <- ggplot(filtered, aes(x = receptor_celltype_selected, y = reorder(LR_pair, mean_score))) +
    geom_point(aes(size = mean_score, color = best_effect_direction), alpha = 0.8) +
    scale_size_continuous(range = c(3, 12), 
                          name = "Interaction\nScore") +
    scale_color_manual(values = effect_colors, 
                       name = "Functional\nEffect",
                       na.value = "gray70") +
    labs(title = paste0(source_celltype, " → All Target Cell Types"),
         subtitle = subtitle_text,
         x = "Receiving Cell Type",
         y = "Ligand - Receptor") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 9),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.x = element_blank())
  
  return(p)
}

my_ligands <-  c("C1QB", "CALR", "PSAP", "SERPING1", "COL1A1", 
                 "C3", "HLA-DRA ", "HLA-DRB1", "FN1", "HLA-DPA1")

tiff("InflammatoryFAPstoall_plot.tiff", 
     height=20, width=20, units="cm", 
     compression = "lzw", res=300)

plot_celltype_to_all_bubble(
  master_expanded,
  source_celltype = "Inflammatory FAPs",
  target_celltypes = c("Inflammatory Cells","Pre-Adipogenic", "Inflammatory FAPs", "Endothelial Cells","Satellite Cells","MME+ FAPs","Smooth Muscle Cells"),
  lr_pairs_of_interest = my_ligands,
  min_score = 0.05
)

dev.off()

tiff("Inflammatorycellstoall_plot.tiff", 
     height=20, width=20, units="cm", 
     compression = "lzw", res=300)

plot_celltype_to_all_bubble(
  master_expanded,
  source_celltype = "Inflammatory Cells",
  target_celltypes = c("Inflammatory Cells","Pre-Adipogenic", "Inflammatory FAPs", "Endothelial Cells","Satellite Cells","MME+ FAPs","Smooth Muscle Cells"),
  lr_pairs_of_interest = my_ligands,
  min_score = 0.05
)

dev.off()

my_ligands <- c("EFEMP2", "LAMA2", "COL6A1", "COL3A1", "DCN", 
                "COL1A2", "MMP2", "COL1A1","SERPING1")

tiff("MMEFAPstoall_plot.tiff", 
     height=20, width=20, units="cm", 
     compression = "lzw", res=300)

plot_celltype_to_all_bubble(
  master_expanded,
  source_celltype = "MME+ FAPs",
  target_celltypes = c("Inflammatory Cells","Pre-Adipogenic", "Inflammatory FAPs", "Endothelial Cells","Satellite Cells","MME+ FAPs","Smooth Muscle Cells"),
  lr_pairs_of_interest = my_ligands,
  min_score = 0.03
)

dev.off()




