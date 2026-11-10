library(tidyverse)

rank_agg    <- read.csv("rank_aggr_forNiche0n1.csv", row.names = 1) #filtered data with source=Niche0 and source=Niche1 from liana_rank_aggregate_all_samples.csv
robust_para <- read.csv("A_robust_para_interactions.csv", row.names = 1)
lr_robust   <- read.csv("B_lr_robust_interactions.csv",row.names = 1)
overlap     <- read.csv("rankagg_misty_overlap.csv",row.names = 1)

# ------------- cross reference rank aggregate and overlap ---------------------

#join overlap back to rank aggregate to get scores
overlap_scored <- overlap %>%
  inner_join(rank_agg, 
             by = c("ligand" = "ligand_complex", 
                    "receptor" = "receptor_complex")) %>%
  select(ligand, receptor, source, target, score, 
         magnitude_rank, specificity_rank, donor) %>%
  group_by(ligand, receptor, source, target) %>%
  summarise(
    mean_score        = mean(score, na.rm = TRUE),
    mean_magnitude    = mean(magnitude_rank, na.rm = TRUE),
    mean_specificity  = mean(specificity_rank, na.rm = TRUE),
    n_donors          = n_distinct(donor),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

print("Cross-ref 1: Overlap pairs scored from rank aggregate")
print(overlap_scored, n = 30)

# ----------- Rank aggregate ligands → MISTy A cell type directionality --------

#get unique ligands from rank aggregate
rank_ligands <- rank_agg %>%
  distinct(ligand_complex) %>%
  pull(ligand_complex)

#find these ligands in MISTy A robust_para
directionality <- robust_para %>%
  filter(predictor %in% rank_ligands) %>%
  select(predictor, target, mean_importance, 
         n_donors, n_donors_positive, n_donors_negative) %>%
  arrange(desc(mean_importance))

print("Cross-ref 2: Rank aggregate ligands that predict cell type abundance (MISTy A)")
print(directionality, n = 30)

#from directionality filter only the ones that predict inflammatory faps and cells and MME+ faps

niche_directionality <- directionality %>% filter(target %in% c("meanscell_abundance_w_sf_Inflammatory.FAPs", "meanscell_abundance_w_sf_Inflammatory.Cells",
                                                                'meanscell_abundance_w_sf_MME..FAPs'))


#----------------- snRNA-seq data ----------------------------------------------------------------

#get ligand and receptor expression per cell type from snRNA-seq (the final DMD+faps subtypes object)

dmd <- readRDS("~/dmd.integrated_withfaps.rds")

#unique ligands and receptors from overlap table
ligands <- unique(overlap$ligand)
receptors <- unique(overlap$receptor)
all_genes <- unique(c(ligands, receptors))

#average expression per cell type
DefaultAssay(dmd)
avg_exp <- AverageExpression(
  dmd, 
  features = all_genes,
  group.by = "merge_type", 
)

avg_exp <- avg_exp$RNA

#for each gene, find the cell type with max expression
gene_celltype_map <- data.frame(
  gene = rownames(avg_exp),
  top_celltype = colnames(avg_exp)[apply(avg_exp, 1, which.max)],
  max_expression = apply(avg_exp, 1, max)
)

#add the ligand assignments
overlap$ligand_celltype <- gene_celltype_map$top_celltype[match(overlap$ligand, gene_celltype_map$gene)]

#receptor assignments  
overlap$receptor_celltype <- gene_celltype_map$top_celltype[match(overlap$receptor, gene_celltype_map$gene)]

write.csv(overlap, "rankagg_mistyoverlap_with_celltypes.csv", row.names = FALSE)


## combinging everything together
master_table <- overlap_scored %>%
  left_join(overlap %>% 
              select(ligand, receptor, ligand_celltype, receptor_celltype),
            by = c("ligand", "receptor")) %>%
  rename(source_niche = source,
         target_niche = target) %>%
  arrange(desc(mean_score))

#save
write.csv(master_table, "master_communication_table.csv", row.names = FALSE)


##combining the above with niche composition
niche_comp <- read.csv("~/niche_composition.csv", row.names = 1)

colnames(niche_comp) <- gsub("\\.+", " ", colnames(niche_comp))
colnames(niche_comp)[colnames(niche_comp) == "MME FAPs"] <- "MME+ FAPs"
colnames(niche_comp)[colnames(niche_comp) == "Pre Adipogenic"] <- "Pre-Adipogenic"

#get cell type percentage in a niche
get_celltype_pct <- function(niche, celltype, comp_df) {
  if (niche %in% rownames(comp_df)) {
    if (celltype %in% colnames(comp_df)) {
      return(comp_df[niche, celltype])
    }
  }
  return(0)
}

#add niche composition percentages to master_table
master_table <- master_table %>%
  rowwise() %>%
  mutate(
    ligand_celltype_pct = get_celltype_pct(source_niche, ligand_celltype, niche_comp),
    receptor_celltype_pct = get_celltype_pct(target_niche, receptor_celltype, niche_comp)
  ) %>%
  ungroup() %>%
  #add validation flags (threshold >5% presence)
  mutate(
    ligand_present_in_source = ligand_celltype_pct > 5,
    receptor_present_in_target = receptor_celltype_pct > 5,
    high_confidence = ligand_present_in_source & receptor_present_in_target
  ) %>%
  arrange(desc(mean_score))


#adding confidence
master_table <- master_table %>%
  mutate(
    validation_status = case_when(
      ligand_celltype_pct > 10 & receptor_celltype_pct > 10 ~ "High confidence",
      ligand_celltype_pct > 5 & receptor_celltype_pct > 5 ~ "Medium confidence",
      ligand_celltype_pct > 1 | receptor_celltype_pct > 1 ~ "Low abundance sender/receiver",
      TRUE ~ "Check expression"
    )
  )


## ------------- adding everything with Misty A results ------------------

#add directionality
misty_a <- robust_para %>%
  filter(view == "para") %>%  
  mutate(
    effect_direction = case_when(
      n_donors_positive > n_donors_negative ~ "Promotes accumulation",
      n_donors_negative > n_donors_positive ~ "Promotes depletion",
      TRUE ~ "Mixed/Unclear"
    ),
    net_donors = n_donors_positive - n_donors_negative 
  ) %>%
  rename(ligand = predictor,
         target_celltype = target)

#clean the MISTy A target_celltype names
misty_a <- misty_a %>%
  mutate(
    target_celltype_clean = target_celltype %>%
      str_remove("meanscell_abundance_w_sf_") %>%  
      gsub("\\.", " ", .)  
  )

misty_a <- misty_a %>%
  mutate(
    target_celltype_clean = case_when(
      target_celltype_clean == "MME  FAPs" ~ "MME+ FAPs",
      target_celltype_clean == "Pre Adipogenic" ~ "Pre-Adipogenic",
      TRUE ~ target_celltype_clean
    )
  )


#adding the top 3 receptor cell types to the master table

#filter avg_exp to only receptors, then create rankings
receptor_rankings <- data.frame(avg_exp) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% receptors) %>%   
  pivot_longer(-gene, names_to = "celltype", values_to = "expression") %>%
  group_by(gene) %>%
  arrange(desc(expression)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  rename(receptor = gene)

#extract top 3 cell types for each receptor
receptor_top3 <- receptor_rankings %>%
  filter(rank <= 3) %>%
  select(receptor, celltype, expression, rank) %>%
  pivot_wider(
    names_from = rank,
    values_from = c(celltype, expression),
    names_glue = "{.value}_{rank}"
  ) %>%
  rename(
    receptor_celltype_1st = celltype_1,
    receptor_expr_1st = expression_1,
    receptor_celltype_2nd = celltype_2,
    receptor_expr_2nd = expression_2,
    receptor_celltype_3rd = celltype_3,
    receptor_expr_3rd = expression_3
  )

#join with master table
master_expanded <- master_complete %>%
  left_join(receptor_top3, by = "receptor") %>%
  mutate(
    receptor_celltype_rank = case_when(
      receptor_celltype == receptor_celltype_1st ~ 1,
      receptor_celltype == receptor_celltype_2nd ~ 2,
      receptor_celltype == receptor_celltype_3rd ~ 3,
      TRUE ~ NA_real_
    )
  )


#re-arranging the columns
master_expanded <- master_expanded[ , c("ligand", "receptor", "source_niche","target_niche","mean_score","mean_magnitude",
                                        "mean_specificity","n_donors", "ligand_celltype","receptor_celltype","receptor_celltype_1st",
                                        "receptor_celltype_2nd","receptor_celltype_3rd","ligand_celltype_pct","receptor_celltype_pct",
                                        "ligand_present_in_source","receptor_present_in_target","validation_status","high_confidence",            
                                        "mean_importance", "effect_direction", "net_donors","n_donors_positive","n_donors_negative","donor_samples",
                                        "mean_morans_bivariate","max_morans_bivariate","n_donors_bivariate","donor_morans_details")]


#adding spatial validation for these new columns

#get cell type percentage in a niche
get_celltype_pct <- function(niche, celltype, comp_df) {
  celltype_col <- celltype
  if (niche %in% rownames(comp_df) && celltype_col %in% colnames(comp_df)) {
    return(comp_df[niche, celltype_col])
  }
  return(0)
}

#add percentages and presence for cell types
master_expanded <- master_expanded %>%
  rowwise() %>%
  mutate(
    #2nd cell type validation
    receptor_celltype_2nd_pct = ifelse(
      !is.na(receptor_celltype_2nd),
      get_celltype_pct(target_niche, receptor_celltype_2nd, niche_comp),
      NA_real_
    ),
    receptor_celltype_2nd_present = ifelse(
      !is.na(receptor_celltype_2nd_pct),
      receptor_celltype_2nd_pct > 5,
      NA
    ),
    
    #3rd cell type validation
    receptor_celltype_3rd_pct = ifelse(
      !is.na(receptor_celltype_3rd),
      get_celltype_pct(target_niche, receptor_celltype_3rd, niche_comp),
      NA_real_
    ),
    receptor_celltype_3rd_present = ifelse(
      !is.na(receptor_celltype_3rd_pct),
      receptor_celltype_3rd_pct > 5,
      NA
    )
  ) %>%
  ungroup()

#add a comprehensive validation status
master_expanded <- master_expanded %>%
  mutate(
    validation_status_primary = validation_status,
    
    #is 2nd cell type present?
    validation_status_secondary = case_when(
      is.na(receptor_celltype_2nd) ~ NA_character_,
      receptor_celltype_2nd_pct > 10 ~ "High confidence (2nd)",
      receptor_celltype_2nd_pct > 5 ~ "Medium confidence (2nd)",
      receptor_celltype_2nd_pct > 1 ~ "Low abundance (2nd)",
      TRUE ~ "Absent (2nd)"
    ),
    
    #is 3rd cell type present?
    validation_status_tertiary = case_when(
      is.na(receptor_celltype_3rd) ~ NA_character_,
      receptor_celltype_3rd_pct > 10 ~ "High confidence (3rd)",
      receptor_celltype_3rd_pct > 5 ~ "Medium confidence (3rd)",
      receptor_celltype_3rd_pct > 1 ~ "Low abundance (3rd)",
      TRUE ~ "Absent (3rd)"
    )
  )


##to answer- is any of the top 3 receptor cell types sufficiently present in the target niche?

master_expanded <- master_expanded %>%
  mutate(
    #best receptor pct across top 3 cell types
    best_receptor_pct = pmax(
      receptor_celltype_pct,         # 1st
      receptor_celltype_2nd_pct,     # 2nd
      receptor_celltype_3rd_pct,     # 3rd
      na.rm = TRUE
    ),
    
    #which rank is driving the best receptor evidence
    receptor_evidence_rank = case_when(
      receptor_celltype_pct >= receptor_celltype_2nd_pct & 
        receptor_celltype_pct >= receptor_celltype_3rd_pct ~ "1st",
      receptor_celltype_2nd_pct >= receptor_celltype_3rd_pct ~ "2nd",
      TRUE ~ "3rd"
    ),
    
    validation_status_extended = case_when(
      ligand_celltype_pct > 10 & receptor_celltype_pct > 10 
      ~ "High confidence (primary receiver)",
      
      ligand_celltype_pct > 10 & receptor_celltype_2nd_pct > 10 
      ~ "High confidence (secondary receiver)",
      
      ligand_celltype_pct > 10 & receptor_celltype_3rd_pct > 10 
      ~ "High confidence (tertiary receiver)",
      
      # Medium confidence
      ligand_celltype_pct > 5 & receptor_celltype_pct > 5 
      ~ "Medium confidence (primary receiver)",
      
      ligand_celltype_pct > 5 & receptor_celltype_2nd_pct > 5 
      ~ "Medium confidence (secondary receiver)",
      
      ligand_celltype_pct > 5 & receptor_celltype_3rd_pct > 5 
      ~ "Medium confidence (tertiary receiver)",
      
      # Low abundance — at least something detectable
      ligand_celltype_pct > 1 | best_receptor_pct > 1 
      ~ "Low abundance sender/receiver",
      
      TRUE ~ "Check expression"
    )
  ) %>%
  select(-best_receptor_pct, -receptor_evidence_rank)

##need to add misty A effect direction properly now

#join for 2nd receptor cell type
misty_a_2nd <- misty_a %>%
  select(ligand, target_celltype_clean, mean_importance, 
         effect_direction, net_donors, n_donors_positive, n_donors_negative) %>%
  rename(
    mean_importance_2nd = mean_importance,
    effect_direction_2nd = effect_direction,
    net_donors_2nd = net_donors,
    n_donors_positive_2nd = n_donors_positive,
    n_donors_negative_2nd = n_donors_negative
  )

#join for 3rd receptor cell type  
misty_a_3rd <- misty_a %>%
  select(ligand, target_celltype_clean, mean_importance,
         effect_direction, net_donors, n_donors_positive, n_donors_negative) %>%
  rename(
    mean_importance_3rd = mean_importance,
    effect_direction_3rd = effect_direction,
    net_donors_3rd = net_donors,
    n_donors_positive_3rd = n_donors_positive,
    n_donors_negative_3rd = n_donors_negative
  )

#add to final table
master_expanded <- master_expanded %>%
  left_join(
    misty_a_2nd,
    by = c("ligand" = "ligand", "receptor_celltype_2nd" = "target_celltype_clean")
  ) %>%
  left_join(
    misty_a_3rd,
    by = c("ligand" = "ligand", "receptor_celltype_3rd" = "target_celltype_clean")
  )

master_expanded <- master_expanded %>%
  mutate(
    best_effect_direction = case_when(
      
      #1st tier has clear directional evidence
      !is.na(effect_direction) & 
        effect_direction != "Mixed/Unclear" ~ effect_direction,
      
      #1st tier is Mixed/Unclear but 2nd tier has clear evidence
      (is.na(effect_direction) | effect_direction == "Mixed/Unclear") &
        !is.na(effect_direction_2nd) & 
        effect_direction_2nd != "Mixed/Unclear" ~ effect_direction_2nd,
      
      #1st and 2nd tier unclear/absent but 3rd tier has clear evidence
      (is.na(effect_direction) | effect_direction == "Mixed/Unclear") &
        (is.na(effect_direction_2nd) | effect_direction_2nd == "Mixed/Unclear") &
        !is.na(effect_direction_3rd) & 
        effect_direction_3rd != "Mixed/Unclear" ~ effect_direction_3rd,
      
      #All tiers are Mixed/Unclear (not NA. evidence exists but inconsistent)
      effect_direction == "Mixed/Unclear" | 
        effect_direction_2nd == "Mixed/Unclear" | 
        effect_direction_3rd == "Mixed/Unclear" ~ "Mixed/Unclear",
      
      TRUE ~ NA_character_
    ),
    
    #which tier is driving the best_effect_direction
    best_effect_tier = case_when(
      !is.na(effect_direction) & 
        effect_direction != "Mixed/Unclear" ~ "1st",
      
      (is.na(effect_direction) | effect_direction == "Mixed/Unclear") &
        !is.na(effect_direction_2nd) & 
        effect_direction_2nd != "Mixed/Unclear" ~ "2nd",
      
      (is.na(effect_direction) | effect_direction == "Mixed/Unclear") &
        (is.na(effect_direction_2nd) | effect_direction_2nd == "Mixed/Unclear") &
        !is.na(effect_direction_3rd) & 
        effect_direction_3rd != "Mixed/Unclear" ~ "3rd",
      
      TRUE ~ NA_character_
    ),
    
    best_effect_summary = case_when(
      !is.na(best_effect_direction) & !is.na(best_effect_tier) ~
        paste0(best_effect_direction, " (", best_effect_tier, " receiver)"),
      !is.na(best_effect_direction) ~ best_effect_direction,
      TRUE ~ NA_character_
    )
  )

# Save
write.csv(master_expanded, "master_expanded_with_rankings.csv", row.names = FALSE)
