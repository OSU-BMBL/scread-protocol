
calc_overlap_list <- function(this_disease_ids, this_control_ids) {
  
  # Options based on current settings
  this_ids <- dataset %>%
    filter(region == this_region & species == this_species) %>%
    pull(data_id)
  
  this_control_ids <- dataset %>%
    filter(region == this_region & species == this_species & condition == 'Control') %>%
    pull(data_id)
  
  this_disease_ids <- dataset %>%
    filter(region == this_region & species == this_species & condition == 'Disease') %>%
    pull(data_id)
  
  overlap_dataset_thres <- function(n) {
    half <- floor(n/2) + 1
    #half <- length(n) - 1
    return(half)
  }
  
  
  this_overlap_result <- list()
  this_deg_rank_result <- data.frame()
  # iterate cell types
  ctIndex = 1
  for (ctIndex in 1:length(CT_LIST)) {
    this_ct <- CT_LIST[ctIndex]
    this_ct_short <- CT_SHORT_LIST[ctIndex]
    #dIndex = 2
    this_deg_rank <- data.frame()
    this_deg <- list()
    #iterate disease IDs
    for (dIndex in 1:length(this_disease_ids)) {
      this_disease_id <- this_disease_ids[dIndex]
      this_control_id <- de_meta %>% 
        filter(data_id == this_disease_id & description == "Disease vs control (same region)") %>%
        pull(b_data_id)
      this_comparison_deg <- de %>%
        filter(cluster == this_ct & ct == this_ct_short & a_data_id == this_disease_id & b_data_id %in% this_control_id) %>%
        arrange(p_val_adj) %>%
        filter(p_val_adj < 0.05) %>%
        filter(case_when(
          this_direction == 'up' ~ avg_logFC > 0.5,
          this_direction == 'down' ~ avg_logFC < -0.5
        )) %>%
        head(TOP) %>% 
        rownames_to_column('rank') %>%
        mutate(
          total_comparison = length(this_disease_ids)
        ) %>%
        select(ct, gene, avg_logFC, a_data_id, b_data_id,rank, total_comparison)
      
      this_deg_rank <- rbind(this_deg_rank, this_comparison_deg)
      this_deg <- list.append(this_deg, this_comparison_deg$gene)
    }
    ### Test a gene in deg list
    
    #sapply(this_deg_result, function(x){
    #  'LINC01481' %in% x
    #})
    
    
    ########### Calculate unioned overlaps from one cell type in all datasets
    this_deg_overlap <- vector()
    half <- OVERLAP_THRES
    # Combined
    data_combines <- combn(length(this_disease_ids),half)
    
    #cIndex <- 1
    for (cIndex in 1:ncol(data_combines)) {
      this_combine <- as.vector(data_combines[,cIndex])
      venn_set <- Venn(this_deg)
      tmp_venn <- overlap(venn_set, this_combine)
      this_deg_overlap <- union(this_deg_overlap, tmp_venn)
    }
    ########### Calculate unioned overlaps from one cell type in all datasets
    this_deg_rank_result <- this_deg_rank %>%
      filter(gene %in% this_deg_overlap) %>%
      rbind(this_deg_rank_result)
    this_overlap_result <- list.append(this_overlap_result, this_deg_overlap)
    
  }
  
  names(this_overlap_result) <- CT_SHORT_LIST
  
  freq_deg <- this_deg_rank_result %>%
    group_by(ct) %>%
    count(gene) %>%
    group_by(gene)%>%
    count(gene,sort = T)
  
  #i=1
  freq_comparisons <- this_deg_rank_result %>%
    filter(gene %in% freq_deg$gene) %>%
    group_by(ct,gene) %>%
    count(gene,sort = T,name = "overlapping_comparison")
  
  overlap_list_result <- tibble()
  for (i in 1:nrow(freq_deg)) {
    this_gene <- freq_deg %>%
      pull(gene) %>%
      nth(i)
    overlap_list_result <- this_deg_rank_result %>%
      filter(gene == this_gene) %>%
      pull(ct) %>%
      unique() %>%
      sort() %>% 
      paste(collapse = ",") %>%
      tibble(ct=.,gene=this_gene) %>%
      rbind(overlap_list_result) %>%
      arrange(ct) %>%
      arrange(desc(nchar(ct))) %>%
      mutate(across(where(is_character),as_factor))
  }
  
  overlap_rank_result <- this_deg_rank_result %>%
    filter(gene %in% freq_deg$gene) %>%
    group_by(gene) %>%
    mutate(
      mean_rank = mean(as.numeric(rank))
    ) %>%
    arrange(gene) %>%
    left_join(freq_comparisons, by="gene") %>%
    tibble() %>%
    mutate(across(where(is_character),as_factor)) %>%
    select(ct.x, gene, avg_logFC, a_data_id, b_data_id, rank, overlapping_comparison, total_comparison, mean_rank) %>%
    rename(ct = ct.x, disease_id = a_data_id, control_id = b_data_id) %>%
    filter(!duplicated(avg_logFC))
  
  return(list(list=overlap_list_result, rank=overlap_rank_result))
}
