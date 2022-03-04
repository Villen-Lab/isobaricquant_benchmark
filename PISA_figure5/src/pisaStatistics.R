pisaStatistics <- function(data,
                           proteinsToUse = c("all", "overlap"),
                           minPeptides = 2,
                           minPSMs = 2,
                           fcCutoff = log2(1)){

  # step 1: determine which proteins were identified by a minimum number of unique peptides
  protein_id <- 
    data %>% 
    ungroup %>% 
    filter(filter == "none") %>% 
    filter(n_pep >= minPeptides)
  
  # then, filter for proteins identified by the required number of peptides and
  # quantified by the required number of PSMs (regardless of number of peptides)
  data %<>% 
    ungroup %>%
    filter(reference %in% protein_id$reference) %>% 
    filter(n_psm >= minPSMs)
  
  # optional: if test done on just overlap, then trim down protein dataframe to proteins
  # identified and quantified in all conditions being compared (i.e. MS2 vs MS3 vs filtered)
  # otherwise, use all proteins available
  if (proteinsToUse == "overlap"){
    
    data %<>% 
      mutate(groupID = paste0(replicate, filter, quantification)) %>% 
      ungroup %>% 
      mutate(uniqueID = n_distinct(groupID)) %>% 
      group_by(reference) %>%
      filter(n() == uniqueID) %>% 
      ungroup
    
  }
  
  
  # step 3: statistical test for proteins with significant change  
  data %<>% 
    ungroup %>%
    group_by(quantification, filter, reference) %>%
    nest() %>% 
    mutate(ttest = map(.x = data, .f = tTest)) %>%
    select(-data) %>%
    unnest(ttest) %>% 
    ungroup %>% 
    group_by(quantification, filter) %>%
    mutate(padj_01 = p.adjust(pval_01, "BH"),
           padj_20 = p.adjust(pval_20, "BH")) %>% 
    pivot_longer(cols = c(delta_01, delta_20, sumSD_01, sumSD_20, pval_01, pval_20, padj_01, padj_20), 
                 names_to = "concentration",
                 values_to = "value") %>% 
    separate(concentration, c("group", "concentration"), sep = "_") %>% 
    pivot_wider(names_from = group, values_from = value)
  
  testResults <- 
    data %>% 
    mutate(kinase = ifelse(reference %in% human_kinases$reference, "kinase", "other")) %>% 
    mutate(impact = case_when(delta < -fcCutoff & padj < 0.05 ~ "destabilized",
                              delta > fcCutoff & padj < 0.05 ~ "stabilized",
                              TRUE ~ "n.s."))
  
  return(testResults)
  
}
 