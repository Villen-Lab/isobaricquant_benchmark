processPisaData <- function(data,
                            filterMetric = c("pepTopXIntFromTotal"),
                            filterCutoffMethod = c("defined", "percentile"),
                            filterCutoffValue = NA,
                            filterPercentile = NA){
  
  # step 1: generate modified dataframe with renamed TMT columns
  # modified_rd <- 
    data %<>% 
    
    # rename TMT quantification columns to reflect temperature treatment
    dplyr::rename(reference = Reference,
                  sequence = Sequence,
                  vehicle_00um_rep1 = Intensity.126,
                  vehicle_00um_rep2 = Intensity.127N,
                  vehicle_00um_rep3 = Intensity.127C,
                  vehicle_00um_rep4 = Intensity.128N,
                  vehicle_00um_rep5 = Intensity.128C,
                  staurosporine_01um_rep1 = Intensity.129N,
                  staurosporine_01um_rep2 = Intensity.129C,
                  staurosporine_01um_rep3 = Intensity.130N,
                  staurosporine_01um_rep4 = Intensity.130C,
                  staurosporine_01um_rep5 = Intensity.131N,
                  staurosporine_20um_rep1 = Intensity.131C,
                  staurosporine_20um_rep2 = Intensity.132N,
                  staurosporine_20um_rep3 = Intensity.132C,
                  staurosporine_20um_rep4 = Intensity.133N,
                  staurosporine_20um_rep5 = Intensity.133C,
                  blank = Intensity.134N)
  
  
  # step 2: sample loading normalization
  
  ### step 2a: determine correction factors
  sl_correction_factors <- 
    data %>% 
    ungroup %>% 
    gather(tmt, int, vehicle_00um_rep1:staurosporine_20um_rep4) %>% 
    group_by(quantification, tmt) %>% 
    summarise(int = sum(int)) %>% 
    mutate(target_intensity = mean(int, na.rm = T),
           correction_factor = target_intensity / int,
           corrected_intensity = int * correction_factor) %>% 
    ungroup %>% 
    select(quantification, tmt, correction_factor)
  
  ### step 2b: apply correction factors
  # modified_rd_sl <- 
  data %<>% 
    ungroup %>% 
    gather(tmt, int, vehicle_00um_rep1:staurosporine_20um_rep4) %>% 
    left_join(y = sl_correction_factors) %>% 
    mutate(intensity_sl = int * correction_factor)
  
  # step 3: determine filterint cutoffs
  if (filterCutoffMethod == "percentile"){
    
    filter_cutoffs <-
      data %>% 
      gather(filter, value, !!as.symbol(filterMetric)) %>% 
      group_by(quantification, filter) %>% 
      summarise(cutoff = round(c(quantile(value, probs = c(filterPercentile))), 2))
    
  } else if (filterCutoffMethod == "defined"){
    
    filter_cutoffs <-
      data %>% 
      gather(filter, value, !!as.symbol(filterMetric)) %>% 
      group_by(quantification, filter) %>% 
      summarise(cutoff = filterCutoffValue)
    
  }
  
  
  all_data <- data
  
  # all_df <- 
  data %<>% 
    gather(filter, value, !!as.symbol(filterMetric)) %>% 
    left_join(y = filter_cutoffs) %>% 
    filter(value >= cutoff) %>% 
    bind_rows(all_data %>% mutate(filter = "none", cutoff = 0))
    
  
    
    # step 4: consolidate
    
    ### step 4a: remove poor quants
  data %<>% 
      ungroup %>% 
      lazy_dt() %>% 
      filter(intensity_sl > 0) %>%
      group_by(quantification, filter, Run.ID, MSN.scan) %>% 
      filter(n() == 15) %>% 
      ungroup %>% 
      separate(tmt, c("treatment", "concentration", "replicate"), sep = "_") %>% 
    as_tibble()
    
    data %<>% 
      lazy_dt() %>% 
      group_by(reference, quantification, filter, concentration, replicate) %>% 
      summarise(intensity_sl = log2(sum(intensity_sl)),
                n_psm = n(),
                n_pep = n_distinct(sequence)) %>% 
      ungroup %>% 
      pivot_wider(names_from = concentration, values_from = intensity_sl) %>% 
      mutate(delta0100 = `01um` - `00um`,
             delta2000 = `20um` - `00um`) %>% 
      as_tibble()
      
    
    
    return(data)
    
}
