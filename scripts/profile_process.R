#### Jinxin Meng, 20240308, 20250227 v0.3 ####

# 2023-11-01: add all_group parameter in profile_filter.
# 2023-12-19: add function profile_replace.
# 2025-02-27: add function profile_top_n/frac.

library(dplyr) 
library(tibble)
library(tidyr)
library(purrr)


#### LOG ####
# LOG transformation method in MaAsLin2;
# The default log transformation incorporated into MaAsLin does add a pseudo-count;
# As is best-known practice currently, the pseudo-count is half the minimum feature; 
# x [numeric]: a numeric vector.
LOG2 <- function(x) log2(replace(x, x == 0, min(x[x > 0]) / 2))
LOG10 <- function(x) log10(replace(x, x == 0, min(x[x > 0]) / 2))


#### profile_transLOG ####
# conduct LOG transformation for a otu_table, rows represent features, 
# and columns represent samples.
profile_transLOG2 <- function(profile) {
  apply(profile, 1, LOG2) %>% 
    t() %>% 
    data.frame(check.names = F)
  }
  
profile_transLOG10 <- function(profile) {
  apply(profile, 1, LOG10) %>% 
    t() %>% 
    data.frame(check.names = F)
  }

profile_transSqrt <- function(profile) {
  sqrt(profile)
  }

#### profile_transRA ####
# Jinxin Meng, 20231119
# replace some value meeting parameter to specified value.
profile_transRA <- function(profile) {
  apply(data.frame(profile, check.names = F), 2, \(x) x/sum(na.omit(x))*100) %>%
    data.frame(check.names = F)
  }
  
#### profile_filter ####
# profile: input a data.frame of relative abundance profile.
# group: mapping (sample|group), also specify by map_names parameter.
# by_group: filter feature in each group.
# all_group prevalence only meet all group are outputted.
# min_prevalence: threshold of prevalence of features in all sample.
# min_abundance: threshold of abundance in a sample is considered a feature presenting in the sample.
profile_filter <- function(profile, group, group_rename = NULL, 
                           by_group = F, all_group = F, n_group = 1, 
                           min_prevalence = 0.1, min_n = NULL,
                           min_abundance = 0.0) {
  profile <- data.frame(profile, check.names = F)
  
  if (!is.null(min_n) & is.numeric(min_n)) 
    min_prevalence <- NULL
  
  if (isTRUE(all_group)) 
    n_group <- NULL
 
  # 按照分组过滤，每个组内计算流行率或者计数，然后进行过滤
  if (isTRUE(by_group) & !missing(group)) { 
    if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
      stop('group field (sample|group)')
    
    if (!is.null(group_rename)) 
      group <- data.frame(group, check.names = F) %>% 
        dplyr::rename(all_of(group_rename))
    
    if (!is.null(min_prevalence) & is.numeric(min_prevalence)) {
      prevalence <- profile %>% 
        t() %>% 
        data.frame(check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% 
        group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))) %>%
        ungroup() %>% 
        dplyr::select(-group)
      
      if (isTRUE(all_group))
        .vec <- map_vec(prevalence, \(x) all(x > min_prevalence))
      
      if (isFALSE(all_group))
        .vec <- map_vec(prevalence, \(x) sum(x > min_prevalence) >= n_group)
      
    } else if (!is.null(min_n) & is.numeric(min_n)) {
      
      count <- profile %>% 
        t() %>% 
        data.frame(check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% 
        group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance))) %>%
        ungroup() %>% 
        dplyr::select(-group)
      
      if (isTRUE(all_group))
        .vec <- map_vec(count, \(x) all(x >= min_n))
        
      if (isFALSE(all_group))
        .vec <- map_vec(count, \(x) sum(x >= min_n) >= n_group)
    }
    
    .vec <- names(.vec[.vec])
    profile <- profile[.vec, ]
  } 
  
  # 不按照分组。计算整体的流行率
  if (!isTRUE(by_group)) {
    if (!is.null(min_prevalence) & is.numeric(min_prevalence))
      .vec <- data.frame(t(profile), check.names = F) %>% 
        purrr::map_vec(\(x) sum(x > min_abundance)/length(x) > min_prevalence)
    
    if (!is.null(min_n) & is.numeric(min_n))
      .vec <- data.frame(t(profile), check.names = F) %>% 
        purrr::map_vec(\(x) sum(x > min_abundance) >= min_n)
    
    .vec <- names(.vec[.vec])
    profile <- profile[.vec,]
    
  }
  return(profile)
}
#### profile_smp2grp ####
# Jinxin Meng, 20231118
# method: merge column method. please get help in summarise_all() function.
profile_smp2grp <- function(profile, group, group_rename = NULL, method = 'mean'){
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  profile <- data.frame(t(profile), check.names = F) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    aggregate(. ~ group, ., method) %>% 
    ungroup() %>% 
    column_to_rownames('group') %>%
    t() %>% 
    data.frame(check.names = F)
  return(profile)
}


#### profile_smp2grp_1 ####
# Jinxin Meng, 20240814, for profile containing NA or NAN .....
profile_smp2grp_1 <- function(profile, group, group_rename = NULL, method = "mean"){
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) 
    stop("group field (sample|group)")
  if (!is.null(group_rename)) group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  data <- data.frame(t(profile), check.names = F) %>% 
    mutate(group = group$group[match(rownames(.), group$sample)])
  
  if (method == "mean") out <- data %>% 
    group_by(group) %>% 
    group_modify(~purrr::map_df(.x, \(x) ifelse(sum(!is.na(x)) == 0, 0, mean(na.omit(x))))) %>% 
    ungroup()

  if (method == "median") out <- data %>% 
      group_by(group) %>% 
      group_modify(~purrr::map_df(.x, \(x) ifelse(sum(!is.na(x)) == 0, 0, median(na.omit(x))))) %>% 
      ungroup()
  out <- out %>% 
    column_to_rownames("group") %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column("name")
  return(out)
}

#### profile_top_n ####
profile_top_n <- function(profile, n = 12, out_other = F, other_name = 'Other', 
                          sort_method = 'mean') {
  .feats <- apply(profile, 1, sort_method) %>% 
    sort(decreasing = T) %>% 
    head(n = n) %>% 
    names()

  if (isFALSE(out_other)) profile <- filter(profile, rownames(profile) %in% .feats)
  
  if(isTRUE(out_other)) {
    profile <- profile %>% 
      mutate(name = rownames(.),
             name = replace(name, !name %in% .feats, other_name)) %>% 
      aggregate(. ~ name, ., sum) %>% 
      column_to_rownames('name')
    }
  
  return(profile)
  }

#### profile_top_frac ####
profile_top_frac <- function(profile, frac = .1, out_other = F,
                             other_name = 'Other', sort_method = 'mean') {
  .feats <- apply(profile, 1, sort_method) %>% 
    sort(decreasing = T) %>% 
    head(n = floor(nrow(profile) * frac)) %>% 
    names()
  
  if (isFALSE(out_other)) profile <- filter(profile, rownames(profile) %in% .feats)
  
  if(isTRUE(out_other)) {
    profile <- profile %>% 
      mutate(name = rownames(.),
             name = replace(name, !name %in% .feats, other_name)) %>% 
      aggregate(. ~ name, ., sum) %>% 
      column_to_rownames('name')
    }
  
  return(profile)
  }

#### profile_replace ####
# Jinxin Meng, 20231119
# replace some value meeting parameter to specified value.
profile_replace <- function(profile, min_value = 1, fill_value = 0, transRA = F){
  profile <- data.frame(profile, check.names = F)
  
  if (isTRUE(transRA)) 
    profile <- apply(profile, 2, x/sum(x) * 100) %>% 
      data.frame(check.names = F)
  
  profile <- apply(profile, 2, \(x) ifelse(x > min_value, x, fill_value)) %>% 
    data.frame(check.names = F)
  
  profile <- profile[rowSums(profile) != 0, colSums(profile) != 0]
  return(profile)
}

#### profile_adjacency ####
# Jinxin Meng, 20231119
# profile to adjacent matrix
profile_adjacency <- function(profile){
  profile <- apply(profile, 2, \(x) ifelse(x > 0, 1, 0)) %>% data.frame(check.names = F)
  profile <- profile[rowSums(profile) != 0, colSums(profile) != 0]
  return(profile)
}
  
#### profile_prevalence ####
# prevalence of each feature in all samples or samples belonged to each group
profile_prevalence <- function(profile, group, by_group = T, min_abundance = 0, count = F) {
  if (isTRUE(by_group) & missing(group)) stop("if by_group = TRUE, group (field: sample|group) should be provided ..")
  if (isFALSE(count)) {
    if (isTRUE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x))*100) %>%
        column_to_rownames(var = "group") %>% t %>% data.frame(check.names = F)
    } else if (isFALSE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        map_vec(\(x) sum(x > min_abundance)/length(x)*100) %>% 
        data.frame(prevalence = .) %>% rownames_to_column(var = "name")
    }
  } else {
    if (isTRUE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
        group_by(group) %>% group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance))) %>%
        column_to_rownames(var = "group") %>% t %>% data.frame(check.names = F)
    } else if (isFALSE(by_group)) {
      prevalence <- data.frame(t(profile), check.names = F) %>% 
        map_vec(\(x) sum(x > min_abundance)) %>% 
        data.frame(prevalence = .) %>% rownames_to_column(var = "name")
    }
  }
  return(prevalence)
}

#### profile_statistics ####
profile_statistics <- function(profile, group, group_rename = NULL, 
                               by_group = T) {
  profile <- data.frame(profile, check.names = F)
  
  if (isTRUE(by_group) & !missing(group)) { 
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) stop("group field (sample|group)")
    if (!is.null(group_rename)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_rename))
    data <- profile %>% 
      t %>% 
      data.frame(check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
    var = unique(group$group)[1]
    
    data <- cbind(
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) mean(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_mean")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) sd(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_sd")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) median(x))) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_median")),
      data %>% 
        group_by(group) %>% 
        group_modify(~map_df(.x, \(x) sum(x>0)/length(x)*100)) %>% 
        column_to_rownames(var = "group") %>% 
        t %>% 
        data.frame(check.names = F) %>% 
        dplyr::rename_with(~paste0(.x, "_prev"))
    ) %>% 
      relocate(starts_with(var))
  } else if (!isTRUE(by_group)) {
    
  } else {
    stop("if by_group = TRUE, group (field: sample|group) should be provided ..")
  }
  
}