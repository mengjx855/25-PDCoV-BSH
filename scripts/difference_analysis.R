#### Jinxin Meng, 20231101, 20250501 v1.3.4 ####

# 20231025: add xlim parameter in plot_roc.
# 20231025: add trans parameter in get_feature_diff.
# 20231025: add map_name parameter in the get_feature_diff.
# 20231025: add names in the plot_volcano.
# 20231101: add prevalence statistics in output of get_feature_diff, min_abundance as a threshold of presence.
# 20231102: modify group_color in plot_volcano using function structure
# 20240401: remove plot_volcano, and rename get_feature_diff as different_analysis.
# 20250115: modify the style of sys.times() using format().
# 20250211: replace 'group_pair' to 'comparison'
# 20250327: fix some bug.
# 20250405: update some parameter.
# 20250501: update some parameter.

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
source("F:/code/R_func/utilities.R")
source("F:/code/R_func/profile_process.R")

#### difference_analysis ####
# group表中分组列必须为group，只支持两组之间比较
# comparison传入一个向量；例如c("A","B")，之后的分析即为A_vs_B
# group[data.frame]: metadata contain sample and group column, also specify by map_names parameter.
# trans[log]: used in difference analysis to transfrom profile, only including LOG in current version. [NULL, RA, LOG]
# min_abundace[num]: the threshold of a feature presencing when evaluate prevalence rate (%).
# 计算丰度这些按照原始表，计算差异考虑进行转换
difference_analysis <- function(profile, group, group_rename = NULL, prefix = NULL,
                                comparison = NULL, trans = NULL, progress = TRUE,
                                method = 'wilcox', exact = NULL, var_equal = FALSE, 
                                min_abundance = 0, add_enriched = FALSE, 
                                log2fc = 0, pvalue = 0.05, ... ) {
  
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))

  if (is.null(comparison)) 
    comparison <- unique(as.character(group$group))

  if (!all(as.character(group$group) %in% comparison)) {
    group <- dplyr::filter(group, group %in% comparison)
    profile <- dplyr::select(profile, all_of(group$sample))
  }
  
  profile <- data.frame(profile, check.names = F) %>% 
    filter(rowSums(.) != 0)
  
  .stat <- table(as.character(group$group))
  
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), "] Sample info: ", 
                 purrr::map2_vec(names(.stat), .stat, ~ paste0(.x, " (n=", .y, ").")) %>% 
                   paste0(collapse = ', ')))

  # 数据转换
  if (is.null(trans)) {
    
    data <- profile %>% 
      t() %>% 
      data.frame(check.names = F) 
    
  } else if (trans == 'RA') {
    
    data <- profile_transRA(profile) %>% 
      t() %>% 
      data.frame(check.names = F)
    
  } else if (trans == 'LOG10') {
    
    data <- profile_transLOG10(profile) %>% 
      t() %>% 
      data.frame(check.names = F)
    
  } else if (trans == 'LOG2') {
    
    data <- profile_transLOG2(profile) %>% 
      t() %>% 
      data.frame(check.names = F)
    
  } else if (trans == 'SQER') {
    
    data <- profile_transSqrt(profile) %>% 
      t() %>% 
      data.frame(check.names = F)
    
  }
  
  # 计算P值可以选择转化的数据，以适用参数检验和非参数检验
  data <- dplyr::mutate(data, group = group$group[match(rownames(data), group$sample)])
  
  # Calculate P-value.
  # pb <- txtProgressBar(style = 3, width = 50, char = "#")
  # setTxtProgressBar(pb, i/length(name)) 
  # close(pb)
  
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), 
                 "] Calculate P-value using ", method, '.'))
  
  .names <- rownames(profile)
  
  if (method == 'wilcox')
    diff <- map_dfr(.names, ~ {
      .data <- dplyr::select(data, value = all_of(.x), group)
      .wilcox <- stats::wilcox.test(value ~ group, .data, exact = exact) %>% 
        suppressWarnings()
      data.frame(name = .x,
                 pval = .wilcox$p.value,
                 method = "Wilcoxon rank-sum test") }, .progress = progress)
  
  if (method == 't') {
    
    if (isTRUE(var_equal)) 
      method_name = "student's t-test"
    
    if (isFALSE(var_equal)) 
      method_name = "Welch t-test"
    
    diff <- map_dfr(.names, ~ {
      .data <- dplyr::select(data, value = all_of(.x), group)
      .t <- stats::t.test(value ~ group, .data, var.equal = var_equal) %>% 
        suppressWarnings()
      data.frame(name = .x,
                 pval = .t$p.value,
                 method = method_name) }, .progress = progress)
  }

  diff <- dplyr::mutate(diff, padj = p.adjust(pval, method = "BH"), .after = 'pval')
  
  # Calculate feature abundance and fold-change
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), "] Calculate abundance."))
  
  # 计算丰度和变化倍数使用原始数据
  data <- data.frame(t(profile), check.names = F) %>% 
    dplyr::mutate(group = group$group[match(rownames(.), group$sample)])

  abundance <- aggregate(. ~ group, data, mean) %>% 
    tibble::column_to_rownames('group') %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    dplyr::select(X1_ab = all_of(comparison[1]), 
                  X2_ab = all_of(comparison[2])) %>% 
    dplyr::mutate(comparison = paste0(comparison, collapse = '_vs_'),
                  FC = X1_ab / X2_ab, 
                  log2FC = log2(FC + 0.001)) %>% 
    tibble::rownames_to_column('name')
  
  # Calculate feature prevalence
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), "] Calculate prevalence."))
  
  prevalence <- data %>% 
    dplyr::group_by(group) %>% 
    dplyr::group_modify(~ purrr::map_df(.x, \(x) sum(x > min_abundance)/length(x) * 100)) %>% 
    tibble::column_to_rownames('group') %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    dplyr::select(X1_pvl = all_of(comparison[1]), 
                  X2_pvl = all_of(comparison[2])) %>% 
    tibble::rownames_to_column("name")
  
  # reshape output
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), "] Output result."))
  
  out <- merge(abundance, prevalence, by = 'name', all = T) %>% 
    dplyr::left_join(diff, by = 'name') %>% 
    dplyr::relocate(X1_pvl, .after = X1_ab) %>% 
    dplyr::relocate(X2_pvl, .after = X2_ab) %>% 
    dplyr::filter(!is.nan(FC)) %>% 
    dplyr::mutate(FC = replace(FC, FC == Inf, max(FC[!is.infinite(FC)])),
                  log2FC = replace(log2FC, log2FC == Inf, 
                                   max(log2FC[!is.infinite(log2FC)])))
  
  # 替换log2FC极大值和极小值
  # log2FC = replace(log2FC, log2FC == -Inf, min(log2FC[!is.infinite(log2FC)]))
  
  if(isTRUE(add_enriched)) {
    message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), 
                   "] Default setting: log2FC = ", log2fc, "."))
    message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), 
                   "] Default setting: pvalue = ", pvalue, "."))
    out$enriched <- ifelse(out$log2FC > log2fc & out$pval < pvalue, comparison[1], 
                           ifelse(out$log2FC < -log2fc & out$pval < pvalue, 
                                  comparison[2], "none"))
  }
  
  # rename the output
  if (is.null(prefix)) 
    .renames <- structure(c('X1_ab', 'X2_ab', 'X1_pvl', 'X2_pvl'),
                          names = c(paste0(comparison, '_ab'), 
                                    paste0(comparison, '_pvl')))
  
  if (!is.null(prefix) & length(prefix) == 2)
    .renames <- structure(c('X1_ab', 'X2_ab', 'X1_pvl', 'X2_pvl'),
                          names = c(paste0(prefix, '_ab'), 
                                    paste0(prefix, '_pvl')))
  
  if (!is.null(prefix) & length(prefix) != 2)
    .renames <- structure(c('X1_ab', 'X2_ab', 'X1_pvl', 'X2_pvl'),
                          names = c(paste0(c('group1','group2'), '_ab'), 
                                    paste0(c('group1','group2'), '_pvl')))

  out <- dplyr::rename(out, all_of(.renames))
  
  message(paste0("[", format(Sys.time(), "%Y-%m-%d %X"), "] Program end."))
  return(out)
}
