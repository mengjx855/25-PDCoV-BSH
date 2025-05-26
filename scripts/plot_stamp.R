#### Jinxin Meng, 20250308, 20250327 v.1 ####

# 20250327: update some parameter.
# 20250327: add mean abundance of feature as output in calcu_stamp() with option 'method=wilcox'.
# 20250517: rename file 'calcu_stamp' as 'plot_stamp'.

library(tidyr)
library(dplyr)
library(tibble)
library(rstatix)
library(ggplot2)
library(patchwork)

#### calcu_stamp ####
calcu_stamp <- function(profile, group, group_rename = NULL, comparison = NULL,
                        pval = 0.05, method = 'wilcox', var_equal = FALSE, exact = NULL,
                        add_enriched = TRUE) {
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(comparison)) 
    comparison <- unique(group$group)
  
  if (!is.null(comparison))
    group <- data.frame(group, check.names = F) %>% 
      dplyr::filter(group %in% comparison)
  
  profile <- data.frame(profile, check.names = F) %>% 
    dplyr::select(any_of(group$sample))
  
  group <- mutate(group, group = factor(group, comparison))
  
  data <- data.frame(t(profile), check.names = F) %>% 
    rownames_to_column('sample') %>% 
    gather(key = 'name', value = 'value', -sample) %>% 
    left_join(group, by = 'sample') %>% 
    group_by(name)
  
  if (method == 't') {
    .rename_vars <- structure(c('estimate1', 'estimate2'), names = comparison)

    diff <- rstatix::t_test(data, value ~ group, var.equal = var_equal,
                            detailed = T) %>%
      dplyr::rename(.rename_vars, pval = 'p') %>%
      rstatix::add_significance(p.col = 'pval', output.col = 'plab',
                                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                symbols = c('***', '**', '*', 'ns')) %>%
      mutate(padj = p.adjust(pval, method = 'BH')) %>%
      dplyr::relocate(padj, plab, .after = pval) %>%
      dplyr::select(-.y.) %>%
      dplyr::rename_with(~ gsub('\\.', '_', .x))

    if (isTRUE(var_equal))
      diff$method <- 'Two Sample t-test'

    if (isFALSE(var_equal))
      diff$method <- 'Welch Two Sample t-test'
  }
  
  if (method == 'wilcox') {
    .ab <- profile %>%
      t %>%
      data.frame(check.names = F) %>%
      mutate(group = group$group[match(rownames(.), group$sample)]) %>%
      aggregate(. ~ group, ., mean) %>%
      column_to_rownames('group') %>%
      t %>%
      data.frame() %>%
      rownames_to_column('name')

    diff <- rstatix::wilcox_test(data, value ~ group, exact = exact,
                                 detailed = T) %>%
      dplyr::rename(pval = 'p') %>%
      rstatix::add_significance(p.col = 'pval', output.col = 'plab',
                                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                symbols = c('***', '**', '*', 'ns')) %>%
      mutate(padj = p.adjust(pval, method = 'BH')) %>%
      dplyr::relocate(padj, plab, .after = pval) %>%
      dplyr::select(-.y.) %>%
      dplyr::rename_with(~ gsub('\\.', '_', .x)) %>%
      left_join(.ab, by = 'name') %>%
      relocate(all_of(comparison), .after = 'estimate')
  }
  
  if (isTRUE(add_enriched)) {
    diff <- diff %>%
      mutate(enriched = ifelse(estimate > 0 & pval < 0.05, comparison[1],
                               ifelse(estimate < 0 & pval < 0.05,
                                      comparison[2], 'none' )))
  }
  return(diff)
}

#### plot_stamp ####
plot_stamp <- function(data, top_n = 10, comparison = NULL, 
                       palette = c('#E69F00', '#56B4E9'),
                       left_xlab = 'Mean proportion (%)', 
                       left_title = 'Relative abundance (%)', 
                       middle_xlab = 'Difference in mean proportions (%)', 
                       middle_title = '95% confidence intervals', 
                       right_xlab = '', right_title = 'P-value') {

  if (is.null(comparison)) 
    comparison <- as.character(unlist(data[1, c('group1', 'group2')]))
  
  palette <- structure(palette, names = comparison)
  
  data <- data %>% 
    filter(pval < 0.05) %>% 
    mutate(.estimate = abs(estimate)) %>% 
    arrange(desc(.estimate)) %>% 
    head(top_n) %>% 
    arrange(desc(estimate)) %>% 
    mutate(name = factor(name, name))
  
  ab <- data %>% 
    select(name, all_of(comparison)) %>% 
    gather(key = 'group', value = 'value', -name)
  
  p1 <- ggplot(ab, aes(value, name, fill = group)) +
    scale_y_discrete(limits = levels(data$name), 
                     labels = function(x) str_wrap(x, width = 40)) +
    labs(y = '', x = left_xlab, title = left_title) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4, 'lines'),
          axis.ticks = element_line(color = 'black'),
          axis.line = element_line(color = 'black'),
          axis.title.x = element_text(color = 'black', size = 12),
          axis.text = element_text(color = 'black', size = 12),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, color = 'black', vjust = .5, hjust = .5, 
                                     margin = margin(t = 0, unit = 'cm')),
          legend.position = 'bottom',
          legend.background = element_rect(fill = 'transparent'),
          legend.direction = 'horizontal',
          legend.key.width = unit(0.7, 'cm'),
          legend.key.height = unit(0.4, 'cm'),
          plot.title = element_text(size = 13, color = 'black', hjust = 0.5, 
                                    margin = margin(0, 0, .3, 0, 'cm')))
  
  for (i in 1:(nrow(data) - 1) ) {
    p1 <- p1 + annotate('rect', ymin = i + 0.5, ymax = i + 1.5, 
                        xmin = -Inf, xmax = Inf, 
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  
  p1 <- p1 + geom_bar(stat = 'identity', position = 'dodge', width = 0.7, 
                      color = 'black', show.legend = F) +
    scale_fill_manual(values = palette)
  
  p2 <- ggplot(data, aes(estimate, name, fill = enriched)) +
    scale_y_discrete(limits = levels(diff$name)) +
    labs(y = '', x = middle_xlab, title = middle_title) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4, 'lines'),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(color = 'black'),
          axis.title.x = element_text(color = 'black', size = 12),
          axis.text = element_text(color = 'black', size = 12),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = 13, color = 'black', hjust = 0.5,
                                    margin = margin(0, 0, .3, 0, 'cm')))
  
  for (i in 1:(nrow(data) - 1) ) {
    p2 <- p2 + annotate('rect', ymin = i + 0.5, ymax = i + 1.5, 
                        xmin = -Inf, xmax = Inf, 
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
  }
  
  p2 <- p2 + geom_errorbar(aes(xmin = conf_low, xmax = conf_high), 
                           position = position_dodge(0.8), width = 0.3, 
                           linewidth = 0.3) +
    geom_point(shape = 21, size = 4) +
    geom_vline(aes(xintercept = 0), linetype = 'dashed', color = 'black', 
               linewidth = .3) +
    scale_fill_manual(values = palette)
  
  p3 <- ggplot(data, aes(1, name, fill = enriched)) +
    geom_text(aes(x = 0, y = name, label = plab), color = 'black',
              vjust = .5, fontface = 'bold', inherit.aes = F, size = 4) +
    labs(x = right_xlab, title = right_title) +
    theme_void() +
    theme(plot.title = element_text(size = 13, color = 'black', hjust = 0.5, 
                                    margin = margin(0, 0, .3, 0, 'cm')))
  p <- p1 + p2 + p3 + plot_layout(widths = c(5, 5, 1))
  return(p) 
} %>% 
  suppressMessages()

