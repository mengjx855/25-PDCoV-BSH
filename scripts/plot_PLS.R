#### Jinxin Meng, 20250328, 20250417, v0.2.1 ####
# 20250416: update function.
# 20250417: plot function pass to plot_dim()
# 20250419: add options sub_sample, sub_group for plot_PLS()

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ropls)
source("/code/R_func/plot_dim.R")

#### plot_PLS ####
plot_PLS <- function(profile, group, group_level = NULL, group_color = NULL, 
                     group_rename = NULL, display_type = 'line',
                     sub_sample = NULL, sub_group = NULL,
                     conf_type = 'ellipse', ellipse_level = .75, 
                     title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, 
                     legend_title = NULL, add_group_label = FALSE, 
                     add_sample_label = FALSE, label_size = 1.5, 
                     show_legend = TRUE, show_grid = FALSE, show_line = TRUE, 
                     aspect_ratio = 3/4, theme = 'default', ...) {
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (!is.null(sub_sample) & is.vector(sub_sample))
    group <- dplyr::filter(group, sample %in% sub_sample)
  
  if (!is.null(sub_group) & is.vector(sub_group))
    group <- dplyr::filter(group, group %in% sub_group)
  
  if (is.null(group_level)) 
    group_level <- unique(group$group)
  
  if (is.null(group_color)) {
    .color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3",
                "#a6d854","#ffd92f","#e5c494","#b3b3b3")
    .len_color <- length(.color)
    .len_level <- length(group_level)
    group_color <- rep(.color, times = ceiling(.len_level/.len_color))[1:.len_level]
  }
  
  group_color <- structure(group_color, names = group_level)
  
  profile <- data.frame(profile, check.names = F) %>% 
    select(all_of(group$sample)) %>% 
    filter(rowSums(.) != 0)
  
  pls <- opls(data.frame(t(profile)), y = as.character(pull(group, group)), 
              orthoI = 0, predI = 2, info.txtC = 'none', fig.pdfC = 'none')
  
  pls_points <- data.frame(pls@scoreMN) %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    rownames_to_column('sample')
  
  plot_data <- pls_points %>% 
    dplyr::left_join(dplyr::select(group, sample, group), by = 'sample') %>% 
    dplyr::mutate(group = factor(group, group_level))

  if (is.null(xlab)) 
    xlab <- paste0("Dim1 (", round(pls@modelDF[1, 1] * 100, 2), "%)")
  
  if (is.null(ylab)) 
    ylab <- paste0("Dim2 (", round(pls@modelDF[2, 1] * 100, 2), "%)")
  
  if (is.null(legend_title))
    legend_title <- 'Group'
  
  if(is.null(subtitle)) 
    subtitle <- substitute(R^2 * X == a ~~ R^2 * Y == b ~~ Q^2 == c ~~ RMSEE == d, 
                           list(a = round(pls@summaryDF[1, 1], 3), 
                                b = round(pls@summaryDF[1, 2], 3),
                                c = round(pls@summaryDF[1, 3], 3), 
                                d = round(pls@summaryDF[1, 4], 3)))

  if(is.null(title)) 
    title <- 'Partial least aquares discriminant analysis'
  
  p <- plot_dim(data = plot_data, group_level = group_level, group_color = group_color, 
                display_type = display_type, conf_type = conf_type, 
                ellipse_level = ellipse_level, title = title, subtitle = subtitle, 
                xlab = xlab, ylab = ylab, legend_title = legend_title, 
                add_group_label = add_group_label, add_sample_label = add_sample_label, 
                label_size = label_size, show_legend = show_legend, 
                show_grid = show_grid, show_line = show_line, 
                aspect_ratio = aspect_ratio, theme = theme)

  message("  ggsave(file = \"PLS.pdf\", width = 6, height = 4.5)")
  return(p)
}

#### plot_OPLS ####
plot_OPLS <- function(profile, group, group_level = NULL, group_color = NULL, 
                     group_rename = NULL, display_type = 'line',
                     sub_sample = NULL, sub_group = NULL, top_frac = .2,
                     conf_type = 'ellipse', ellipse_level = .75, 
                     title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, 
                     legend_title = NULL, add_group_label = FALSE, 
                     add_sample_label = FALSE, label_size = 1.5, 
                     show_legend = TRUE, show_grid = FALSE, show_line = TRUE, 
                     aspect_ratio = 3/4, theme = 'default', ...) {
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (!is.null(sub_sample) & is.vector(sub_sample))
    group <- dplyr::filter(group, sample %in% sub_sample)
  
  if (!is.null(sub_group) & is.vector(sub_group))
    group <- dplyr::filter(group, group %in% sub_group)
  
  if (length(unique(as.character(group$group))) != 2)
    stop('OPLS-DA only available for binary classification (use PLS-DA for multiple classes)')
  
  if (is.null(group_level)) 
    group_level <- unique(group$group)
  
  if (is.null(group_color)) {
    .color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3",
                "#a6d854","#ffd92f","#e5c494","#b3b3b3")
    .len_color <- length(.color)
    .len_level <- length(group_level)
    group_color <- rep(.color, times = ceiling(.len_level/.len_color))[1:.len_level] 
  }
  
  group_color <- structure(group_color, names = group_level)
  
  # top_frac 选择方差最大的top比例的代谢物进行分析
  .names <- apply(profile, 1, sd) %>% 
    sort(decreasing = T) %>% 
    head(length(.) * top_frac) %>% 
    names()
  
  profile <- data.frame(profile, check.names = F) %>% 
    select(all_of(group$sample)) %>% 
    filter(rowSums(.) != 0) %>% 
    filter(rownames(.) %in% .names)
  
  opls <- opls(data.frame(t(profile)), y = as.character(pull(group, group)), 
               orthoI = NA, predI = 5, info.txtC = 'none', fig.pdfC = 'none') %>% 
    suppressWarnings()
  
  if (nrow(opls@modelDF) == 0)
    stop("No model was built because the first predictive component was already not significan")
  
  opls_points <- cbind(opls@scoreMN, opls@orthoScoreMN) %>% 
    data.frame() %>% 
    dplyr::select(p1, o1) %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    rownames_to_column('sample')
  
  plot_data <- opls_points %>% 
    dplyr::left_join(dplyr::select(group, sample, group), by = 'sample') %>% 
    dplyr::mutate(group = factor(group, group_level))
  
  if (is.null(xlab)) 
    xlab <- paste0("Dim1 (", round(opls@modelDF[1, 1] * 100, 2), "%)")
  
  if (is.null(ylab)) 
    ylab <- paste0("Dim2 (", round(opls@modelDF[2, 1] * 100, 2), "%)")
  
  if (is.null(legend_title))
    legend_title <- 'Group'
  
  if(is.null(subtitle)) 
    subtitle <- substitute(R^2 * X == a ~~ R^2 * Y == b ~~ Q^2 == c ~~ RMSEE == d, 
                           list(a = round(opls@summaryDF[1, 1], 3), 
                                b = round(opls@summaryDF[1, 2], 3),
                                c = round(opls@summaryDF[1, 3], 3), 
                                d = round(opls@summaryDF[1, 4], 3)))
  
  if(is.null(title))
    title <- 'Orthogonal partial least squares discriminant analysis'
  
  p <- plot_dim(data = plot_data, group_level = group_level, group_color = group_color, 
                display_type = display_type, conf_type = conf_type, 
                ellipse_level = ellipse_level, title = title, subtitle = subtitle, 
                xlab = xlab, ylab = ylab, legend_title = legend_title, 
                add_group_label = add_group_label, add_sample_label = add_sample_label, 
                label_size = label_size, show_legend = show_legend, 
                show_grid = show_grid, show_line = show_line, 
                aspect_ratio = aspect_ratio, theme = theme)
  
  message("  ggsave(file = \"OPLS.pdf\", width = 6, height = 4.5)")
  return(p)
}