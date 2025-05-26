#### Jinxin Meng, mengjx855@163.com, 20220529, 20250418, v0.6.4 ####
# 20221001: 添加选择不同距离尺度的参数 dis_method
# 20230101: update function
# 20231204: update function: check_file_name was deprecated.
# 20250107: add parameter add_lab_to_plot, show_legend, lab_size, show_grid in plot_PCoA function
# 20250224: update function
# 20250417: plot function pass to plot_dim()
# 20250418: add plot_PCoA_box().

library(dplyr)
library(tibble)
library(tidyr)
library(rstatix)
library(ggplot2)
library(vegan)
source("/code/R_func/utilities.R")
source("/code/R_func/plot_dim.R")

#### calcu_PCoA ####
calcu_PCoA <- function(profile, distance, group, group_rename = NULL, 
                       dim = 2, cumulative_eig = NULL, dis_method = "bray", 
                       prefix = NULL, adonis2 = FALSE, ...) {
  if (isTRUE(adonis2) & missing(group)) 
    stop("Error: Need group infomation when perform adonis2 analysis")
  
  if (!missing('profile') & missing('distance'))
    distance <- t(data.frame(profile, check.names = F)) %>% 
      vegdist(method = dis_method, na.rm = T)
  
  if (is.null(cumulative_eig)) {
    PCoA <- cmdscale(distance, k = dim, eig = T)
    PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
    PCoA_points <- rownames_to_column(data.frame(PCoA$points), 'sample') 
  }
  
  if (!is.null(cumulative_eig)) {
    .dim <- ncol(profile) - 1
    PCoA <- cmdscale(distance, k = .dim, eig = T) %>% suppressWarnings()
    PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
    dim <- which(cumsum(PCoA_eig) >= cumulative_eig)[1]
    PCoA_points <- dplyr::select(data.frame(PCoA$points), all_of(seq_len(dim))) %>% 
      rownames_to_column('sample')
  }
  
  if (!is.null(prefix)) 
    colnames(PCoA_points)[2:(dim + 1)] <- paste0(prefix, seq_len(dim)) 
  
  if (is.null(prefix)) 
    colnames(PCoA_points)[2:(dim + 1)] <- paste0('PC', seq_len(dim))
  
  data <- list()
  data[['point']] <- PCoA_points
  data[['dist']] <- distance
  data[['eigenvalue']] <- PCoA_eig
  data[['eig_lab']] <- paste0(colnames(PCoA_points)[2:(dim + 1)],"(", 
                              PCoA_eig[1:dim], "%)")
  
  if (isTRUE(adonis2)) {
    if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename))
      stop('group field (sample|group)')
    
    if (!is.null(group_rename)) 
      group <- data.frame(group, check.names = F) %>% 
        dplyr::rename(all_of(group_rename))
    
    group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
      as.data.frame(row.names = NULL)
    
    adonis <- adonis2(distance ~ group, group, permutations = 999)
    
    label <- paste0("'R'^2~'='~'", round(adonis[1,3], digits = 4),
                    "'~~italic('p')~'<'~'", adonis[1,5], "'") %>% 
      as.formula() %>% 
      eval()
    
    calcu_adjusted_r2 <- function(adonis_object) {
      n_observations <- adonis_object$Df[3]+1
      d_freedom <- adonis_object$Df[1]
      r2 <- adonis_object$R2[1]
      adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
      return(adjusted_r2)
    } 
    r2adj <- calcu_adjusted_r2(adonis)
    
    data[['group']] <- group
    data[['adonis2']] <- adonis
    # data[["adonis2_r2"]] <- round(adonis[1,3], digits = 4)
    # data[["adonis2_r2adj"]] <- round(r2adj, digits = 4)
    # data[["adonis2_pval"]] <- adonis[1,5]
    data[["adonis2_lab"]] <- label
  }
  
  return(data)
}

#### plot_PCoA ####
plot_PCoA <- function(profile, group, distance, group_level = NULL, 
                      group_color = NULL, group_rename = NULL,
                      dis_method = "bray", display_type = 'line', 
                      conf_type = 'ellipse', ellipse_level = .75, 
                      title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, 
                      legend_title = NULL, add_group_label = FALSE, 
                      add_sample_label = FALSE, label_size = 2, 
                      show_legend = TRUE, show_grid = FALSE, show_line = TRUE, 
                      aspect_ratio = 3/4, theme = 'default', ...) {
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(group_level)) 
    group_level <- unique(group$group)
  
  if (is.null(group_color)) {
    .color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3",
                "#a6d854","#ffd92f","#e5c494","#b3b3b3")
    .len_color <- length(.color)
    .len_level <- length(group_level)
    group_color <- rep(.color, times = ceiling(.len_level/.len_color))[1:.len_level] 
  }
  
  if (!missing('profile') & missing('distance'))
    distance <- data.frame(profile, check.names = F) %>% 
      t() %>% 
      vegan::vegdist(method = dis_method, na.rm = T)
  
  PCoA <- cmdscale(distance, k = 2, eig = T)
  PCoA_points <- data.frame(PCoA$points) %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    tibble::rownames_to_column('sample')
  PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
  plot_data <- PCoA_points %>% 
    dplyr::left_join(dplyr::select(group, sample, group), by = 'sample') %>% 
    dplyr::mutate(group = factor(group, group_level))
  
  if (is.null(xlab)) 
    xlab <- paste('PCoA1 (', PCoA_eig[1], '%)', sep = '')
  
  if (is.null(ylab)) 
    ylab <- paste('PCoA2 (', PCoA_eig[2], '%)', sep = '')
  
  if (is.null(legend_title))
    legend_title <- 'Group'
    
  if (is.null(title))
    title <- paste0(stringr::str_to_sentence(string = dis_method), '-distance PCoA')
  
  # adonis
  group <- group[match(rownames(as.matrix(distance)), group$sample),] %>%
    as.data.frame(row.names = NULL)
  adonis <- adonis2(distance ~ group, group, permutations = 999)
  
  if (is.null(subtitle))
    subtitle <- substitute("PERMANOVA: " * R^2 == a ~ ", " ~ italic(p) < b, 
                           list(a = round(adonis[1,3], 4), 
                                b = adonis[1,5]))
    # subtitle <- paste0("'R'^2~'='~'", round(adonis[1,3], digits = 4), 
    #                    "'~~italic('p')~'<'~'", adonis[1,5], "'") %>% 
    # as.formula() %>% eval()
  
  p <- plot_dim(data = plot_data, group_level = group_level, group_color = group_color, 
                display_type = display_type, conf_type = conf_type, 
                ellipse_level = ellipse_level, title = title, subtitle = subtitle, 
                xlab = xlab, ylab = ylab, legend_title = legend_title, 
                add_group_label = add_group_label, add_sample_label = add_sample_label, 
                label_size = label_size, show_legend = show_legend, 
                show_grid = show_grid, show_line = show_line, 
                aspect_ratio = aspect_ratio, theme = theme)
  
  message("  ggsave(file = \"PCoA.pdf\", width = 6, height = 4.5)")
  return(p)
}

#### plot_PCoA_box ####
plot_PCoA_box <- function(profile, group, distance, group_level = NULL, 
                          group_color = NULL, group_rename = NULL,
                          dis_method = "bray", display_type = 'line', 
                          conf_type = 'ellipse', ellipse_level = .75, 
                          title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, 
                          legend_title = NULL, add_group_label = FALSE, 
                          add_sample_label = FALSE, label_size = 1.5, 
                          show_legend = TRUE, show_grid = FALSE, show_line = TRUE, 
                          theme = 'default', ...) {
  if (!all(c('sample', 'group') %in% colnames(group)) & is.null(group_rename)) 
    stop('group field (sample|group)')
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(group_level)) 
    group_level <- unique(group$group)
  
  if (is.null(group_color)) {
    .color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3",
                "#a6d854","#ffd92f","#e5c494","#b3b3b3")
    .len_color <- length(.color)
    .len_level <- length(group_level)
    group_color <- rep(.color, times = ceiling(.len_level/.len_color))[1:.len_level] 
  }
  
  if (!missing('profile') & missing('distance'))
    distance <- data.frame(profile, check.names = F) %>% 
      t() %>% 
      vegan::vegdist(method = dis_method, na.rm = T)
  
  PCoA <- cmdscale(distance, k = 2, eig = T)
  PCoA_points <- data.frame(PCoA$points) %>% 
    dplyr::rename_with(~ c('X1', 'X2')) %>% 
    tibble::rownames_to_column('sample')
  PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
  plot_data <- PCoA_points %>% 
    dplyr::left_join(dplyr::select(group, sample, group), by = 'sample') %>% 
    dplyr::mutate(group = factor(group, group_level))
  
  if (is.null(xlab)) 
    xlab <- paste('PCoA1 (', PCoA_eig[1], '%)', sep = '')
  
  if (is.null(ylab)) 
    ylab <- paste('PCoA2 (', PCoA_eig[2], '%)', sep = '')
  
  if (is.null(legend_title))
    legend_title <- 'Group'
  
  if (is.null(title))
    title <- paste0(stringr::str_to_sentence(string = dis_method), '-distance PCoA')
  
  # adonis
  group <- group[match(rownames(as.matrix(distance)), group$sample),] %>%
    as.data.frame(row.names = NULL)
  adonis <- adonis2(distance ~ group, group, permutations = 999)
  
  if (is.null(subtitle))
    subtitle <- substitute("PERMANOVA: " * R^2 == a ~ ", " ~ italic(p) < b, 
                           list(a = round(adonis[1,3], 4), 
                                b = adonis[1,5]))

  p_main <- plot_dim(data = plot_data, group_level = group_level, group_color = group_color, 
                     display_type = display_type, conf_type = conf_type, 
                     ellipse_level = ellipse_level, xlab = xlab, ylab = ylab, 
                     legend_title = legend_title, add_group_label = add_group_label, 
                     add_sample_label = add_sample_label, label_size = label_size, 
                     show_legend = FALSE, show_grid = show_grid, show_line = show_line, 
                     theme = theme) +
    scale_x_continuous(expand = c(.01, .01)) +
    scale_y_continuous(expand = c(.01, .01))
  
  x_limits <- ggplot_build(p_main)$layout$panel_params[[1]]$x.range
  y_limits <- ggplot_build(p_main)$layout$panel_params[[1]]$y.range
  
  .top_data <- dplyr::select(plot_data, group, X1)
  .top_diff <- .top_data %>% 
    rstatix::pairwise_wilcox_test(X1 ~ group) %>% 
    dplyr::mutate(comparison = paste0(group1, '-', group2)) %>% 
    dplyr::pull(p, name = comparison) %>% 
    multcompView::multcompLetters()
  .top_label <- data.frame(label = .top_diff$Letters) %>% 
    rownames_to_column('group') %>% 
    left_join(group_by(.top_data, group) %>% 
                slice_min(order_by = X1), by = 'group')

  p_top <- ggplot(.top_data, aes(X1, group, color = group)) +
    geom_boxplot(fill = NA, outlier.shape = NA, show.legend = F, width = .61) +
    geom_jitter(size = 1.2, height = .2, show.legend = F) +
    geom_text(aes(x = X1, group, label = label), .top_label, inherit.aes = F, 
              vjust = .5, size = 4) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values = group_color) +
    scale_x_continuous(expand = c(.01, .01), limits = x_limits) +
    labs(title = title, subtitle = subtitle, y = '', x = '') +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .5, color = 'black'),
          axis.ticks.length = unit(2, 'mm'), 
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 12, color = 'black'),
          axis.text.x = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(size = 12, color = 'black'),
          plot.margin = unit(c(0, 0, 0, 0), 'mm'),
          plot.subtitle = element_text(size = 12, color = 'black'),
          panel.border = element_rect(linewidth = .5, color = 'black', fill = NA),
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          legend.background = element_blank(),
          legend.text = element_text(size = 10, color = 'black'),
          legend.title = element_text(size = 10, color = 'black'))
  
  .right_data <- dplyr::select(plot_data, group, X2)
  .right_diff <- .right_data %>% 
    rstatix::pairwise_wilcox_test(X2 ~ group) %>% 
    dplyr::mutate(comparison = paste0(group1, '-', group2)) %>% 
    dplyr::pull(p, name = comparison) %>% 
    multcompView::multcompLetters()
  .right_label <- data.frame(label = .right_diff$Letters) %>% 
    rownames_to_column('group') %>% 
    left_join(group_by(.right_data, group) %>% 
                slice_max(order_by = X2), by = 'group')
  
  p_right <- ggplot(.right_data, aes(group, X2, color = group)) +
    geom_boxplot(fill = NA, outlier.shape = NA, show.legend = F, width = .61) +
    geom_jitter(size = 1.2, width = .2, show.legend = show_legend) +
    geom_text(aes(group, X2, label = label), .right_label, inherit.aes = F, 
              vjust = .5, size = 4) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values = group_color) +
    scale_y_continuous(expand = c(.01, .01), limits = y_limits) +
    labs(y = '', x = '', color = legend_title) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .5, color = 'black'),
          axis.ticks.length = unit(2, 'mm'), 
          axis.ticks.y = element_blank(),
          axis.text = element_text(size = 12, color = 'black'),
          axis.text.y = element_blank(),
          axis.line = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'mm'),
          plot.title = element_text(size = 12, color = 'black'),
          plot.subtitle = element_text(size = 12, color = 'black'),
          panel.border = element_rect(linewidth = .5, color = 'black', fill = NA),
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          legend.background = element_blank(),
          legend.text = element_text(size = 10, color = 'black'),
          legend.title = element_text(size = 10, color = 'black'))
  
  p <- ggarrange(p_top, NULL, p_main, p_right, ncol = 2, nrow = 2,
                 heights = c(1.7, 3), widths = c(3, 1.5), align = 'hv')
  
  message("  ggsave(file = \"PCoA.pdf\", width = 8, height = 6.5)")
  return(p)
}


#### calcu_pairwise_adonis ####
calcu_pairwise_adonis <- function(profile, group, group_level = NULL,
                                  group_names = NULL, group_rename = NULL,
                                  dis_method = "bray", permutations = 999,
                                  add_plab = T, ...) {
  profile <- data.frame(profile, check.names = F)
  
  if (!is.null(group_rename)) 
    group <- data.frame(group, check.names = F) %>% 
      dplyr::rename(all_of(group_rename))
  
  if (is.null(group_level)) 
    group_level <- as.character(unique(group$group))
  
  calcu_adjusted_r2 <- function(adonis_object) {
    n_observations <- adonis_object$Df[3]+1
    d_freedom <- adonis_object$Df[1]
    r2 <- adonis_object$R2[1]
    adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
    return(adjusted_r2)
  }
  
  data <- map_dfr(combn(group_level, m = 2, simplify = F), ~ {
    .group <- dplyr::filter(group, group %in% .x)
    .data <- dplyr::select(profile, any_of(.group$sample))
    .adonis <- adonis2(t(.data) ~ group, .group, permutations = permutations, 
                       distance = dis_method)
    .r2adj <- calcu_adjusted_r2(.adonis)
    data.frame(comparison = paste0(.x[1], '_vs_', .x[2]),
               r2 = .adonis$R2[1], r2adj = .r2adj, 
               pval = .adonis$`Pr(>F)`[1]) } )
  
  if (isTRUE(add_plab))
    data <- dplyr::mutate(data, 
                          plab = cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                     labels = c('***', '**', '*', 'ns')))
  
  return(data)
}


#### calcu_distance ####
calcu_distance <- function(profile, dis_method = "bray", tree, weighted = T, ...) {
  if (dis_method == "unifrac") {
    if (missing(tree)) 
      stop("need phylogentic tree if calcu unifrac-based distance, 🥰 please input a tree: 
           ape::read.tree(\"phenogenetics.tre\")")
    
    ps <- phyloseq::phyloseq(phyloseq::otu_table(profile, taxa_are_rows = T), tree)
    distance <- phyloseq::UniFrac(ps, weighted = weighted)
    
    if (isTRUE(weighted)) 
      message("calcu Weighted-unifrac-based distance")
    if (isFALSE(weighted)) 
      message("calcu Unweighted-unifrac-based distance.")
    
  } else {
    distance <- vegdist(t(data.frame(profile, check.names = F)), method = dis_method, na.rm = T)
    message(paste0("calcu ", dis_method, "-based distance"))
  }
  return(distance)
}

#### calcu_betadisper ####
calcu_betadisper <- function(profile, group, distance, group_level = NULL, group_color = NULL, 
                             group_rename = NULL, dis_method = "bray") {
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) stop("group field (sample|group)")
  if (!is.null(group_rename)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_rename))
  if (is.null(group_level)) group_level <- unique(group$group)
  if (is.null(group_color)) {
    color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3")
    group_color <- rep(color, times = ceiling(length(group_level)/length(color)))[1:length(group_level)]
  }
  
  if (!missing(profile) & missing(distance)) {
    distance <- vegdist(t(data.frame(profile, check.names = F)), method = dis_method, na.rm = T)
  } else if (missing(profile) & !missing(distance)) {
    distance <- distance
  } else {
    stop("Error in distance compute: cannot open profile or distance data.")
  }
  
  group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
    as.data.frame(row.names = NULL)
  dispersion <- betadisper(distance, group$group)
  
  centroids <- dispersion$centroids
  points <- dispersion$vectors
  distance_to_centroid = map_dfr(group_level, \(x) 
                                 apply(points[group[group$group == x, "sample"],], 1, \(y) 
                                       sqrt(sum((y - centroids[x,])^2)) ) %>% 
                                   data.frame(value = .) ) %>% 
    rownames_to_column("sample") %>% 
    left_join(group, by = "sample") %>% 
    mutate(group = factor(group, levels = group_level))
  
  plot <- ggpubr::ggboxplot(distance_to_centroid, "group", "value", fill = "group", 
                            palette = group_color, legend = "none", 
                            xlab = "", ylab = "Distance to centroid", outlier.size = 1) +
    theme(aspect.ratio = 1)
  
  out <- list(
    dispersion = dispersion,
    distance = distance,
    group = group,
    test = list(anova = anova(dispersion),
                permutest = permutest(dispersion),
                TukeyHSD = TukeyHSD(dispersion)),
    distance_to_centroid = distance_to_centroid,
    plot = plot
  )
  
  return(out)
}
