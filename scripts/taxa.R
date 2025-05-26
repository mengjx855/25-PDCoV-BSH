#### Jinxin Meng, 20220918, 20250416, v0.1.1 ####

# 20250206: fix some bug.
# 20250224: update function.
# 20250417: update function.

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
source("F:/code/r_func/profile_process.R")
source("F:/code/r_func/utilities.R")

#### taxa_split ####
# taxonomy: a data.frame containning two field, name|taxonomy..., also rename by taxonomy_name.
# taxa_seq = "d:s", 指定层级区域
taxa_split <- function(taxonomy, taxonomy_rename = NULL, sep = ";", 
                       taxa_seq = "d:s", na_fill = "Unknown", 
                       rm_suffix = F) {
  
  if (!all(c("name", "taxonomy") %in% colnames(taxonomy)) & is.null(taxonomy_rename)) 
    stop("taxonomy field (name|taxonomy)")
  
  if (!is.null(taxonomy_rename)) 
    taxonomy <- dplyr::rename(data.frame(taxonomy, check.names = F), all_of(taxonomy_rename))

  data <- data.frame(taxonomy, check.names = F)
  
  # 判断层级
  taxa_info <- structure(c("d", "p", "c", "o", "f", "g", "s", "t"),
                         names = c("domain", "phylum", "class", "order", 
                                   "family", "genus", "species", "strain")) 

  cut_off <- match(stringr::str_split_1(taxa_seq, ":"), taxa_info)
  cut_off <- cut_off[1]:cut_off[2]
  taxa_name <- names(taxa_info)[cut_off]
  taxa_level <- taxa_info[cut_off] %>% 
    paste0(., "__")
  
  out <- matrix(nrow = nrow(data), ncol = length(cut_off))
  
  for (i in 1:nrow(data)) {
    .taxa <- data$taxonomy[i]
    if (is.na(.taxa))
      out[i,] <- purrr::map_vec(taxa_level, \(x) ifelse(grepl("__$", x), 
                                                        paste0(x, na_fill), x))
    
    if (!is.na(x))
      out[i,] <- purrr::map_vec(stringr::str_split_1(.taxa, sep), \(x) 
                                ifelse(grepl("__$", x), 
                                       paste0(x, na_fill), x))[cut_off]
  }
  
  colnames(out) <- taxa_name
  out <- data.frame(out) %>% 
    add_column(name = data$name, .before = 1)
  
  # 去除门的子群标记，例如，p__Firmicutes_A, 修改为p__Firmicutes
  if (isTRUE(rm_suffix))
    out <- map_df(out, \(x) gsub("_\\w$", "", x) %>% gsub("_\\w ", " ", ., perl = T))
  
  return(out)
}

#### taxa_trans ####
# 转换物种表的层级关系
# profile: input a profile table
# taxonomy: input a data.frame, colnames following: name|phylum|family..., also rename by taxonomy_name parameter.
# taxonomy_name = a vector used to rename taxon table.
# to: 指定转变为哪个级别的
# top_n: select top n name, including Other.
# top_list: select top n name by specify name.
# other_name: a name to rename other name.
# out_all: output all result.
# na_fill: if some name a not corresponding front level, rename it.
# smp2grp: merge sample to group.
# group: input a data.frame, colnames following: sample|group, also rename by group_name parameter.
# method: merge method, normal by mean or median.
# reture a profile table.
taxa_trans <- function(profile, taxonomy, group, to = "family", top_n = 12, top_list = NULL, 
                       other_name = "Other", out_all = F, na_fill = "Unclassified",
                       transRA = F, smp2grp = F, method = "mean", 
                       taxonomy_rename = NULL, group_rename = NULL) {
  
  profile <- data.frame(profile, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  
  if (!all(c("name", to) %in% colnames(taxonomy)) & is.null(taxonomy_rename)) 
    stop(paste0("taxonomy field must have \"name\" and \"", to, 
                "\", and other level optional input."))
  
  if (!is.null(taxonomy_rename)) 
    taxonomy <- dplyr::rename(taxonomy, all_of(taxonomy_rename))
  
  if(isTRUE(out_all)) 
    top_n <- 0
  
  data <- profile %>% 
    mutate(taxa = taxonomy[match(rownames(.), taxonomy$name), to],
           taxa = ifelse(is.na(taxa), na_fill, taxa)) %>% 
    group_by(taxa) %>% 
    summarise_all(sum) %>% 
    ungroup() %>% 
    filter(taxa != '') %>% 
    column_to_rownames(var = "taxa") %>% 
    data.frame(check.names = F)
  
  # trans taxon
  if(top_n > 0 & is.null(top_list)) {
    .taxa_name <- data.frame(value = rowSums(data)) %>%
      arrange(desc(value)) %>% 
      head(n = top_n - 1) %>% 
      rownames(.) %>% 
      as.character()
    
    data <- rownames_to_column(data, "name") %>% 
      mutate(name = ifelse(name %in% .taxa_name, name, other_name)) %>% 
      group_by(name) %>% 
      summarise_all(sum) %>% 
      ungroup() %>% 
      data.frame(check.names = F) %>% 
      column_to_rownames("name")
  } 
  
  
  if (!is.null(top_list)) {
      data <- rownames_to_column(data, "name") %>% 
        mutate(name = ifelse(name %in% top_list, name, other_name)) %>% 
        group_by(name) %>% 
        summarise_all(sum) %>% 
        ungroup() %>% 
        data.frame(check.names = F) %>% 
        column_to_rownames('name')
  }
  
  # sample to group
  if (isTRUE(smp2grp)) {
    if (missing(group)) 
      stop("missing group files")
    
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) 
      stop("group field (sample|group)")
    
    if (!is.null(group_rename)) 
      group <- dplyr::rename(data.frame(group, check.names = F), 
                             all_of(group_rename))
    
    data <- data.frame(t(data), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
      group_by(group) %>% 
      summarise_all(method) %>% 
      ungroup() %>% 
      column_to_rownames("group") %>% 
      t() %>% 
      data.frame(check.names = F)
  }
  
  # transRA
  if(isTRUE(transRA)) {
    data <- apply(data, 2, \(x) x/sum(x)*100) %>%
      data.frame(check.names = F)
  }
  
  data <- data[rowSums(data) != 0, colSums(data) != 0]
  return(data)
}

#### plot_compos ####
# 物种组成,
# taxonomy 输入taxa.split()函数输出的表格，或者是指定新表格，第一列为name，随后是不同层级的分类信息
# display 指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# to 显示哪个水平组成, 可选phylum, class, family, genus, species
# top_n 输入显示多少个物种，不足数量的按最大数量算
# top_list 指定物种列表
# group_level/sample_level，分别之指定sample和group的顺序
# taxa_level/taxa_color 指定物种的排序和颜色
plot_compos <- function(profile, taxonomy, group, display = 'group', 
                        to = 'family', top_n = 12, top_list = NULL, 
                        group_level = NULL, sample_level = NULL, width = .75,
                        taxa_level = NULL, taxa_color = NULL, plot_title = NULL,
                        taxonomy_rename = NULL, group_rename = NULL) {
  if (missing(profile) & missing(taxonomy)) 
    stop("Error: missing taxonomy info.")
  
  if (display == "group" & missing(group)) 
    stop("Error: missing group info.")
  
  profile <- data.frame(profile, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  
  colors <- c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F",
              "#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6",
              "#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295",
              "#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6")
  
  if(display == 'group') {
    other_name <- paste0(stringr::str_to_lower(stringr::str_sub(to, 1, 1)), "__Other")
    
    data <- taxa_trans(profile, taxonomy, group, to = to, top_n = top_n,
                       top_list = top_list, other_name = other_name, smp2grp = T, 
                       transRA = T,group_rename = group_rename, 
                       taxonomy_rename = taxonomy_rename)
    
    if(is.null(taxa_level)) 
      taxa_level <- names(sort(rowSums(data), decreasing = T))
    
    if(is.null(taxa_color)) 
      taxa_color <- rep(colors, time = ceiling(nrow(data)/20))[1:nrow(data)]
    
    if(is.null(group_level)) 
      group_level <- names(sort(unlist(data[taxa_level[1],])))
    
    if(is.null(plot_title)) 
      plot_title <- paste0(stringr::str_to_sentence(to), ' level composition')
    
    fill_title <- paste0(stringr::str_to_sentence(to), ' taxa')
    
    plot_data <- rownames_to_column(data, "name") %>%  
      gather(., key = "group", value = "value", -name) %>% 
      mutate(group = factor(group, group_level), 
             name = factor(name, taxa_level))
    
    p <- ggbarplot(plot_data, 'group', 'value', fill = 'name', color = "#000000",
                   position = position_stack(), size = .4, width = width, 
                   palette = taxa_color, x.text.angle = 90, legend = 'right') + 
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 10, color = "#000000"),
            axis.title = element_text(size = 10, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 10, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 6, height = 4)")
    
  } else if(display == "sample") {
    other_name <- paste0(stringr::str_to_lower(stringr::str_sub(to, 1, 1)), "__Other")
    
    data <- taxa_trans(profile, taxonomy, group, to = to, top_n = top_n,
                      top_list = top_list, other_name = other_name, smp2grp = F, transRA = T,
                      group_rename = group_rename, taxonomy_rename = taxonomy_rename)
    
    if(is.null(taxa_level)) 
      taxa_level <- names(sort(rowSums(data), decreasing = T))
    
    if(is.null(taxa_color)) 
      taxa_color <- rep(colors, time = ceiling(nrow(data)/20))[1:nrow(data)]
    
    if(is.null(sample_level)) 
      sample_level <- names(sort(unlist(data[taxa_level[1],])))
    
    if(is.null(plot_title)) 
      plot_title <- paste0(stringr::str_to_sentence(to), ' Level')
    
    fill_title <- paste0(plot_title, " taxa")
    
    plot_data <- rownames_to_column(data, "name") %>%  
      gather(., key = "sample", value = "value", -name) %>% 
      mutate(sample = factor(sample, sample_level),
             name = factor(name, taxa_level))
    
    p <- ggbarplot(plot_data, 'sample', 'value', fill = 'name', color = "#000000",
                   position = position_stack(), size = .4, width = width, 
                   palette = taxa_color, x.text.angle = 90, legend = 'right') + 
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 10, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 10, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 8, height = 4)")
  }
  return(p)
}

#### plot_compos_multiple ####
# taxonomy输入taxa.split()函数输出的表格，或者是指定新表格，第一列为name，随后是不同层级的分类信息
# display指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# top_n输入显示多少个分类，不足数量的按最大数量算
plot_compos_multiple <- function(profile, taxonomy, group, display = 'group', 
                                 top_n = 12, top_list = NULL, group_level = NULL, 
                                 sample_level = NULL, taxonomy_rename = NULL, 
                                 taxa_color = NULL, group_rename = NULL, 
                                 width = .75, ...) {
  if(missing(profile) & missing(taxonomy)) 
    stop("missing taxonomy info.")
  
  .taxa_level <- setdiff(colnames(taxonomy), c('name', 'domain', 'kingdom'))
  
  p <- map(.taxa_level, ~ 
             plot_compos(profile, taxonomy, group, display = display,
                         top_n = top_n, top_list = top_list, to = .x,
                         group_level = group_level, sample_level = sample_level, 
                         taxonomy_rename = taxonomy_rename, taxa_color = taxa_color,
                         group_rename = group_rename, width = width) %>% 
             suppressMessages() ) %>% 
    cowplot::plot_grid(plotlist = ., nrow = 2, align = "v")
  
  cat("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 18, height = 8)")
  
  return(p)
}

#### plot_compos_manual ####
plot_compos_manual <- function(profile, group, display = "group", top_n = 12, 
                               out_all = F, group_level = NULL, width = .75,
                               sample_level = NULL, taxa_level = NULL, 
                               taxa_color = NULL, plot_title = NULL, 
                               fill_title = NULL, group_rename = NULL) {
  
  if(display == "group" & missing(group))
    stop("missing group info.")
  
  if (isTRUE(out_all)) 
    top_n <- 9999
  
  # profile <- apply(profile, 2, \(x) ifelse(is.nan(x), 0, x)) %>% 
  #   apply(., 2, \(x) ifelse(is.na(x), 0, x))

  profile <- profile[rowSums(profile) != 0, colSums(profile) != 0] %>% 
    data.frame(check.names = F)
  
  # colors <- c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F",
  #             "#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6",
  #             "#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295",
  #             "#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6")
  colors <- c("#80b1d3","#b3de69","#fdb462","#8dd3c7","#bc80bd","#fb8072",
              "#ffed6f","#fccde5","#bebada","#ccebc5","#ffffb3","#d9d9d9")
  
  if(display == "group") {
    data <- profile_smp2grp(profile, group, group_rename = group_rename) %>% 
      profile_transRA()
    
    .taxa_name <- names(sort(rowSums(data), decreasing = T))[1:(top_n - 1)]
    
    data <- rownames_to_column(data, "name") %>% 
      mutate(name = ifelse(name %in% .taxa_name, name, "Other")) %>% 
      group_by(name) %>% 
      summarise_all(sum) %>% 
      column_to_rownames("name")
    
    if(is.null(taxa_level)) 
      taxa_level <- names(sort(rowSums(data), decreasing = T))
    
    if(is.null(taxa_color)) 
      taxa_color <- rep(colors, time = ceiling(nrow(data)/20))[1:nrow(data)]
    
    if(is.null(group_level)) 
      group_level <- names(sort(unlist(data[taxa_level[1],])))
    
    if(is.null(plot_title)) 
      plot_title <- "Composition analysis"
    
    plot_data <- rownames_to_column(data, "name") %>%
      gather(., key = "group", value = "value", -name) %>% 
      mutate(group = factor(group, group_level), 
             name = factor(name, taxa_level))
    
    p <- ggbarplot(plot_data, 'group', 'value', fill = 'name', color = "#000000",
                   position = position_stack(), linewidth = .1, width = width, 
                   palette = taxa_color, x.text.angle = 90, legend = 'right') + 
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 10, color = "#000000"),
            axis.title = element_text(size = 10, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 10, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())

  } else if (display == "sample") {
    data <- profile_transRA(profile)
    
    .taxa_name <- names(sort(rowSums(data), decreasing = T))[1:(top_n - 1)]
    
    data <- rownames_to_column(data, "name") %>% 
      mutate(name = ifelse(name %in% .taxa_name, name, "Other")) %>% 
      group_by(name) %>% 
      summarise_all(sum) %>% 
      column_to_rownames("name")
    
    if(is.null(taxa_level))
      taxa_level <- names(sort(rowSums(data), decreasing = T))
    
    if(is.null(taxa_color)) 
      taxa_color <- rep(colors, time = ceiling(nrow(data)/20))[1:nrow(data)]
    
    if(is.null(sample_level)) 
      sample_level <- names(sort(unlist(data[taxa_level[1],])))
    
    if(is.null(plot_title)) 
      plot_title <- "Composition analysis"
    
    fill_title <- paste0(plot_title, " taxa")
    
    plot_data <- rownames_to_column(data, 'name') %>% 
      gather(key = "sample", value = "value", -name) %>% 
      mutate(sample = factor(sample, sample_level), 
             name = factor(name, taxa_level))
    
    p <- ggbarplot(plot_data, 'sample', 'value', fill = 'name', color = "#000000",
                   position = position_stack(), linewidth = .1, width = width, 
                   palette = taxa_color, x.text.angle = 90, legend = 'right') + 
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 10, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 10, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
  }
  message("  ggsave(file = \"compos_barplot.pdf\", width = 8, height = 4)")
  return(p)
}

#### plot_taxa_boxplot ####
plot_taxa_boxplot <- function(profile, group, group_level = NULL, group_color = NULL, 
                              group_rename = NULL, method = "wilcox", xlab = "", 
                              ylab = "Relative Abundance(%)", legend_title = "group", 
                              show_legend = F, aspect_ratio = 1, x_text_angle = 90) {
  profile <- data.frame(profile, check.names = F)
  feats_name <- rownames(profile)
  
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_rename)) 
    stop("group field (sample|group)")
  
  if (!is.null(group_rename)) 
    group <- dplyr::rename(data.frame(group, check.names = F), all_of(group_rename))
  
  data <- data.frame(t(profile), check.names = F) %>% 
    rownames_to_column("sample") %>% 
    dplyr::filter(sample %in% group$sample) %>% 
    left_join(group, by = "sample")
  
  if (is.null(group_level)) 
    group_level <- unique(data$group)

  color <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#ffff99",
             "#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
  
  if (is.null(group_color))
    group_color <- rep(color, times = ceiling(length(group_level)/12))[1:length(group_level)]

  p_list <- map(feats_name, \(x) {
    .data <- data %>% 
      dplyr::select(sample, group, all_of(x)) %>% 
      dplyr::rename(value = all_of(x)) %>% 
      mutate(group = factor(group, group_level))
    
    source("F:/code/R_func/calcu_diff.R")
    comparisons <- calcu_diff(dplyr::select(.data, -group),
                              dplyr::select(.data, -value), 
                              group_level = group_level, method = method) %>% 
      filter(pval < 0.05) %>% 
      pull(comparison) %>% 
      strsplit("_vs_")
    
    p <- ggboxplot(.data, 'group', 'value', fill = 'group', size = .6, 
                   outlier.shape = NA,palette = group_color, 
                   x.text.angle = x_text_angle, legend = 'right') +
      geom_jitter(size = .7, width = .4) +
      labs(x = xlab, y = ylab, color = legend_title, title = x) +
      ggpubr::stat_compare_means(comparisons = comparisons, method = method, size = 3,
                                 method.args = list(exact = F), label = "p.signif", 
                                 tip.length = .01, step.increase = .03, vjust = .9) +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 10, color = "#000000"),
            axis.title = element_text(size = 10, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000", hjust = .5, face = "italic"),
            legend.text = element_text(size = 10, color = "#000000"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank(),
            aspect.ratio = aspect_ratio) +
      guides(color = "none")
    
    if (isFALSE(show_legend))
      p <- p + guides(fill = "none")
    
    return(p)
  } )
  message()
  return(p_list)
}
