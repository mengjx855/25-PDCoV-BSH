#### Jinxin Meng, 20250328, 20250415 ####
setwd("F:/project/20250313_PDCoV_BAC_Bile_Zhangyq/data_available/statistics/MBX-bileTarget/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx)
library(ComplexHeatmap)
source("/code/R_func/profile_process.R")

#### abundance boxplot ####
grp <- c('Mock','PDCoV')
grp_col <- structure(c('#a0a0a4','#d7b0b0'), names = grp)

.level <- c('CA','CDCA','HCA','DCA','LCA','UDCA','HDCA','GCA','TCA','GCDCA',
            'TCDCA','GHCA','THCA','THDCA','GUDCA','TUDCA','GDCA','TDCA',
            'GLCA','TLCA')

profile <- read.delim('BA_data.txt', sep = "\t", row.names = 1)

group <- data.frame(sample = colnames(profile)) %>% 
  mutate(group = gsub('\\d+', '', sample))

data <- profile %>% 
  rownames_to_column('name') %>% 
  gather('sample', 'value', -name) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(.value = log10(value),
         name = factor(name, .level))

data %>% 
  group_by(name) %>% 
  rstatix::anova_test(.value ~ group) %>% 
  rstatix::add_significance('p')

ggboxplot(data, 'name', '.value', fill = 'group', palette = grp_col, width = .6,
          outlier.shape = NA, xlab = '', ylab = 'log10 content', legend = 'none',
          x.text.angle = 45) +
  annotate('text', x = c(1, 4, 5, 8, 19), y = c(3, 2.5, 2.8, 2, 1), 
           label = c('**', '**', '*', '**', '*'), color = 'red') +
  annotate('segment', x = .75, xend = 1.25, y = 2.9) +
  annotate('segment', x = 3.75, xend = 4.25, y = 2.4) +
  annotate('segment', x = 4.75, xend = 5.25, y = 2.7) +
  annotate('segment', x = 7.75, xend = 8.25, y = 1.9) +
  annotate('segment', x = 18.75, xend = 19.25, y = 0.9)
ggsave('biles_ab.pdf', width = 6, height = 3)

#### abundance heatmap #### 
grp <- c('Mock','PDCoV')
grp_col <- structure(c('#a0a0a4','#d7b0b0'), names = grp)

col_data <- data.frame(name = colnames(profile)) %>% 
  left_join(group, by = c('name' = 'sample'))

col_split <- factor(col_data$group)

pdf('biles_ab.heatmap.pdf', width = 7, height = 7)
pheatmap(log10(profile), scale = 'row',
         cellwidth = 10,cellheight = 10,
         colorRampPalette(c("#3288bd", "#ffffff", "#d53e4f"))(100),
         fontsize_col = 8, fontsize_row = 8,
         treeheight_col = 0, treeheight_row = 0,
         border = "#ffffff", border_gp = gpar(col = "#000000"),
         column_split = col_split,
         column_gap = unit(1, "mm"))
dev.off()


#### accociation ####
MGX_profile <- read.delim('species_data.txt', sep = "\t", row.names = 1)
MBX_profile <- read.delim('BA_data.txt', sep = "\t", row.names = 1)

corr <- psych::corr.test(t(MGX_profile), t(MBX_profile), method = "spearman")
r_data <- data.frame(t(corr$r), check.names = F)
p_data <- data.frame(t(corr$p), check.names = F)

plot_label <- apply(p_data, 2, \(x) ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", "")))
plot_label <- plot_label[(apply(plot_label, 1, \(x) sum(x!='')) != 0),]
r_data <- r_data[rownames(plot_label),]

pheatmap((r_data),
         cluster_rows = T, cluster_cols = T,
         border = "#ffffff", border_gp = gpar(col = "#000000"),
         show_rownames = T, show_colnames = T,
         cellwidth = 10,cellheight = 10,
         treeheight_col = 0, treeheight_row = 0,
         fontsize_col = 8, fontsize_row = 8,
         display_numbers = (plot_label),
         fontface_col = "italic")
