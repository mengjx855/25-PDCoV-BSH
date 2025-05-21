#### Jinxin Meng, 20250327, 20250518 ####
setwd('F:/project/20250313_PDCoV_BAC_Bile_Zhangyq/data_available/statistics/MBX-nonTarget/')
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx)
library(ropls)
library(clusterProfiler)
source('/code/R_func/difference_analysis.R')
source('/code/R_func/calcu_diff.R')

#### info ####
group_level <- c('Mock','PDCoV')
group_color <- structure(c('#a0a0a4','#d7b0b0'), names = group_level)

MBT_info <- read.delim('MBT_info.tsv')
group <- read.delim('group.tsv')
profile <- read.delim('MBT_profile.tsv', row.names = 1) %>% 
  select(all_of(group$sample))

#### PLS ####
source('/code/R_func/plot_PLS.R')

plot_PLS(profile, group, group_level = group_level, group_color = group_color,
         add_group_label = T,  label_size = 3, ellipse_level = .9, 
         show_legend = F, show_grid = T, show_line = F)
ggsave('PLS.pdf', width = 4, height = 3)

#### diff ####
# profile %>% 
#   rownames_to_column('name') %>% 
#   gather('sample', 'value', -name) %>% 
#   mutate(.value = log10(value)) %>% 
#   group_by(name) %>% 
#   shapiro_test(.value)

plsda <- opls(data.frame(t(log10(profile))), y = pull(group, group), orthoI = 0)
diff <- difference_analysis(log10(profile), group, comparison = c('PDCoV', 'Mock'), 
                            method = 't')

diff <- data.frame(plsda_vip = plsda@vipVn) %>% 
  rownames_to_column('name') %>% 
  left_join(diff, ., by = 'name') %>% 
  mutate(enriched = ifelse((plsda_vip > 1 & padj < 0.05) & log2FC > 0, 'PDCoV', 
                           ifelse((plsda_vip > 1 & padj < 0.05) & log2FC < 0, 
                                  'Mock', 'none'))) %>% 
  left_join(MBT_info, c('name' = 'cpd_id'))
write.xlsx(diff, 'diff.xlsx')

#### volcano ####
plot_data <- diff %>% 
  select(name, cpd_name, padj, log2FC, enriched, plsda_vip) %>% 
  mutate(.enriched = ifelse(is.na(enriched), 'none', enriched),
         .padj = -log10(padj),
         .log2FC = log2FC,
         .log2FC = ifelse(.log2FC < -5, -5, .log2FC),
         .log2FC = ifelse(.log2FC > 2.5, 2.5, .log2FC),
         .vip = ifelse(plsda_vip > 1, plsda_vip, .8))

ggscatter(plot_data, '.log2FC', '.padj', color = '.enriched', size = '.vip', 
          xlab = 'log2FoldChange', ylab = '-log10 BH-adjusted P-value', legend = 'right',
          palette = c(group_color, 'none' = 'grey78')) +
  scale_size_continuous(range = c(1, 2.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = 'longdash') +
  geom_vline(xintercept = c(0), linetype = 'longdash') +
  theme(aspect.ratio = 1)
ggsave('diff.volcano.pdf', width = 5, height = 5)

#### enrich ####
path_info <- read.delim('cpd2path_enrichment.tsv')

eKEGG <- diff %>% 
  filter(enriched != 'none') %>% 
  filter(KEGG != '') %>% 
  pull(KEGG) %>% 
  map(~ unlist(strsplit(.x, ';'))) %>% 
  flatten_chr() %>% 
  enricher(TERM2GENE = path_info, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame
write.xlsx(eKEGG, 'diff.eKEGG.xlsx')

path_name <- c('map00120', 'map00401', 'map00300', 'map00400', 'map00121', 'map00997', 
               'map00460', 'map00140', 'map00330', 'map00470', 'map00740', 'map00220')

plot_data <- eKEGG %>%
  mutate(path = map_vec(strsplit(ID, ':'), \(x) x[1])) %>% 
  filter(path %in% path_name) %>% 
  arrange(FoldEnrichment)

ggscatter(plot_data, 'ID', 'FoldEnrichment', fill = 'pvalue', 
          rotate = T, size = 'FoldEnrichment',
          shape = 21, legend = 'right', xlab = '', ylab = 'Fold Enrichment') +
  scale_fill_viridis_c(begin = .6) +
  scale_size_continuous(range = c(3, 6)) +
  theme(aspect.ratio = 2, 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'))
ggsave('diff.eKEGG.pdf', width = 8, height = 4.5)
