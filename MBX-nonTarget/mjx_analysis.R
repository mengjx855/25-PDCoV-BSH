#### Jinxin Meng, 20250327, 20250415 ####
setwd("F:/project/20250313_PDCoV_BAC_Bile_Zhangyq/MBX-nonTarget/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx)
library(ropls)
library(clusterProfiler)
source("/code/R_func/difference_analysis.R")
source("/code/R_func/calcu_diff.R")

#### info ####
grp <- c('Mock','PDCoV')
grp_col <- structure(c('#a0a0a4','#d7b0b0'), names = grp)

MBT_info <- read.delim("MBT_info.tsv")
group <- read.delim("group.tsv")
profile <- read.delim("MBT_profile.tsv", row.names = 1) %>% 
  select(all_of(group$sample))

#### PLS ####
source("/code/R_func/plot_PLS.R")

plot_PLS(profile, group, group_order = grp, group_color = grp_col,
         add_group_label = T,  lab_size = 3, ellipse_level = .9, 
         show_legend = F, show_grid = T, show_line = F)
ggsave("PLS.pdf", width = 4, height = 3)

#### diff ####
plsda <- opls(data.frame(t(log10(profile))), y = pull(group, group), orthoI = 0)
diff <- difference_analysis(log10(profile), group, comparison = c('PDCoV', 'Mock'), 
                            method = 't')

diff <- data.frame(plsda_vip = plsda@vipVn) %>% 
  rownames_to_column("name") %>% 
  left_join(diff, ., by = "name") %>% 
  mutate(enriched = ifelse((plsda_vip > 1 | pval < 0.05) & log2FC > 0, "PDCoV", 
                           ifelse((plsda_vip > 1 | pval < 0.05) & log2FC < 0, 
                                  "Mock", "none"))) %>% 
  left_join(MBT_info, c('name' = 'cpd_id'))
write.xlsx(diff, "diff.xlsx")

#### volcano ####
bile_related <- read.delim("/database/KEGG/enrichment_analysis/cpd2path_enrichment.tsv") %>% 
  filter(grepl('bile', pathway)) %>% 
  pull(ID) %>% 
  unique()

plot_label <- diff %>% 
  filter(rowSums(sapply(bile_related, \(x) grepl(x, KEGG))) > 0) %>% 
  filter(enriched != 'none')
  
plot_data <- diff %>% 
  select(name, cpd_name, pval, log2FC, enriched, plsda_vip) %>% 
  mutate(.enriched = ifelse(is.na(enriched), 'none', enriched),
         .pval = -log10(pval),
         .log2FC = log2FC,
         .log2FC = ifelse(.log2FC < -5, -5, .log2FC),
         .log2FC = ifelse(.log2FC > 2.5, 2.5, .log2FC),
         .vip = ifelse(plsda_vip > 1, plsda_vip, .8),
         .label = ifelse(name %in% plot_label$name, cpd_name, ''))

ggscatter(plot_data, ".log2FC", ".pval", color = ".enriched", size = '.vip', 
          xlab = 'log2FoldChange', ylab = "-log10 P-value",
          palette = c(grp_col, 'none' = 'grey78')) +
  scale_size_continuous(range = c(1, 2.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(0), linetype = "dashed") +
  theme(aspect.ratio = 1)
ggsave('diff.volcano.pdf', width = 5, height = 5)

ggscatter(plot_data, ".log2FC", ".pval", color = ".enriched", size = '.vip', 
          xlab = 'log2FoldChange', ylab = "-log10 P-value",
          palette = c(grp_col, 'none' = 'grey78')) +
  scale_size_continuous(range = c(1, 2.5)) +
  geom_text(aes(label = .label), size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(0), linetype = "dashed") +
  theme(aspect.ratio = 1)
ggsave('diff.volcano.bile.pdf', width = 5, height = 5)

#### enrich ####
path_info <- read.delim("cpd2path_enrichment.tsv")

eKEGG <- diff %>% 
  filter(enriched != 'none') %>% 
  filter(KEGG != "") %>% 
  pull(KEGG) %>% 
  map(~ unlist(strsplit(.x, ";"))) %>% 
  flatten_chr() %>% 
  enricher(TERM2GENE = path_info, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame
write.xlsx(eKEGG, "diff.eKEGG.xlsx")

kegg_info <- read.delim('path_description') %>% 
  filter(category == 'Metabolism')

path_name <- c('map00997', 'map01040', 'map01230', 'map00120', 'map00920', 
               'map01200', 'map00400', 'map00121', 'map00310', 'map00100',
               'map00350', 'map00460', 'map00270')

plot_data <- eKEGG %>%
  mutate(path = map_vec(strsplit(ID, ':'), \(x) x[1])) %>% 
  filter(path %in% path_name) %>% 
  arrange(FoldEnrichment)

ggscatter(plot_data, "ID", "FoldEnrichment", fill = "pvalue", 
          rotate = T, size = 'FoldEnrichment',
          shape = 21, legend = 'right', xlab = '', ylab = 'Fold Enrichment') +
  scale_fill_viridis_c(begin = .6) +
  scale_size_continuous(range = c(3, 6)) +
  theme(aspect.ratio = 2, 
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'))
ggsave('diff.eKEGG.pdf', width = 8, height = 4.5)

