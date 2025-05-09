#### Jinxin Meng, 20250325, 20250415 ####
setwd("F:/project/20250313_PDCoV_BAC_Bile_Zhangyq/MGX-16S/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, rstatix)
source("/code/R_func/profile_process.R")
source("/code/R_func/taxa.R")

#### metadata ####
grp <- c('Mock','PDCoV')
grp_col <- structure(c('#a0a0a4','#d7b0b0'), names = grp)

group <- read.delim("group.tsv")
profile <- read.delim("profile.tsv", row.names = 1, check.names = F) %>% 
  profile_transRA()
taxa <- read.delim("taxa.tsv")
tr <- treeio::read.tree('tree.nwk')

#### alpha ####
source("/code/R_func/diversity.R")

map(c('pd', 'shannon', 'richness', 'ace'), ~
  calcu_alpha(profile, method = .x, tree = tr) %>% 
    left_join(group, by = "sample") %>% 
    mutate(group = factor(group, grp)) %>% 
    ggboxplot("group", "value", fill = "group", palette = grp_col, 
            legend = "none", outlier.size = 1, xlab = "", 
            ylab = paste0(.x, ' index')) +
    stat_compare_means(comparisons = list(grp), label = "p.signif", 
                       step.increase = .06, vjust = .7, tip.length = .02) +
    theme(aspect.ratio = 2) ) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1, align = 'v')

#### PCoA ####
source("/code/R_func/plot_PCoA.R")

plot_PCoA(profile, group, group_order = grp, group_color = grp_col, 
          add_group_label = T, lab_size = 3, ellipse_level = .9, 
          show_legend = F, show_grid = T, show_line = F)

#### compos ####
source("F:/code/R_func/taxa.R")

colors <- c("#8dd3c7","#ffffb3","#80b1d3","#b3de69","#fdb462","#bc80bd","#fb8072",
            "#ffed6f","#fccde5","#bebada","#e5c494","#ccebc5","#d9d9d9")

# phylum
data <- taxa_trans(profile, taxa, group, to = 'phylum', out_all = T) %>% 
  filter(rownames(.) != 'p__Unknown')
plot_compos_manual(data, group, taxa_color = colors, group_order = grp,
                   plot_title = 'Phylum level') +
  theme(aspect.ratio = 2.4)

# genus
data <- taxa_trans(profile, taxa, group, to = 'genus', top_n = 12,
                   other_name = 'g__Other') %>% 
  filter(rownames(.) != 'g__Unknown')
plot_compos_manual(data, group, taxa_color = colors, group_order = grp, out_all = T,
                   plot_title = 'Genus level') +
  theme(aspect.ratio = 2.4)

# abundance
data <- map(c('phylum', 'genus'), ~
              taxa_trans(profile, taxa, group, to = .x, out_all = T, transRA = T) %>% 
              profile_smp2grp(group) %>% 
              rownames_to_column('name')) %>% 
  set_names(c('phylum', 'genus'))

# diff
diffs <- map(c('phylum', 'family', 'genus'), ~
               taxa_trans(profile, taxa, group, to = .x, out_all = T) %>%
               profile_transRA() %>% rownames_to_column('name') %>% 
               gather('sample', 'value', -name) %>% 
               left_join(group, by = 'sample') %>% 
               group_by(name) %>% 
               wilcox_test(value ~ group, comparisons = c('PDCoV', 'Mock'), 
                           p.adjust.method = 'BH', detailed = T) %>% 
               select(-.y.) ) %>% 
  set_names(c('phylum', 'family', 'genus'))

#### lefse ####
library(microeco)
library(magrittr)

set.seed(2025)
group <- read.delim("group.tsv") %>% 
  column_to_rownames('sample') %>% 
  mutate(group = factor(group))
profile <- read.delim("profile.tsv", row.names = 1) %>% 
  profile_transRA()
taxa <- read.delim("taxa.tsv") %>% 
  column_to_rownames('name') %>% 
  tidy_taxonomy() %>% 
  rename_all(~ str_to_title(.x))
tr <- treeio::read.tree('tree.nwk') # 有根树

dataset <- microtable$new(sample_table = group, otu_table = profile, 
                          tax_table = taxa, phylo_tree = tr)

trans_abund$new(dataset)

lefse <- trans_diff$new(dataset = dataset, method = 'lefse', group = 'group',
                        alpha = 1)

lefse$plot_diff_bar(use_number = 1:30, width = 0.6, group_order = grp, 
                    color_values = grp_col) + 
  theme_pubr()

lefse$plot_diff_cladogram(use_taxa_num = 160, use_feature_num = 40, 
                          clade_label_level = 6, group_order = grp, 
                          color = grp_col)

#### correlation ####
library(tidygraph)
library(ggraph)

metadata <- read.delim('viral_loads.txt', row.names = 1)

data <- taxa_trans(profile, taxa, to = 'genus', out_all = T) %>% 
  profile_top_n(n = 30)

corr <- psych::corr.test(t(data), metadata)

r_data <- data.frame(t(corr$r), check.names = F) %>% 
  rownames_to_column('name') %>% 
  gather('taxa', 'rval', -name)
p_data <- data.frame(t(corr$p), check.names = F) %>% 
  rownames_to_column('name') %>% 
  gather('taxa', 'pval', -name)

plot_data <- cbind(r_data, select(p_data, pval)) %>% 
  filter(taxa != 'g__Unknown') %>% 
  filter(pval < 0.05)

# graph
edges <- select(plot_data, from = name, to = taxa, rval) %>% 
  mutate(direct = ifelse(rval < 0, 'neg', 'pos'),
         rval = abs(rval))
nodes <- data.frame(node = unique(c(edges$from, edges$to))) %>% 
  mutate(type = ifelse(grepl('g__', node), 'taxa', 'virus'),
         phylum = ifelse(type == 'taxa', taxa$phylum[match(node, taxa$genus)], ''))

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- c("#8dd3c7","#bebada","#fb8072","#80b1d3",
            "#fdb462","#b3de69","#fccde5","#bc80bd","#ffed6f",
            '#1b9e77','#d95f02','#7570b3','#e7298a')

ggraph(g, layout = 'linear', circular = T) +
  geom_edge_arc(aes(edge_width = rval, edge_linetype = direct, 
                    edge_colour = direct), strength = .2) +
  scale_edge_linetype_manual(values = c(2, 1)) +
  scale_edge_width(range = c(.5, .8)) +
  scale_edge_color_manual(values = c("#7fc97f","#fdc086")) +
  geom_node_point(aes(color = phylum, shape = type), size = 3) +
  scale_shape_manual(values = c(16, 15)) +
  scale_color_manual(values = colors) +
  geom_node_text(aes(label = node), size = 2) +
  theme_void() +
  theme(aspect.ratio = 1)
ggsave('cor_network.pdf', width = 5, height = 4)
