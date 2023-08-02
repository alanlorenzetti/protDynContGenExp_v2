# alorenzetti 20200824

# description ####
# this script will take previously
# created objects and parse them
# making them ready to generate
# scatter plots

# spectronaut counts vs. tpms ####
# creating spectronaut counts wide dataset
spectroWide = spectroLong %>% 
  filter(libType != "ribo") %>% 
  select(-se_abundance, -nruns, -libType) %>% 
  pivot_wider(values_from = "mean_abundance",
              names_from = "timepoint",
              names_prefix = "mean_abundance_protein_lysate_") %>% 
  left_join(., dict, by = c("locus_tag" = "query_id")) %>% 
  mutate(locus_tag = subject_id) %>% 
  select(-subject_id) %>%
  filter(!is.na(locus_tag))
  
# joining protein counts to total rna counts
abundDF = left_join(spectroWide, totrnatpm, by = "locus_tag")

# performing quantile normalization
abundDFnorm = abundDF %>%
  select(-locus_tag) %>% 
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble()

abundDFnorm = abundDFnorm %>% 
  mutate(locus_tag = abundDF$locus_tag) %>% 
  select(locus_tag, everything())
colnames(abundDFnorm) = colnames(abundDF)

# making a longer version of this dataset
abundDFnormLong = abundDFnorm %>% 
  pivot_longer(cols = -locus_tag,
               names_pattern = "mean_abundance_(.*_.*)_(.*)",
               names_to = c("libtype", "timepoint"), 
               values_to = "mean_abundance") %>% 
  distinct() %>% 
  pivot_wider(names_from = "libtype",
              values_from = "mean_abundance")

# plotting scatters for each one of the time points
abundCor = abundDFnormLong %>% 
  ggplot(aes(y = log10(protein_lysate), x = log10(rna_total))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = "#4E79A7") +
  stat_cor() +
  facet_wrap(~timepoint) +
  ylab("Log<sub>10</sub>(Protein Abundance)") +
  xlab("Log<sub>10</sub>(TPM)")

# number of entries with protein 
# and mRNA levels per time point
# abundDFnormLong %>% drop_na() %>% group_by(timepoint) %>% summarise(n = n())

# saving previous plot
ggsave("./plot/abundanceCorrelation.svg",
       plot=abundCor,
       width = 5, height = 5.5)

ggsave("./plot/abundanceCorrelation.png",
       plot=abundCor, dpi = 600,
       width = 5, height = 5.5)

ggsave("./plot/abundanceCorrelation.tiff",
       plot=abundCor, dpi = 600, compression = "lzw",
       width = 5, height = 5.5)

ggsave("./plot/abundanceCorrelation.pdf",
       plot=abundCor,
       width = 5, height = 5.5)

# oneomics relative changes vs. deseq2 relative changes ####
# parsing oneomics dataset
oneomicsLong2Parsed = oneomicsLong2 %>% 
  mutate(sigStatus = case_when(sigStatus == "down" ~ "yes",
                               sigStatus == "up" ~ "yes",
                               TRUE ~ "no"),
         timepoint = str_replace(timepoint, " vs.*$", "")) %>% 
  select(locus_tag, timepoint, lfc, padj) %>% 
  rename_with(.cols = starts_with("lfc"),
              .fn = ~ str_replace(., "lfc", "log2FoldChange_protein")) %>% 
  rename_with(.cols = starts_with("padj"),
              .fn = ~ str_replace(., "padj", "padj_protein")) %>% 
  pivot_wider(names_from = c("timepoint"),
              values_from = c("log2FoldChange_protein", "padj_protein")) %>% 
  mutate(sigStatus_protein_TP2 = case_when(abs(log2FoldChange_protein_TP2) >= log2fcthreshold &
                                             padj_protein_TP2 < padjthreshold ~ "yes",
                                           TRUE ~ "no"),
         sigStatus_protein_TP3 = case_when(abs(log2FoldChange_protein_TP3) >= log2fcthreshold &
                                             padj_protein_TP3 < padjthreshold ~ "yes",
                                           TRUE ~ "no"),
         sigStatus_protein_TP4 = case_when(abs(log2FoldChange_protein_TP4) >= log2fcthreshold &
                                             padj_protein_TP4 < padjthreshold ~ "yes",
                                           TRUE ~ "no")) %>% 
  left_join(., dict, by = c("locus_tag" = "query_id")) %>% 
  mutate(locus_tag = subject_id) %>% 
  select(locus_tag,
         everything(),
         -subject_id) %>% 
  filter(!is.na(locus_tag))

# joining datasets
relativeChangeDF = left_join(oneomicsLong2Parsed, unifiedFinFil, by = "locus_tag") %>%
  select(locus_tag,
         starts_with("log2FoldChange"),
         starts_with("padj"))

# performing quantile normalization of log2foldchange values
relativeChangeDFNorm = relativeChangeDF %>%
  select(starts_with("log2FoldChange")) %>%
  as.matrix() %>%
  normalize.quantiles() %>%
  as_tibble()

cidx = grepl("^log2FoldChange", colnames(relativeChangeDF)) %>% which()

relativeChangeDF[,cidx] = relativeChangeDFNorm

# arranging plots
p = list()

# TP2
p[["TP2"]] = relativeChangeDF %>%
  filter(abs(log2FoldChange_protein_TP2) >= log2fcthreshold &
           abs(log2FoldChange_mRNA_TP2) >= log2fcthreshold &
           padj_protein_TP2 < padjthreshold &
           padj_mRNA_TP2 < padjthreshold) %>%
  ggplot(aes(y = log2FoldChange_protein_TP2, x = log2FoldChange_mRNA_TP2)) +
  geom_point(alpha = 0.25) +
  geom_smooth(formula = y ~ x, method = "lm",
              color = "#4E79A7") +
  stat_cor() +
  xlim(c(-6,6)) + ylim(c(-6,6)) +
  xlab("mRNA Log<sub>2</sub>(Fold Change)") +
  ylab("Protein Log<sub>2</sub>(Fold Change)") +
  ggtitle("A")

# TP3
p[["TP3"]] = relativeChangeDF %>%
  filter(abs(log2FoldChange_protein_TP3) >= log2fcthreshold &
           abs(log2FoldChange_mRNA_TP3) >= log2fcthreshold &
           padj_protein_TP3 < padjthreshold &
           padj_mRNA_TP3 < padjthreshold) %>%
  ggplot(aes(y = log2FoldChange_protein_TP3, x = log2FoldChange_mRNA_TP3)) +
  geom_point(alpha = 0.25) +
  geom_smooth(formula = y ~ x, method = "lm",
              color = "#4E79A7") +
  stat_cor() +
  xlim(c(-6,6)) + ylim(c(-6,6)) +
  xlab("mRNA Log<sub>2</sub>(Fold Change)") +
  ylab("Protein Log<sub>2</sub>(Fold Change)") +
  ggtitle("B")

# TP4
p[["TP4"]] = relativeChangeDF %>%
  filter(abs(log2FoldChange_protein_TP4) >= log2fcthreshold &
           abs(log2FoldChange_mRNA_TP4) >= log2fcthreshold &
           padj_protein_TP4 < padjthreshold &
           padj_mRNA_TP4 < padjthreshold) %>%
  ggplot(aes(y = log2FoldChange_protein_TP4, x = log2FoldChange_mRNA_TP4)) +
  geom_point(alpha = 0.25) +
  geom_smooth(formula = y ~ x, method = "lm",
              color = "#4E79A7") +
  stat_cor() +
  xlim(c(-6,6)) + ylim(c(-6,6)) +
  xlab("mRNA Log<sub>2</sub>(Fold Change)") +
  ylab("Protein Log<sub>2</sub>(Fold Change)") +
  ggtitle("C")

relChPlot = ggarrange(plotlist = p, ncol = 3, nrow = 1)

# saving plot
ggsave("./plot/relativeChangeCorrelation.svg",
       plot=relChPlot,
       width = 7, height = 3)

