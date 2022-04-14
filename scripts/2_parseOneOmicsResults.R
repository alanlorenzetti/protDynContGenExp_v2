# alorenzetti 202008

# description ######
# this script will take the proteomics
# files output by OneOmics software
# and perform wrangling and filtering
# in this version
# we are using OneOmics output
# combining all the biological
# replicates, therefore, only one
# input table

# thresholds
padj = 0.05
pthr = padj
lfcthr = 0.75

# loading files #####
oneomics = read_xlsx("./data/20200729_HL_VNGcplte_3BRs_RTCal_heatmap_0conf_0reproducibility_all2070proteins_s.xlsx") %>% 
  select(-secondary_name)

# parsing colnames
colnames(oneomics) = colnames(oneomics)
colnames(oneomics) = gsub(" ", "_", colnames(oneomics))
colnames(oneomics) = gsub("log2_sfc", "lfc", colnames(oneomics))
colnames(oneomics) = sub("_vs._TP1", "", colnames(oneomics))
colnames(oneomics)[-1] = sub("(.*)_(.*)", "\\2_\\1", colnames(oneomics)[-1])
colnames(oneomics)[1] = "locus_tag"

# adding BH adjusted pval
oneomics = oneomics %>% 
  mutate(padj_TP2 = p.adjust(pval_TP2, method = "BH"),
         padj_TP3 = p.adjust(pval_TP3, method = "BH"),
         padj_TP4 = p.adjust(pval_TP4, method = "BH"))

# filtering based on
# c (confidence >= 0.75 )
# mlr >= 0.2
cthr = 0.75
mlrthr = 0.2

oneomics = oneomics %>%
  filter((c_TP2 >= cthr & mlr_TP2 >= mlrthr) |
           c_TP3 >= cthr & mlr_TP3 >= mlrthr |
           c_TP4 >= cthr & mlr_TP4 >= mlrthr)

# pivoting object
oneomicsLong = oneomics %>% 
  pivot_longer(cols = contains("TP"),
               names_to = c("type", "timepoint"),
               names_sep = "_")

# filtering out attributes that are not
# going to be used anymore
oneomicsWide = oneomics %>% 
  select(locus_tag, contains("lfc"), contains("padj"))

colnames(oneomicsWide)[-1] = sub("lfc_", "mean_lfc_protein_lysate_", colnames(oneomicsWide)[-1])
colnames(oneomicsWide)[-1] = sub("padj_", "mean_padj_protein_lysate_", colnames(oneomicsWide)[-1])

# adjusting locus_tags
# according to dictProd
# "VNG5199H"  "VNG0606G"  "VNG0779C"  "VNG1585Cm" "VNG0780H"  "VNG6339H"
# don't have a matching sequence in pfeiLocusTag
# those are going to be removed
oneomicsWide = left_join(oneomicsWide, dict, by = c("locus_tag" = "query_id"))

remove = oneomicsWide$locus_tag[is.na(oneomicsWide$subject_id) %>% which()]
oneomicsWide = oneomicsWide[!(oneomicsWide$locus_tag %in% remove),] %>% 
  mutate(locus_tag = subject_id) %>% 
  select(-subject_id)

oneomicsLong2 = oneomicsLong %>% 
  pivot_wider(names_from = type) %>% 
  mutate(sigStatus = case_when(lfc >= lfcthr & padj < pthr ~ "up",
                               lfc <= -lfcthr & padj < pthr ~ "down",
                               TRUE ~ "no"))

oneomicsLong2$timepoint = str_replace(oneomicsLong2$timepoint, "$", " vs. TP1")

# saving OneOmics long table
# shared with only on Apr 13th 2022
if(!dir.exists("results")){dir.create("results")}

write.xlsx(x = oneomicsLong2,
           file = "results/20220413_oneomics_long.xlsx",
           overwrite = T)

# plots ####
# they will be created and unified
# in the next script
oneomicsplots = list()

# volcano plots
#svglite("plot/oneomicsVolcanoPlot.svg", width = 8, height = 2.5)
oneomicsplots$volcanos = oneomicsLong2 %>%
  ggplot(aes(x = lfc, y = -log10(padj), color = sigStatus)) +
  geom_point2(alpha = 0.25, show.legend = F) +
  facet_wrap(~ timepoint) +
  xlim(c(-6,6)) +
  scale_color_manual(values = c("up" = "#E15759",
                                "down" = "#4E79A7",
                                "no" = "grey30")) +
  xlab("Log<sub>2</sub>(Fold Change)") +
  ylab("-Log<sub>10</sub>(Adjusted *P*)") +
  ggtitle("C")
#dev.off()

# plotting heatmap
heatColors = colorRamp2(c(-4,0,4), colors = c("#4E79A7", "white", "#E15759"))
M = oneomics %>% select(contains("lfc")) %>% as.matrix()
colnames(M) = str_replace(colnames(M), "lfc_", "")

#svglite("plot/oneomicsHeatmap.svg", width = 3, height = 4.5)
oneomicsplots$ht = Heatmap(M,
                           col = heatColors,
                           border = T,
                           heatmap_legend_param = list(
                             at = c(-6, 0, 6),
                             title = "LFC",
                             border = T,
                             title_gp = gpar(fontsize = 10, fontfamily = "Arial"), 
                             labels_gp = gpar(fontsize = 10, fontfamily = "Arial")
                           ),
                           column_names_gp = gpar(fontsize = 10, fontfamily = "Arial"))
#dev.off()

oneomicsplots$drawnht = grid.grabExpr(draw(oneomicsplots$ht,
                                           heatmap_legend_side = "right"))

# getting number of differentially expressed genes
oneomicsLong2 %>% filter(timepoint == "TP2 vs. TP1") %>% select(sigStatus) %>% table()
oneomicsLong2 %>% filter(timepoint == "TP3 vs. TP1") %>% select(sigStatus) %>% table()
oneomicsLong2 %>% filter(timepoint == "TP4 vs. TP1") %>% select(sigStatus) %>% table()

# plotting a venn diagram
venn = list()
venn$`TP2 vs TP1` = oneomicsLong2 %>%
  filter(timepoint == "TP2 vs. TP1") %>% 
  filter(sigStatus == "up" | sigStatus == "down") %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

venn$`TP3 vs TP1` = oneomicsLong2 %>%
  filter(timepoint == "TP3 vs. TP1") %>% 
  filter(sigStatus == "up" | sigStatus == "down") %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

venn$`TP4 vs TP1` = oneomicsLong2 %>%
  filter(timepoint == "TP4 vs. TP1") %>% 
  filter(sigStatus == "up" | sigStatus == "down") %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

plot(venn(venn), quantities = T)
