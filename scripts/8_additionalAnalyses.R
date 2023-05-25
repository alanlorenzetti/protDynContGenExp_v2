# alorenzetti 20230513
# description ####
# in this script, we will
# do a few sanity checks
# before submission
# and comparison of mRNA
# levels of detected vs undetected proteins

# getting started ####
# path to protein database used by proteomics analysis
protFile = "./data/Halobacterium-20080205_VNG_cRAP_TargDecoy_plusRT.fasta"
original_protein_names = tibble(locus_tag = readAAStringSet(protFile) %>% names()) %>% 
  mutate(product = str_replace(locus_tag, "^(.*?) (.*)$", "\\2"),
         locus_tag = str_replace(locus_tag, "^(.*?) (.*)$", "\\1"))
original_protein_names = original_protein_names[1:2646,]

# proteins detected by spectronaut 
proteins_detected_by_spectronaut = strsplit(x = spectro1$PG.ProteinAccessions, split = ",") %>%
  unlist() %>%
  sort() %>%
  unique()
proteins_detected_by_spectronaut = proteins_detected_by_spectronaut[-1]

# writing additional tables
# for manual inspection
ori_prots_w_repr = dict_complete %>%
  filter(locus_tag %in% original_protein_names$locus_tag) 

ori_prots_repetitive = dict_complete %>%
  filter(locus_tag %in% original_protein_names$locus_tag) %>% 
  select(representative, locus_tag) %>% 
  group_by(representative) %>% 
  summarise(n = n(), locus_tag = paste0(locus_tag, collapse = ",")) %>% 
  filter(n >= 2)
  
detected_by_spectro = ori_prots_w_repr %>%
  filter(locus_tag %in% proteins_detected_by_spectronaut)
detected_by_spectro$representative %>% unique() %>% length()

not_detected_by_spectro = ori_prots_w_repr %>%
  filter(!locus_tag %in% proteins_detected_by_spectronaut)
not_detected_by_spectro$representative %>% unique() %>% length()

# creating a dataframe to plot
# mRNA levels of detected vs non detected
abundance_detect_status_ttest = totrnatpmlong %>% 
  filter(locus_tag %in% detected_by_spectro$representative |
           locus_tag %in% not_detected_by_spectro$representative) %>% 
  mutate(detection_status = case_when(locus_tag %in% detected_by_spectro$representative ~ "detected",
                                      locus_tag %in% not_detected_by_spectro$representative ~ "not detected",
                                      TRUE ~ NA_character_))

detected_in_ttest = abundance_detect_status_ttest %>% filter(detection_status == "detected") %>% pull(locus_tag) %>% unique()
detected_in_ttest = ori_prots_w_repr %>% 
  filter(representative %in% detected_in_ttest)
detected_in_ttest$representative %>% unique() %>% length()

not_detected_in_ttest = abundance_detect_status_ttest %>% filter(detection_status == "not detected") %>% pull(locus_tag) %>% unique()
not_detected_in_ttest = ori_prots_w_repr %>% 
  filter(representative %in% not_detected_in_ttest)
not_detected_in_ttest$representative %>% unique() %>% length()

detected_absent = ori_prots_w_repr %>% 
  filter(
    representative %in% unique(detected_by_spectro$representative)[!unique(detected_by_spectro$representative) %in% unique(detected_in_ttest$representative)]
)

not_detected_absent = ori_prots_w_repr %>% 
  filter(
    representative %in% unique(not_detected_by_spectro$representative)[!unique(not_detected_by_spectro$representative) %in% unique(not_detected_in_ttest$representative)]
)

totrnatpmlong %>% 
  filter(locus_tag %in% detected_absent$representative)

totrnatpmlong %>% 
  filter(locus_tag %in% not_detected_absent$representative)

write.xlsx(x = list(ori_prots_w_repr = ori_prots_w_repr,
                    ori_prots_repetitive = ori_prots_repetitive,
                    detected_by_spectro = detected_by_spectro,
                    not_detected_by_spectro = not_detected_by_spectro,
                    detected_in_ttest = detected_in_ttest,
                    not_detected_in_ttest = not_detected_in_ttest,
                    detected_absent = detected_absent,
                    not_detected_absent = not_detected_absent),
           file = "results/sets_for_manual_inspection.xlsx",
           overwrite = T)

# plotting
abundance_detect_status_ttest_plot = abundance_detect_status_ttest %>% 
  ggplot(aes(y = log10(abundance), x = detection_status)) +
  geom_violin() +
  geom_jitter(width = 0.1, size = 0.1, alpha = 0.1) +
  facet_wrap(~ timepoint) +
  stat_compare_means(label.y = 4,
                     label.x.npc = "center",
                     method = "t.test", label = "p.signif") +
  ylab("Log<sub>10</sub>(TPM)") +
  xlab("Detection status") +
  scale_x_discrete(labels = c("detected" = "Detected",
                              "not detected" = "Undetected"))

# saving previous plot
ggsave("./plot/protAbundanceTtest.svg",
       plot=abundance_detect_status_ttest_plot,
       width = 5, height = 5.5)

ggsave("./plot/protAbundanceTtest.png",
       plot=abundance_detect_status_ttest_plot, dpi = 600,
       width = 5, height = 5.5)

# combining previous plot with abundance plot
panel_cor_protabun = ggarrange(plotlist = list(abundCor,
                                               abundance_detect_status_ttest_plot),
                               ncol = 1, nrow = 2, labels = "AUTO")

# saving previous plot
ggsave("./plot/panel_protAbundanceTtest_and_cor.svg",
       plot=panel_cor_protabun,
       width = 5, height = 11)

ggsave("./plot/panel_protAbundanceTtest_and_cor.png",
       plot=panel_cor_protabun, dpi = 600,
       width = 5, height = 11)

# number of entries for each time point
# totrnatpmlong %>%
#   filter(locus_tag %in% proteins_detected_by_spectronaut_in_dict_new_lt |
#            locus_tag %in% proteins_not_detected_by_spectronaut_in_dict_new_lt) %>%
#   mutate(detection_status = ifelse(locus_tag %in% proteins_not_detected_by_spectronaut_in_dict_new_lt, "not detected", "detected")) %>%
#   group_by(timepoint, detection_status) %>%
#   summarise(n = n(), mean = mean(abundance))
