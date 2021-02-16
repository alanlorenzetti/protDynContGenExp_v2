# alorenze 20200824

# description ####
# this script will load deLomana et. al 2020
# total RNA-Seq dataset and perform
# differential expression analysis

# setting thresholds ####
padjthreshold = 0.05
log2fcthreshold = 0.75

# loading files ####
# those have to be generated using a whole distinct pipeline
# https://github.com/alanlorenzetti/runKallisto
# reading totalrna raw estimate counts
totrna=read_delim("./data/tableEstCountsTotalRNA23.tsv",delim="\t")

# reading and parsing totalrna tpms
# exploratory analysis of kallisto TPM datasets
# to manually compute TPMs, follow the instructions on the following
# page: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
totrnatpm=read_delim("data/tableTpmTotalRNA23.tsv",delim="\t")
totrnatpm["mean_abundance_rna_total_TP1"] = totrnatpm %>%
  select(matches("total-RNA-[1-3]-1")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP2"] = totrnatpm %>%
  select(matches("total-RNA-[1-3]-2")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP3"] = totrnatpm %>%
  select(matches("total-RNA-[1-3]-3")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP4"] = totrnatpm %>%
  select(matches("total-RNA-[1-3]-4")) %>% rowMeans()

totrnatpm = totrnatpm %>%
  select(target_id, contains("abundance")) %>% 
  mutate(locus_tag = sub("\\|.*$", "", target_id)) %>% 
  select(-target_id)

# starting DE analysis ####
## building se object for both
# totrna
samples = totrna %>% 
  select(starts_with("total")) %>%
  colnames()
tp = samples %>% 
  sub("^.*-RNA-[1-3]-([1-4])_S.*$", "\\1", .)
bioRep = samples %>% 
  sub("^.*-RNA-([1-3])-[1-4]_S.*$", "\\1", .)

colData = data.frame(row.names = samples,
                     timepoint = tp,
                     replicate = bioRep)

totrnaSE = SummarizedExperiment(assay = list(counts=assay),
                                rowData = totrna[,c(1:3)],
                                colData = colData)
rownames(totrnaSE) = rowData(totrnaSE)$target_id

# creating deseq2 objects
totrnadds = totrnaSE
totrnadds = DESeqDataSet(totrnadds, design = ~ timepoint)

# removing genes with zero counts and performing DESeq2 analysis
totrnadds = totrnadds[rowSums(counts(totrnadds)) > 1, ]
totrnadds = DESeq(totrnadds)

# result tables for contrasts
results = list()
type = "totrna"

for(i in 2:4 %>% as.character()){
  results[[paste0(type,i)]] = results(get(paste0(type,"dds")),
                                      contrast= c("timepoint", i, "1"),
                                      alpha = padjthreshold) %>% 
    as_tibble(rownames = "target_id")
}

# final result tables
# containing sig status and borderline status
resultsFin = list()
for(i in 2:4 %>% as.character()){
  varName = paste0("sig",type,"OnTP",i)
  resultsFin[[paste0(type,i)]] = results[[paste0(type,i)]] %>% 
    mutate(!!varName := case_when(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold ~ "yes",
                                  TRUE ~ "no"))
  resultsFin[[paste0(type,i)]] = resultsFin[[paste0(type,i)]] %>% 
    mutate(target_id = str_replace(target_id, "\\|.*$", "")) %>% 
    dplyr::rename(locus_tag = "target_id")
}

# unifying final results
unifiedFin = full_join(x = resultsFin[["totrna2"]],
                       y = resultsFin[["totrna3"]],
                       by = "locus_tag",
                       suffix = c("_TP2", "_TP3"))

unifiedFin = full_join(x = unifiedFin,
                       y = resultsFin[["totrna4"]],
                       by = "locus_tag")

colnames(unifiedFin)[16:21] = paste0(colnames(unifiedFin)[16:21], "_TP4")

# filtering unifiedFin to include only variables
# we are going to use
unifiedFinFil = unifiedFin %>% 
  select(locus_tag,
         starts_with("log2Fold"),
         starts_with("padj"),
         starts_with("sig")) %>% 
  rename_with(.cols = starts_with("log2FoldChange"),
              .fn = ~ str_replace(., "log2FoldChange", "log2FoldChange_mRNA")) %>% 
  rename_with(.cols = starts_with("padj"),
              .fn = ~ str_replace(., "padj", "padj_mRNA")) %>% 
  rename_with(.cols = starts_with("sig"),
              .fn = ~ sub("sigtotrnaOn", "sigStatus_mRNA_", .))
