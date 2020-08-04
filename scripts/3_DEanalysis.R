# alorenze 20200617
# this script should be run after
# 1_locusTag.R
# 2_parseOneOmicsResults.R

# setting thresholds
padjthreshold = 0.05
log2fcthreshold = 0.75
borderlinezero = 0.25

############ loading packages
# biocmanager is required for loading and installing packages
if(!require("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# required libs
packs = c("tidyverse",
          "DESeq2",
          "tximport")

# loading and installing packs if needed
p_load(char=packs)

############# loading files
# those have to be generated using a whole distinct pipeline
# https://github.com/alanlorenzetti/runKallisto
# reading totalrna counts
totrna=read_delim("./data/tableEstCountsTotalRNA23.tsv",delim="\t")

# reading riboseq counts
riborna=read_delim("./data/tableEstCountsRiboSeqTrim15.tsv",delim="\t")

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

assay = select(totrna, starts_with("total")) %>%
  as.matrix() %>% 
  round(digits = 0)

totrnaSE = SummarizedExperiment(assay = list(counts=assay),
                                rowData = totrna[,c(1:3)],
                                colData = colData)
rownames(totrnaSE) = rowData(totrnaSE)$target_id

# riborna
samples = riborna %>% 
  select(starts_with("ribosomal")) %>%
  colnames()
tp = samples %>% 
  sub("^.*_RNA_[1-3]-([1-4])_S.*$", "\\1", .)
bioRep = samples %>% 
  sub("^.*_RNA_([1-3])-[1-4]_S.*$", "\\1", .)

colData = data.frame(row.names = samples,
                     timepoint = tp,
                     replicate = bioRep)

assay = select(riborna, starts_with("ribosomal")) %>%
  as.matrix() %>% 
  round(digits = 0)

ribornaSE = SummarizedExperiment(assay = list(counts=assay),
                                rowData = riborna[,c(1:3)],
                                colData = colData)
rownames(ribornaSE) = rowData(ribornaSE)$target_id

# creating deseq2 objects
totrnadds = totrnaSE
totrnadds = DESeqDataSet(totrnadds, design = ~ timepoint)

ribodds = ribornaSE
ribodds = DESeqDataSet(ribodds, design = ~ timepoint)

# removing genes with zero counts and performing DESeq2 analysis
totrnadds = totrnadds[rowSums(counts(totrnadds)) > 1, ]
totrnadds = DESeq(totrnadds)

ribodds = ribodds[rowSums(counts(ribodds)) > 1, ]
ribodds = DESeq(ribodds)

# result tables for contrasts
results = list()
for(type in c("totrna", "ribo")){
  for(i in 2:4 %>% as.character()){
    results[[paste0(type,i)]] = results(get(paste0(type,"dds")),
                                        contrast= c("timepoint", i, "1"),
                                        alpha = padjthreshold) %>% 
      as_tibble(rownames = "target_id")
  }
}

# final result tables
# containing sig status and borderline status
resultsFin = list()
for(type in c("totrna", "ribo")){
  for(i in 2:4 %>% as.character()){
    varName = paste0("sig",type,"OnTP",i)
    varName2 = paste0("borderlineZero",type,"OnTP",i)
    resultsFin[[paste0(type,i)]] = results[[paste0(type,i)]] %>% 
      mutate(!!varName := case_when(abs(log2FoldChange) >= log2fcthreshold & padj < padjthreshold ~ "yes",
                                    TRUE ~ "no")) %>% 
      mutate(!!varName2 := case_when(abs(log2FoldChange) < borderlinezero ~ "yes", 
                                     TRUE ~ "no"))
  }
}

# unifying final results
unifiedFin = list()
for(i in 2:4 %>% as.character()){
  unifiedFin[[paste0("TP",i)]] = inner_join(resultsFin[[paste0("totrna",i)]],
                                            resultsFin[[paste0("ribo",i)]],
                                            by = "target_id") %>% 
    rename("id" = target_id) %>% 
    mutate(locus_tag = sub("\\|.*$", "", id))
}
