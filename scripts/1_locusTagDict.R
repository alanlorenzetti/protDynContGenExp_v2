# alorenzetti 20200617

# description #####
# this script will generate a dictionary of locus_tag
# using sequence similarity search (rbh)
# the goal is to convert names from original halo genome annotation 
# to the nonredundant transcriptome (part of runKallisto scripts)
# derived from the third party annotation effort by pfeiffer et al 2019
# please check https://github.com/alanlorenzetti/runKallisto to see
# the workflow used to generate countTables

# following will be done
# 1. get halo proteins used in the SWATH-MS search database
# 2. get proteins from transcript database created using pfeiffer2019 annotation
# 3. run blast to generate the correspondence dictionary

# for the record: pfeifer2019 third party annotation
# https://www.ncbi.nlm.nih.gov/nuccore/BK010829
# https://www.ncbi.nlm.nih.gov/nuccore/BK010830
# https://www.ncbi.nlm.nih.gov/nuccore/BK010831

# loading libs #####
source("./scripts/0_loadingLibs.R")

# our in-house non redundant transcriptome is being loaded here
pfeiFile = "./data/pfeiGen.fa"

# path to protein database used by proteomics analysis
protFile = "./data/Halobacterium-20080205_VNG_cRAP_TargDecoy_plusRT.fasta"

# reading the protein database used in proteomics analysis
protSeqsAA = readAAStringSet(protFile)

# removing all decoys and prot seqs
# keeping only halo ORFs
protSeqsAA = protSeqsAA[1:2646]

# adjusting names and writing sequences
names(protSeqsAA) = sub(" .*$", "", names(protSeqsAA), perl = T)
writeXStringSet(protSeqsAA, "./data/protGen.fa", format = "fasta")

##############PFEI
# reading non redundant transcriptome file
pfeiSeqs = readDNAStringSet(filepath = pfeiFile, format = "fasta", use.names = T)

# removing non coding entries
pfeiSeqs = pfeiSeqs[!str_detect(names(pfeiSeqs), "VNG_t|VNG_r|VNG_s")]

# converting to aa based on 11 table and writing
arc = getGeneticCode(id_or_name2 = "11")
pfeiSeqsAA = translate(pfeiSeqs, genetic.code = arc)
writeXStringSet(pfeiSeqsAA, "./data/pfeiGenAA.fa", format = "fasta")

# RBBH between datasets
if(!file.exists("./data/res.RData")){
res = blast_best(query_file = "./data/protGen.fa",
                  subject_file = "./data/pfeiGenAA.fa",
                  seq_type = "protein")
save(res, file = "./data/res.RData")
} else {
  load("./data/res.RData")
}

# allowing correspondence only if
# identity >= .98 and subsetting
# locus_tag cols
resFil = res %>%
  filter(perc_identity >= 98) %>% 
  select(query_id, subject_id) %>% 
  ungroup()

# adjusting dict names
dict = resFil %>% 
  mutate(query_id = sub("\\|.*$", "", query_id),
         subject_id = sub("\\|.*$", "", subject_id)) %>% 
  mutate(query_id = sub("-.*$", "", query_id)) %>% 
  mutate(query_id = sub("_", "", query_id))
