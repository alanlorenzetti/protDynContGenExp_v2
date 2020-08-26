# alorenzetti 202008

# description ####
# this script will take oneomics PCA dataset
# with the first two principal components
# (those with the highest variances)

# loading libs #####
source("./scripts/0_loadingLibs.R")

# loading file ####
pcaOneOmics = read_csv("./data/20200729_score_plot_PCA-PCVG.csv")
colnames(pcaOneOmics) = c("Time Point",
                          "sample",
                          "PC1",
                          "PC2")

# wrangling to get info from sample names
pcaOneOmics = pcaOneOmics %>% 
  mutate(`Bio. Replicate` = str_replace(sample, "^.*(BR[1-3]).*", "\\1"),
         `Run` = str_replace(sample, "^.*r[0]([1-3]).wiff.*$", "\\1"))

# plotting
svglite("plot/oneomicsPCA.svg", width = 6, height = 3.5)
pcaOneOmics %>% 
  ggplot(aes(x=PC1, y=PC2,
             color = `Time Point`,
             shape = `Bio. Replicate`)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("TP1" = "#E15759",
                                "TP2" = "#F28E2B",
                                "TP3" = "#4E79A7",
                                "TP4" = "#59A14F")) +
  guides(colour = guide_legend(order = 1), 
         shape = guide_legend(order = 2)) +
  ggtitle("A")
dev.off()