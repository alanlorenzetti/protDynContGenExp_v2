# alorenze 20200617
# this script should be run after locusTagDictGbk2Tpa.R
# 1_locusTag.R
# 2_parseOneOmicsResults.R
# 3_DEanalysis.R

# setting thresholds
log2fcthreshold = 0.75
borderlinezero = 0.25

# biocmanager is required for loading and installing packages
if(!require("BiocManager")){install.packages("BiocManager"); library("BiocManager")}

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# required libs
packs = c("tidyverse",
          "ggpubr",
          "DescTools",
          "see",
          "grid",
          "ggthemes")

# loading packs
p_load(char = packs)

# setting ggplot pubr theme
theme_set(theme_pubr())

# unifying datasets
alldfs = list()
for(i in 2:4 %>% as.character()){
  varName[1] = paste0("sigProtOnTP",i)
  varName[2] = paste0("borderlineZeroProtOnTP",i)
  varName[3] = paste0("sigtotrnaOnTP",i)
  varName[4] = paste0("sigriboOnTP",i)
  varName[5] = paste0("borderlineZerototrnaOnTP",i)
  varName[6] = paste0("borderlineZeroriboOnTP",i)
  alldfs[[paste0("TP",i)]] = inner_join(tpivstp1[[paste0("TP",i)]],
                                        unifiedFin[[paste0("TP",i)]],
                                        by = c("locus_tag" = "locus_tag")) %>% 
    select(locus_tag,
           log2FoldChange.x,
           log2FoldChange.y,
           lfc, 
           !!varName[1],
           !!varName[2],
           !!varName[3],
           !!varName[4],
           !!varName[5],
           !!varName[6]) %>%
    rename("mRNA" = log2FoldChange.x,
           "RPF" = log2FoldChange.y,
           "protein" = lfc) %>% 
    mutate(`RPF-mRNA-quad` = case_when(mRNA >= log2fcthreshold & RPF >= log2fcthreshold ~ "Q1",
                                       mRNA <= -log2fcthreshold & RPF >= log2fcthreshold ~ "Q2",
                                       mRNA <= -log2fcthreshold & RPF <= -log2fcthreshold ~ "Q3",
                                       mRNA >= log2fcthreshold & RPF <= -log2fcthreshold ~ "Q4",
                                       mRNA >= borderlinezero & abs(RPF) < borderlinezero ~ "BL14",
                                       mRNA <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL23",
                                       RPF >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                       RPF <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                       TRUE ~ "center")) %>% 

    mutate(`protein-mRNA-quad` = case_when(mRNA >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                           mRNA <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                           mRNA <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                           mRNA >= log2fcthreshold & protein <= -log2fcthreshold ~ "Q4",
                                           mRNA >= borderlinezero & abs(protein) < borderlinezero ~ "BL14",
                                           mRNA <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                           protein >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                           protein <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                           TRUE ~ "center")) %>%
    
    mutate(`protein-RPF-quad` = case_when(RPF >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                          RPF <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                          RPF <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                          RPF >= log2fcthreshold & protein <= -log2fcthreshold  ~ "Q4",
                                          RPF >= borderlinezero & abs(protein) < borderlinezero ~ "BL14",
                                          RPF <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                          protein >= borderlinezero & abs(RPF) < borderlinezero ~ "BL12",
                                          protein <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL34",
                                          TRUE ~ "center")) %>% 
    
    mutate(`RPF-mRNA-quad` = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "yes" & get(varName[3]) == "no" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "no" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       abs(RPF) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                       TRUE ~ "center"),
           `protein-mRNA-quad` = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "yes" & get(varName[3]) == "no" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "no" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           abs(protein) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                           TRUE ~ "center"),
           `protein-RPF-quad` = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "yes" & get(varName[4]) == "no" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "no" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          abs(protein) < borderlinezero & abs(RPF) < borderlinezero ~ "core",
                                          TRUE ~ "center")) %>% 
    
    mutate(RPF_mRNA_sig = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                    get(varName[4]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                    get(varName[4]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                    TRUE ~ "faint"),
           protein_mRNA_sig = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                        get(varName[1]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                        get(varName[1]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                        TRUE ~ "faint"),
           protein_RPF_sig = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ "bold",
                                       get(varName[1]) == "yes" & get(varName[4]) == "no" ~ "mid",
                                       get(varName[1]) == "no" & get(varName[4]) == "yes" ~ "mid",
                                       TRUE ~ "faint")) %>% 
    
    mutate(timePoint = i %>% as.numeric())
}

# function to plot tpi vs tpi-1
plotShift = function(df, type, timepointContrast){
  
  typey=sub("^(.*)_.*$", "\\1", type)
  typex=sub("^.*_(.*)$", "\\1", type)
  
  time = timepointContrast
  
  ihigh = sub("^(.*)_vs_.*$","\\1",time)
  ilow = sub("^.*_vs_(.*$)","\\1",time)
  
  xlabel = paste0(typex,"_",ilow)
  ylabel = paste0(typey,"_",ihigh)
  colour = paste0("sigProtOn",ihigh)
  
  model = lm(formula = df[[ihigh]][,typey] %>% unlist() %>% unname() ~ df[[ilow]][,typex] %>% unlist() %>% unname())
  cortext = paste0("R^2 = ", summary(model)$r.squared %>% round(digits = 3))
  pval = pf(summary(model)$fstatistic[1],
            summary(model)$fstatistic[2],
            summary(model)$fstatistic[3],
            lower.tail = FALSE) %>% 
    formatC(format = "e", digits = 2)
  pvaltext = paste0("p = ", pval)
  
  grob = grobTree(textGrob(paste0(cortext, "; ", pvaltext),
                           gp=gpar(fontsize=8),
                           x=0.075, y=0.95, hjust=0))
  
  p = ggplot() +
    geom_point2(mapping = aes(x=df[[ilow]][,typex] %>% unlist() %>% unname(),
                             y=df[[ihigh]][,typey] %>% unlist() %>% unname()),
               alpha=0.05, size = 2, 
               show.legend = F) +
    geom_smooth(inherit.aes = F,
                mapping = aes(x=df[[ilow]][,typex] %>% unlist() %>% unname(),
                              y=df[[ihigh]][,typey] %>% unlist() %>% unname()),
                show.legend = F,
                formula = "y ~ x",
                method = "lm") +
    xlim(c(-10,10)) + ylim(c(-10,10)) +
    xlab(xlabel) + ylab(ylabel) +
    annotation_custom(grob = grob)
  
  return(p)
}

# gen for every type
types = c("RPF_mRNA",
          "protein_mRNA",
          "protein_RPF")

timepoints = c("TP3_vs_TP2",
               "TP4_vs_TP3",
               "TP2_vs_TP2",
               "TP3_vs_TP3",
               "TP4_vs_TP4")

# creating plots for alldfs
shiftPlots = list()
for(type in types){
  for(time in timepoints){
    shiftPlots[[type]][[time]] = plotShift(alldfs, type, time)
  }
}

# arranging plots
pl = ggarrange(plotlist = list(shiftPlots$protein_mRNA$TP2_vs_TP2,
                               shiftPlots$protein_mRNA$TP3_vs_TP3,
                               shiftPlots$protein_mRNA$TP4_vs_TP4,
                               shiftPlots$protein_mRNA$TP3_vs_TP2,
                               shiftPlots$protein_mRNA$TP4_vs_TP3),
               ncol = 3, nrow = 2)

# saving panel protein in function of mRNAs in TP4 vs TP1
ggsave("./plot/tp4_vs_tp3_prot_mrna_cor.svg",
       plot=shiftPlots$protein_mRNA$TP4_vs_TP3,
       width = 5, height = 5)

ggsave("./plot/prot_mrna_cor_panel.svg",
       plot=pl, 
       width = 10, height = 6)
 