# alorenzetti 20200617
# run after
# 1_locusTag.R
# 2_parseOneOmicsResults.R
# 3_DEanalysis.R
# 4_dfJoiningAndLM.R

# this script will take proteomics data objects
# and perform PCA to group proteins samples
# biological replicates and time points taken into account

# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# required libs
packs = c("tidyverse",
          "broom",
          "plotly")

# loading packs
p_load(char = packs)

# the following are not available at CRAN
# library(devtools) ; install_github("allydunham/tblhelpr")
# install_github("AckerDWM/gg3D")
library(tblhelpr)

# setting ggplot2 theme
theme_set(theme_bw())

# filtering in log2foldchange data contained in dfpca
# and dropping lines containing NA
dfpca = dfpca %>%
  select(feature, contains("log2_sfc")) %>% 
  drop_na()

# transposing dataset
dfpca = transpose_tibble(dfpca, feature, id_col = "sample")

# adapted from:
# https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
# creating a dataframe containing all pca info we need
dfpcaComplete = dfpca %>% 
  nest(data=everything()) %>% 
  mutate(pca = map(data, ~ prcomp(.x %>% select(-sample), 
                                  center = T, scale = T)),
         pca_aug = map2(pca, data, ~augment(.x, data = .y)))

# computing the explained variance for
# each of the principal components
expVar = dfpcaComplete %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains(".fittedPC")), .funs = list(~var(.))) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# preparing df to plot PCA
pcaAugDf = dfpcaComplete %>%
  unnest(pca_aug) %>%
  mutate(br = paste0("BR ", c(rep(1,3), rep(2,3), rep(3,3))),
         tp = rep(c("TP2 vs TP1","TP3 vs TP1","TP4 vs TP1"),3))
  
# 3D scatter with plotly
# setting up label titles
# rstudio is recommended for interactive visualization
f = list(
  family = "Courier New, monospace")
x = list(
  title = paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)"),
  titlefont = f
)
y = list(
  title = paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)"),
  titlefont = f
)
z = list(
  title = paste0("PC3", " (", round(expVar$var_exp[3] * 100, digits = 2), "%)"),
  titlefont = f
)

# plotting
plot_ly(data=pcaAugDf,
        x=~.fittedPC1,
        y=~.fittedPC2,
        z=~.fittedPC3,
        type="scatter3d",
        mode="markers",
        color=~tp,
        symbol=~br,
        symbols = c("circle", "square", "diamond")) %>% 
  layout(scene = list(xaxis=x,
                      yaxis=y,
                      zaxis=z))

# finding top three loadings in PC1
dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = T),1] %>% head(3)

# finding bottom three loadings in PC1
dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = F),1] %>% head(3)

# plotting 2D scatter (alternative version to the 3D scatter)
twodscatterplot = pcaAugDf %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, colour = tp, shape = br)) +
  geom_point() +
  xlab(paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)")) +
  scale_color_discrete(name="Time Point Contrast") +
  scale_shape_discrete(name="Biological Replicate")

# saving previous plot
ggsave("./plot/twodscatterplot.svg",
       plot=twodscatterplot,
       width = 7, height = 5)

# plotting explained variance
cumuvarplot = expVar %>%
  gather(key = key, value = value, var_exp:cum_var_exp) %>% 
  ggplot(aes(pc, value, group = key)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~key, scales = "free_y") +
  lims(y = c(0, 1)) +
  labs(y = "Variance",
       x = "Principal Components")

# saving previous plot
ggsave("./plot/cumuvarplot.svg",
       plot=cumuvarplot,
       width = 7, height = 5)
