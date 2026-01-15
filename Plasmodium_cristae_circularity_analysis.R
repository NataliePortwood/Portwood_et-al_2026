---
setwd("C:/Users/Natalie/Desktop")
title: "R Notebook"
output: html_notebook
---
```{r}
library(tidyverse)
library(cowplot)

EMImages <- read.csv("C:/Users/Natalie/Desktop/R studio/Book2.csv")
cols <- c("ImageNum", "Length", "Width", "Ratio", "Odd", "MitoNum", "Remarks")
names(EMImages) <- cols

```

# Data Overview - Without the Key
```{r}

# Plot the distribution of the ratios 
ggplot(EMImages, aes(x = Ratio)) + geom_histogram(binwidth = 0.1, fill = "white", color = "black") + labs(title = "Plasmodium Mitochondrial Cristae Circularity (Treatment & Control)", x = "Ratio (Long side / Short side)", y = "Frequency") +
  theme_cowplot()
```

```{r}
# Plot as a scatter plot, length on x-axis, width on y-axis
ggplot(EMImages, aes(x = Length, y = Width)) + geom_point() + labs(title = "Cristae Length vs Width", x = "Length (um)", y = "Width (um)") +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
```
```{r}
# remove the one outlier (image 34, length = 0.4)
EMImages_outlierRemoved <- EMImages[EMImages$Length != 0.4,]
# Plot again
ggplot(EMImages_outlierRemoved, aes(x = Length, y = Width)) + geom_point() + labs(title = "Cristae Length vs Width", x = "Length (um)", y = "Width (um)") +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
```

```{r}
# Now with the image key:
treatmentImages <- c(1,3,4,5,6,9,10,11,12,13,14,16,17,19,20,21,24,26,28,30,32,
                     34,34,35,36,39,42,50,51,52,53,54,56,57,58,60,61)

EMImages$Treatment <- ifelse(EMImages$ImageNum %in% treatmentImages, "Treatment", "Control")
```


```{r}
# Plot, segregating data based on treatment
ggplot(EMImages, aes(x = Length, y = Width, color = Treatment)) + geom_point() + labs(title = "Cristae Length vs Width", x = "Length (um)", y = "Width (um)") +
  theme_cowplot() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("black", "#009E73"))
```

```{r}
# Plot the distributions of the ratios overlapping. Make transparent to see overlap
# Plto density instead of histogram
density <- ggplot(EMImages, aes(x = Ratio, fill = Treatment)) + geom_density(alpha = 0.4) + labs(title = "Plasmodium Mitochondrial Cristae Circularity", x = "Ratio (Long side / Short side)", y = "Density") +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#009E73")) +
  scale_x_continuous(limits = c(0.5, 5))

density

ggsave("densityPlot.png", width = 6, height = 4, dpi = 300)
```

```{r}
# Kolmogorov-Smirnov test
ks.test(EMImages$Ratio[EMImages$Treatment == "Treatment"], EMImages$Ratio[EMImages$Treatment == "Control"])
```

# Log Normalization - All data
```{R}
# Try normalizing the data - log transform
# Exclude NA values
EMImages <- EMImages[!is.na(EMImages$Ratio),]
EMImages$LogRatio <- log(EMImages$Ratio)

ggplot(EMImages, aes(x = LogRatio, fill = Treatment)) + geom_density(alpha = 0.35) + labs(title = "Mitochondrial Cristae Circularity", x = "Log Ratio (Long side / Short side)", y = "Density") +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#009E73")) +
  scale_x_continuous(limits = c(-0.3, 1.7))

ggsave("densityPlot_log.png", width = 6, height = 4, dpi = 300)
```
# Log Normalized - Remove 'Fused' Cristae
```{r}
# Pull out data points where remarks regex matches 'Fused'
EMImages_fusedremoved <- EMImages[!grepl("Fused|fused", EMImages$Remarks),]

#relabel the fac

# Show the log transformed data
cristae_plot <- ggplot(EMImages_fusedremoved, aes(x = Ratio, fill = Treatment)) + 
  geom_density(alpha = 0.35) + 
  labs(
    title = "Mitochondrial Cristae Circularity", 
    x = "Ratio (Long side / Short side)", 
    y = "Density",
    fill = "tralopyril"   # This changes the legend title
  ) +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#76069A")) +
  scale_x_continuous(
    breaks = seq(0.75, 2, by = 0.25),
    limits = c(0.75, 2),
    expand = c(0, 0)
  )

print(cristae_plot)
```

## Also plot the non-normalized data
```{r}
# two-sided, unpaired t-test
tRes <- t.test(EMImages_fusedremoved[EMImages_fusedremoved$Treatment == "Tralopyril",]$Ratio, EMImages_fusedremoved[EMImages_fusedremoved$Treatment == "Control",]$Ratio)

fusedRemovedNonNorm <- ggplot(EMImages_fusedremoved, aes(x = Ratio, fill = Treatment)) + geom_density(alpha = 0.35) + labs(title = "Mitochondrial Cristae Circularity", x = "Ratio (Long side / Short side)", y = "Density") +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#760"))

pval      <- tRes$p.value
sig_label <- ifelse(pval < 0.001, "***",
              ifelse(pval < 0.01,  "**",
              ifelse(pval < 0.05,  "*",  "n.s.")))
fusedRemovedNonNorm

fusedRemovedNonNorm +
  annotate("text",
           x     = 1.25,
           y     = 2.5 * 1.06,  
           label = sig_label,
           size  = 8)  +
  annotate("segment",
           x     = 1.13,
           xend  = 1.37,
           y     = 2.5 * 1.05,
           yend  = 2.5 * 1.05,
           size  = 0.5) +
  scale_x_continuous(limits = c(0.75, 2)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
  
  
  


fusedRemovedNonNorm
ggsave('fusedRemovedNonNorm.png', width = 6, height = 4, dpi = 300)
```
# Number of cristae per mitochondrion
```{r}
# Within each image, count the number of cristae per mitochondrion based off
# the MitoNum value
cristaeCounts <- EMImages %>% group_by(ImageNum, Treatment, MitoNum) %>% summarise(nCristae = n()) %>% filter(ImageNum > 12)

cristaeCounts <- cristaeCounts %>% group_by(ImageNum, Treatment) %>% summarise(meanCristae = mean(nCristae)) 

cristaePlot <- ggplot(cristaeCounts, aes(x = meanCristae, fill = Treatment)) + geom_density(alpha = 0.35) + labs(title = "Number of Cristae per Mitochondrion (including fused)", x = "Number of Cristae", y = "Density") +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#76069A")) +
  scale_x_continuous(limits = c(0, 7))

cristaePlot
```

```{r}
ks.test(cristaeCounts[cristaeCounts$Treatment == "Treatment",]$meanCristae, cristaeCounts[cristaeCounts$Treatment == "Control",]$meanCristae)
```

# Number of fused cristae control vs treatment
```{r}
# Count the number of fused cristae per mitochondrion per image
fusedCounts <- EMImages %>% group_by(ImageNum, Treatment) %>% summarise(nFused = sum(grepl("Fused|fused", Remarks))) %>% filter(ImageNum > 12)

fusedPlot <- ggplot(fusedCounts, aes(x = nFused, fill = Treatment)) + geom_density(alpha = 0.35) + labs(title = "Number of Fused Cristae per Image", x = "Number of Fused Cristae", y = "Density") +
  theme_cowplot() +
  scale_fill_manual(values = c("black", "#009E73")) +
  scale_x_continuous(limits = c(0, 3))

fusedPlot
```

```{r}
ks.test(fusedCounts[fusedCounts$Treatment == "Treatment",]$nFused, fusedCounts[fusedCounts$Treatment == "Control",]$nFused)
```
