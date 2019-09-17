t-SNE Analysis for Flow Cytometry Data
======================================

t-Distributed Stochastic Neighbor Embedding (**t-SNE**) is a non-linear technique for dimensionality reduction that is particularly well suited for the visualization of high-dimensional datasets.

This tutorial will use [FlowJo](https://www.flowjo.com/) to load FSC file(s), gate alive population and export CSV file containing compensated values. We will load the CSV file and use R packages to perform tSNE analysis.

This totorial will use 7-color staining of myeloid cells in liver of Clec4f-Cre.

Marker | Fluorescence/Channel| 
:------:|:----------------:|
DAPI (cell viability) | BUV395 |
CD45 | PerCP-Cy5.5 |
CD11b | APC-Cy7 |
F4/80 | BV510 |
Ly6C | PE-Cy7 |
Ly6G | APC |
dTomato | PE |


## FlowJo

To compare samples, we will merge the samples together into a single FCS file first (a process called concatenation), and then perform dimensionality reduction on the concatenated data set.

During the concatenation step, we will create new derived parameters based on keyword/value pairs (such as timepoints, treatment groups, or simulations).

* Clean up the data

We will first apply manual gates to exclude doublets and debris by FSC and SSC channels. Then, we will use DAPI to eliminate dead cells from each sample. This step reduces noise in the data and can improve the tSNE algorithm output.

* Downsample

It is recommended to initiate the calculation on a gated population containing no more than 20,000 events.

Choose **Workspace** -> **Plugins** -> **DownSampleV3**

* Export Data

Choose **File** -> **Export/Concatenate** -> **Export / Concatenate Populations**.

![Export](FACS_Files/Export_sample.png)

## tSNE Analysis Using R

### tSNE

```{r}
# install.packages("Rtsne")
library(Rtsne)
data2 <- data[, colnames_proj]
class(data2) <- "numeric"
set.seed(123)  # set random seed
rtsne_out <- Rtsne(data2, pca = FALSE, verbose = TRUE, perplexity = 100) ## The larger perplexity, the longer it takes to finished analysis. But the ponts will be more separated.
names(rtsne_out)
data_df <- data.frame(cellType=data[, 1], tSNE1=rtsne_out$Y[, 1], tSNE2=rtsne_out$Y[, 2])
```

### Dimention Reduction

```{r}
library(ggplot2)
library(RColorBrewer)
data_df_sample <- data_df[sample(1:nrow(data_df), 15000), ]

ggplot(data_df_sample) + geom_point(aes(x=tSNE1, y=tSNE2, fill=cellType), alpha=0.32, shape=21, stroke=0) + scale_fill_manual(values=brewer.pal(n=6, name="Dark2"))+ theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"), legend.title=element_text(colour="black"))
```

### Pseudo-time Trajectory Analysis

```{r}
# source("https://bioconductor.org/biocLite.R")
# biocLite("destiny")
library(destiny)
data_CD4 <- data[data[, 1]=="CD4T", -c(1,2,3)]
class(data_CD4) <- "numeric"

sigmas <- find_sigmas(log2(data_CD4 + 1))
CD4.destiny <- DiffusionMap(log2(data_CD4 + 1), optimal_sigma(sigmas))
CD4.destiny.df <- data.frame(DC1=eigenvectors(CD4.destiny)[,1], DC2=eigenvectors(CD4.destiny)[,2], DC3=eigenvectors(CD4.destiny)[,3])
CD4.destiny.df <- cbind(CD4.destiny.df, log2(data_CD4 + 1))

colnames(CD4.destiny.df)

ggplot(CD4.destiny.df, aes(x=DC1, y=DC2)) + geom_point(aes(fill=`CD45RA-153(Eu153)Dd`), alpha=0.32, shape=21, stroke=0) + scale_fill_gradientn(colours=rev(brewer.pal(n=9, name="RdYlBu"))) + theme_bw() + theme(axis.title.x = element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"))

# install.packages("scatterplot3d")
library(scatterplot3d)

ii <- cut(CD4.destiny.df[, "CD45RA-153(Eu153)Dd"], breaks = seq(min(CD4.destiny.df[, "CD45RA-153(Eu153)Dd"]), max(CD4.destiny.df[, "CD45RA-153(Eu153)Dd"]), len = 1000), include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(brewer.pal(n=9, name="RdYlBu"))(1000)[ii]

scatterplot3d(x=CD4.destiny.df$DC1, y=CD4.destiny.df$DC2, z=CD4.destiny.df$DC3, angle=190, pch=21, xlab="DC1", ylab="DC2", zlab="DC3", bg=colors)
```



