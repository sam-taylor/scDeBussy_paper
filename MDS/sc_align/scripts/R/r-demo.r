# Basic analysis demo
library(ggplot2)
library(httpgd)

hgd() # Start graphics device

# Simple plot
ggplot(mtcars, aes(mpg, wt)) +
    geom_point(color = "steelblue") +
    ggtitle("R Plotting in VSCode via WSL")

# Verify package management
if (!requireNamespace("dplyr")) install.packages("dplyr")
library(dplyr)
mtcars %>% summarize(avg_mpg = mean(mpg))

# In R console - update repos first
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install base dependencies
install.packages(c("textshaping", "ragg", "pkgdown", "xml2", "fs"))

# If pkgdown fails, install via remotes
if (!require("pkgdown")) remotes::install_github("r-lib/pkgdown")

# Install Bioconductor dependencies
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("BiocGenerics", "S4Vectors"))

## make sure you remove the old version of TSCAN if you've downloaded it before
if ("TSCAN" %in% rownames(installed.packages())) remove.packages("TSCAN")
## download dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    BiocManager::install("ComplexHeatmap")
}
if (!require("devtools")) {
    install.packages("devtools")
}
## downlaod the most updated TSCAN from Github
devtools::install_github("zji90/TSCAN")
## downlaod the most updated Lamian from Github
devtools::install_github("Winnie09/Lamian")
