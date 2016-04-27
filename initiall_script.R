# RScript to analyse the data by LLRS

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")

thca <- readRDS("seTHCA.rds")
meta.info <- colData(thca)

filter.info <- function(x){
    # Function to filter the
    if (sum(is.na(x))/length(x) > 0.5){
        FALSE
    } else {
        TRUE
    }
}
sample.info <- sapply(meta.info, filter.info)
filtered <- meta.info[, sample.info]

## Normalization: CPM scaling
dge <- DGEList(counts = assays(thca)$counts, genes = mcols(thca)$symbol)
logCPM <- cpm(dge, log = TRUE, prior.count = 0.25)

# Plot options
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "",
             cex.axis = 1.2, cex.lab = 1.5, las = 1)
