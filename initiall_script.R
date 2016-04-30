# RScript to analyse the data by LLRS

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("plyr")
library("ggplot2")
# Install dplyr or something

thca <- readRDS("seTHCA.rds")
meta.info <- colData(thca)

filter.info <- function(x){
    # Function to filter the NA, [Not Available], [Not Applicable] and [Unknown]
    # Return the proportion of information on x
    nas <- sum(is.na(x), summary(x)["[Not Available]"],
               summary(x)["[Unknown]"], summary(x)["[Not Applicable]"],
               summary(x)["[Not Available]|[Not Available]|[Not Available]"],
               summary(x) ["[Not Evaluated]"],
               na.rm = TRUE)
    return(1- nas/length(x))
}

sample.info <- sapply(meta.info, filter.info)

# Plotting the histogram of the information
hist(sample.info)
# Seeing the most informative to set a threshold,
# frequency are the number of columns with such % of information
hist(sample.info, xlim=c(0.3, 1), ylim=c(0, 35))

# Representing the information by column vs total information
relative.info <- sample.info[order(sample.info)]/sum(sample.info)
plot(relative.info, main="Information brought by column")
# Very few columns bring all the meaning

filtered <- meta.info[, sample.info > 0.6]
# We end up with 45 columns

## Normalization: CPM scaling
dge <- DGEList(counts = assays(thca)$counts, genes = mcols(thca)$symbol)
logCPM <- cpm(dge, log = TRUE, prior.count = 0.25)

# Plot options
pars <- par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
# We store the options of the plotting
# Warning this takes time: analyzing 452 RNA-seq is slow
multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM", legend = NULL, main = "",
             cex.axis = 1.2, cex.lab = 1.5, las = 1)

#TODO: Convert all the not wanted categories (Not Available, Not Evaluated... to NA)
# for a better (more informative) plotting
fplot <- ggplot(as.data.frame(filtered))

## ggplot options for axes and background color
# Options of the labels
l <- theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 1))
# Removing the background
b <-  theme(axis.line = element_line(colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey", size = 0.5, linetype = 2),
            panel.grid.minor.y = element_line(colour = "grey", size = 0.25, linetype = 3),
            panel.border = element_blank(),
            panel.background = element_blank())

# Comparing control vs cases and gender
fplot + geom_bar(aes(type, fill=gender)) + b

# Comparing type of tumor and tumor status
fplot + geom_bar(aes(type, fill=tumor_status)) + b

# Comparing type, and gender with ethnicity
fplot + facet_wrap(~ type) + geom_bar(aes(gender, fill=ethnicity)) + b
fplot + facet_wrap(~ gender) + geom_bar(aes(type, fill=ethnicity)) + b

# Compare age to gender and tumor
fplot + geom_bar(aes(age_at_diagnosis, fill=gender), stat="count") + b

# Compare race and ethnicity
fplot + facet_wrap(~ ethnicity) + geom_bar(aes(race)) + b + l

# Compare race and gender

# Seeing if the type of patient is correlated with age of diagnosis
fplot + geom_bar(aes(age_at_diagnosis, fill=type), stat="count") + b

# Date relationship with samples
fplot + geom_bar(aes(form_completion_date)) + b + l

# Explore where is the tumor
fplot + geom_bar(aes(tumor_tissue_site)) + b +l
# Some tumors are not know for sure it is in the thyroid

# History of the thyroid and type and gender
fplot + geom_bar(aes(history_thyroid_disease, fill=type)) + b +l
fplot + geom_bar(aes(history_thyroid_disease, fill=gender)) + b +l

# To be the proportion of gender on the samples
fplot+geom_bar(aes(type, fill=gender)) + b + l

# Compare the source site of the samples
fplot+geom_bar(aes(tissue_source_site, fill=gender)) + b + l
fplot+geom_bar(aes(tissue_source_site, fill=type)) + b + l

# compare
fplot + geom_bar(aes(lymph_nodes_examined)) + b + l
fplot + facet_wrap(~type)+geom_bar(aes(tumor_focality, fill=gender))+b+l

fplot + facet_wrap(~tumor_focality)+geom_bar(aes(age_at_diagnosis, fill=gender))+b+l
fplot +geom_bar(aes(age_at_diagnosis, fill=tumor_focality), na.rm=T)+ b + l
fplot+geom_bar(aes(gender, fill=residual_tumor))+b
fplot+geom_bar(aes(type, fill=gender))+b
fplot+geom_bar(aes(type, fill=history_other_malignancy))
