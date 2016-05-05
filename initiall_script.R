# RScript to analyse the data by LLRS

library("SummarizedExperiment")
library("edgeR")
library("geneplotter")
library("plyr")
library("ggplot2")
library("scales")
# Install dplyr or something

thca <- readRDS("seTHCA.rds")
meta.info <- colData(thca)
meta.info.summary <- mcols(meta.info, use.names=TRUE)


filter.info <- function(x){
    # Function to filter the NA, [Not Available], [Not Applicable] and [Unknown]
    # Return the proportion of information on x
    nas <- sum(is.na(x))
#     , summary(x)["[Not Available]"],
#                summary(x)["[Unknown]"], summary(x)["[Not Applicable]"],
#                summary(x)["[Not Available]|[Not Available]|[Not Available]"],
#                summary(x) ["[Not Evaluated]"],
#                na.rm = TRUE)
    return(1 - (nas/length(x)))
}

na.replace <- function(x){
    # Remove the unwanted factors
    lev <- levels(x)
    lev <- lev[!(lev %in% "[Not Available]")]
    lev <- lev[!(lev %in% "[Unknown]")]
    lev <- lev[!(lev %in% "[Not Applicable]")]
    lev <- lev[!(lev %in% "[Not Available]|[Not Available]|[Not Available]")]
    lev <- lev[!(lev %in% "[Not Evaluated]")]

    x <- factor(x, levels=lev)
}

meta.info2 <- as.data.frame(lapply(meta.info, na.replace))
sample.info <- sapply(meta.info2, filter.info)

# Plotting the histogram of the information
hist(sample.info)
# Seeing the most informative to set a threshold,
# frequency are the number of columns with such % of information
hist(sample.info, xlim=c(0.3, 1), ylim=c(0, 35))

# Representing the information by column vs total information
relative.info <- sample.info[order(sample.info)]/sum(sample.info)
plot(relative.info, main="Information brought by column")
# Very few columns bring all the meaning


# Filter those variable with more than 60% not being NA
filtered <- meta.info2[, sample.info > 0.6]
meta.info.summary2 <- meta.info.summary[sample.info > 0.6,]

# We end up with 45 columns
selected.meta <- meta.info.summary[sample.info > 0.6,]

# Date
filtered$form_completion_date <- as.Date(filtered$form_completion_date)

# Numeric variables
filtered$lymph_nodes_examine <- as.numeric(filtered$lymph_nodes_examine)
filtered$lymph_nodes_examined_he_count <- as.numeric(filtered$lymph_nodes_examined_he_count)
filtered$days_to_initial_pathologic_diagnosis <- as.numeric(filtered$days_to_initial_pathologic_diagnosis)
filtered$birth_days_to <- as.numeric(filtered$birth_days_to)
filtered$last_contact_days_to <- as.numeric(filtered$last_contact_days_to)
filtered$age_at_diagnosis <- as.numeric(filtered$age_at_diagnosis)
filtered$tumor_size_width <- as.numeric(filtered$tumor_size_width)
filtered$tumor_size_width.1 <- as.numeric(filtered$tumor_size_width.1)
filtered$tumor_size_width.2 <- as.numeric(filtered$tumor_size_width.2)

# All the other variables are categorical

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
# Make lables as percentatge
p <- scale_y_continuous(labels=percent)

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
fplot + geom_bar(aes(race, fill=gender))+b+l
fplot + geom_bar(aes(x=gender, y=(..count..)/sum(..count..), fill=race))+b+l

# Seeing if the type of patient is correlated with age of diagnosis
fplot + geom_bar(aes(age_at_diagnosis, fill=type), stat="count") + b

# Date relationship with samples
fplot + geom_bar(aes(form_completion_date)) + b + l

# Explore where is the tumor
fplot + geom_bar(aes(tumor_tissue_site, fill=type)) + b +l
fplot + geom_bar(aes(tumor_tissue_site, fill=tissue_source_site)) + b +l
# Some tumors are not know for sure it is in the thyroid

# History of the thyroid and type and gender
fplot + geom_bar(aes(history_thyroid_disease, fill=type)) + b +l
fplot + geom_bar(aes(history_thyroid_disease, fill=gender)) + b +l

# To be the proportion of gender on the samples
fplot+geom_bar(aes(type, fill=gender)) + b + l

# Compare the source site of the samples
fplot+geom_bar(aes(tissue_source_site, fill=gender)) + b + l
fplot+geom_bar(aes(tissue_source_site, fill=type)) + b + l

# See how did they take care of homogenity in the tissue number extraction
fplot + geom_bar(aes(lymph_nodes_examined)) + b + l
fplot + geom_bar(aes(lymph_nodes_examined, fill=gender)) + b + l
fplot + geom_bar(aes(lymph_nodes_examined, fill=type)) + b + l
hist(as.numeric(filtered$lymph_nodes_examined), breaks=24) # Every 2 nodes

fplot + facet_wrap(~type)+geom_bar(aes(tumor_focality, fill=gender))+b+l

# Histogram of gender and tumor_focality over age
fplot + geom_bar(aes(as.numeric(age_at_diagnosis), fill=gender))+b+l
fplot +geom_bar(aes(as.numeric(age_at_diagnosis), fill=tumor_focality), na.rm=T)+ b + l

# Gender on residual tumor
fplot+geom_bar(aes(residual_tumor, fill=gender))+b

fplot+geom_bar(aes(type, fill=gender))+b

# Did he/she had  other malignancy?
fplot+geom_bar(aes(type, fill=history_other_malignancy))+b

# How was the diagnostic done
fplot+geom_bar(aes(histologic_diagnosis, fill=gender))+b+l

# Where on the thyroid was placed
fplot+geom_bar(aes(laterality, fill=type))+b+l

fplot+geom_point(aes(age_at_diagnosis, lymph_nodes_examined))+b
fplot+geom_point(aes(as.numeric(age_at_diagnosis), as.numeric(lymph_nodes_examined)))+b

# Correlating size
fplot+geom_violin(aes(tumor_size_width.1, tumor_size_width.2))+b
fplot+geom_point(aes(tumor_size_width.1, tumor_size_width.2))+geom_density2d(aes(tumor_size_width.1, tumor_size_width.2))+b
fplot+geom_violin(aes(tumor_size_width.1, tumor_size_width))+b
fplot+geom_violin(aes(tumor_size_width.2, tumor_size_width))+b

# fplot+geom_point(aes(history_radiation_exposure, age_at_diagnosis))+b+l
fplot+geom_bar(aes(age_at_diagnosis, fill=history_radiation_exposure))+b+l

fplot+geom_bar(aes(age_at_diagnosis, fill=histologic_diagnosis))+b+l

# Exploring information about the length and GC content of the transcripts
read.info <- as.data.frame(mcols(thca))
rownames(read.info) <- read.info$symbol
ggplot(read.info, aes(txgc, log10(txlen))) + geom_point()+geom_density2d() + b
ggplot(read.info, aes(txgc))+geom_density()+b
ggplot(read.info, aes(log10(txlen)))+geom_density()+b

## Normalization: CPM scaling
dge <- DGEList(counts = assays(thca)$counts, group = filtered$type,
               remove.zeros = TRUE) # Removing all the genes not expressed
logCPM <- cpm(dge, log = TRUE, prior.count = 4) # 4 for normal analysis

# Plot options
pars <- par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
# We store the options of the plotting
# Warning this takes time: analyzing 452 RNA-seq is slow
multidensity(as.list(as.data.frame(logCPM)), xlab = "log2 CPM",
            legend = NULL, main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

plotSmear(dge, lowess = TRUE)
abline(h = 0, col = "blue", lwd = 2) # Almost perfect with all the samples
dgenorm <- calcNormFactors(dge) # Normalization TMM within samples

# Too much
# plotMDS(dgenorm, col = c("red", "blue")[as.integer(dgenorm$samples$group)], cex = 0.7)
library("cqn")
# Within normalization
stopifnot(all(rownames(dge$counts) == rownames(read.info)))
# We should be able to map between number of read and symbol! we can't


# TODO: continue here
cqn.subset <- cqn(dge$counts, lengths = read.info$txlen, x = read.info$txgc,
                  sizeFactors = colSums(dge$counts))
# Now the filtering should be done