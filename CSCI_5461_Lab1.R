## Load all necessary packages
library(dplyr)

## 1. Load all necessary data
count_data <- read.csv(file = "HW1-GSE62944-count.csv", header = TRUE)
clinical_data <- read.csv(file = "HW1-GSE62944-clinical.csv", header = TRUE)
DESeq <- read.csv(file = "HW1-DESeq2.csv", header = TRUE)

###############################################################################

## 2. Explore the datasets

## a) How many genes are included in the dataset?
# Calculate how many genes in the datasets:
nrow(count_data)
nrow(DESeq)

# Calculate how many patient samples in the dataset:
nrow(clinical_data)

# Calculate how many are short- and long-term:
unique(clinical_data$Group)
sum(clinical_data$Group=="short")
sum(clinical_data$Group=="long")

###############################################################################

## 3. Data processing and normalization

## a) Compute the total number of reads per sample:

read_per_sample <- colSums(count_data[, c(2:ncol(count_data))])

# Plot a bar graph that summarizes these values:

pdf("Bar_plot_read_counts.pdf")
barplot(read_per_sample, xaxt = "n", xlab = "Samples", ylab = "Read counts", main = "Read counts in all samples", ylim = c(0, 8e+07))
dev.off()

# b) Perform a "total count" normalization for read-depth in each sample:

## **Define a Normalization Function**
normalize_function <- function(read_count_file) {
  # Calculate total number of reads for each sample:
  read_count_per_sample <- colSums(read_count_file)
  # Find the median number of total read count across all samples:
  N_med <- median(read_count_per_sample)
  # Generate a empty dataframe with the same dimension as the original data: 
  norm_data <- data.frame(matrix(0, 
                                 nrow = nrow(read_count_file), 
                                 ncol = ncol(read_count_file)))
  # Use each column as a unit to normalize the data.
  for (i in 1:ncol(read_count_file)) {
    norm_data[i] = read_count_file[i] * (N_med/read_count_per_sample[i])
  }
  # Set the row and column's names
  rownames(norm_data) = rownames(read_count_file)
  colnames(norm_data) = colnames(read_count_file)
  
  return (norm_data)
}

# Clean up the original count dataset
cleaned_count_file <- count_data[2:ncol(count_data)]
rownames(cleaned_count_file) <- count_data$X

# Use self-defined function to normalize the cleaned-up read count dataset
normalized_dataset <- normalize_function(cleaned_count_file)

# Compute the total number of reads per sample after normalization.
norm_read_per_sample <- colSums(normalized_dataset)

# Plot a bar graph that summarizes these values:
pdf("bar_plot_norm.pdf")
barplot(norm_read_per_sample, xaxt = "n", xlab = "Samples", ylab = "Normalized read counts", main = "Reads counts in all samples after normalization")
dev.off()

# c) Log-transform the normalized count data.

log_norm_data <- log(normalized_dataset+1)

## Plot the histogram.
pdf("histogram_norm_log.pdf")
hist(as.matrix(log_norm_data), breaks = 100, xlab = "log-transformed expression levels", main = "Histogram of log-transformed expression levels in all samples")
dev.off()

## d) Plot individual histograms of the log-transformed data for the first five samples.
## Plot the histogram.
pdf("histogram_norm_log_all.pdf")
par(mfrow = c(3,2))
hist(log_norm_data[, 1], breaks = 100, xlab = "log-transformed expression levels", xlim = c(0,14), main = "Histogram of log-transformed\n expression levels in sample 1")

hist(log_norm_data[, 2], breaks = 100, xlab = "log-transformed expression levels", xlim = c(0,14), main = "Histogram of log-transformed\n expression levels in sample 2")

hist(log_norm_data[, 3], breaks = 100, xlab = "log-transformed expression levels", xlim = c(0,14), main = "Histogram of log-transformed\n expression levels in sample 3")

hist(log_norm_data[, 4], breaks = 100, xlab = "log-transformed expression levels", xlim = c(0,14), main = "Histogram of log-transformed\n expression levels in sample 4")

hist(log_norm_data[, 5], breaks = 100, xlab = "log-transformed expression levels", xlim = c(0,14), main = "Histogram of log-transformed\n expression levels in sample 5")
dev.off()

# (e) Implement and perform quantile normalization on your log-transformed data across all samples such that each has the same empirical distribution. Plot a histogram of the quantile-normalized data for each of the first five samples.

# **Define a Quantile Normalization Function**

quan_normalize_function <- function(data_need_quan_normalize) {
  # Get the rank of original data so that the calculated ranked average data can be put back:
  rank_of_data <- data.frame(apply(data_need_quan_normalize, 2, function(x) rank(x, ties.method = "average")))
  
  # Order the data in each column. Return a dataframe, each column is ordered:
  sorted_data <- data.frame(apply(data_need_quan_normalize, 2, sort))
  
  # Calculate the average of each level of ordered data, i.e., calculate the average of each row of newly ordered data:
  average <- apply(sorted_data, 1, mean)
  
  # Based on the rank of original dataset, put the average of ordered data back: 
  quan_normalized <- data.frame(apply(rank_of_data, 2, function(x) average[x]))
  
  # Set the row names:
  rownames(quan_normalized) <- rownames(data_need_quan_normalize)
  
  return(quan_normalized)
}

# Use the self-defined function to quantile-normalize the log-transformed data:
quan_norm <- quan_normalize_function(log_norm_data)

# Compare the results from self-defined quantile normalization function (quan_norm) to the built-in function from "preprocessCore":

# Install and load the package:
## BiocManager::install("preprocessCore")
library(preprocessCore)

# Convert the log-transformed data into matrix:
mat <- as.matrix(log_norm_data)
# Quantile normalize data.
quan_norm_from_library <- data.frame(normalize.quantiles(mat))

## Conclusion: The results from the self-defined function (quan_norm) and "normalize.quantiles" function (quan_norm_from_library) are exact the same. 

# Plot a histogram of the quantile-normalized data for each of the first five samples:
pdf("histogram_norm_log_quan.pdf")
par(mfrow=c(3,2))
hist(quan_norm[, 1], breaks = 100, xlab = "quantile-normalized expression levels", main = "Histogram of quantile-normalized\n expression levels in sample 1")

hist(quan_norm[, 2],breaks = 100,  xlab = "quantile-normalized expression levels", main = "Histogram of quantile-normalized\n expression levels in sample 2")

hist(quan_norm[, 3], breaks = 100, xlab = "quantile-normalized expression levels", main = "Histogram of quantile-normalized\n expression levels in sample 3")

hist(quan_norm[, 4], breaks = 100, xlab = "quantile-normalized expression levels", main = "Histogram of quantile-normalized\n expression levels in sample 4")

hist(quan_norm[, 5], breaks = 100, xlab = "quantile-normalized expression levels", main = "Histogram of quantile-normalized\n expression levels in sample 5")
dev.off()

###############################################################################

## 4. Analysis of differential expression

## a) Find the top 10 gene.

# Separate the quantile-normalized count data based on "short" and "long" term patients.

# Transpose the quantile-normalized data:
quan_norm_transpose <- data.frame(t(quan_norm))

# Change the sample name in clinical data and make them as row names, so they can be aligned with what the row names of transposed quantile-normalized data:
clinical_data$sampleName <- gsub("-", ".", clinical_data$sampleName)
rownames(clinical_data) <- clinical_data$sampleName

# Based on the row names, combine both clinical data and quantile-normalized and transposed data together:
merge_data <- merge(clinical_data, quan_norm_transpose, by = "row.names", all = TRUE)

# Change the row names of merged data into all gene names:
rownames(merge_data) <- merge_data$Row.names

# Subtract the first 3 columns, since they are not useful in further analysis:
cleaned_data <- merge_data[, -c(1:3)]

# Perform the Wilcoxon test to compare the differential expression of each gene between short- and long-termed patients. All the p-values will be stored in a dataframe.
wilcox_result <- data.frame(apply(cleaned_data[2:(ncol(cleaned_data)-1)], 2, function(x) wilcox.test(x ~ cleaned_data$Group)[["p.value"]]))

# Rename the column name of the dataframe:
colnames(wilcox_result) <- "p_value"

# Remove all NA values in the wilcoxon test result:
wilcox_result <- na.omit(wilcox_result)

# Sort the data in an increasing order.
ordered_wilcox_result <- arrange(wilcox_result, p_value)

# Get the top 10 genes that have the lowest p-values.
top10_genes <- ordered_wilcox_result[1:10, , drop = FALSE]
print(top10_genes)

## b) Report the number of significant genes based on the Wilcoxon rank-sum test at an uncorrected p < 0.05 cutoff. Report the number of significant genes from the DESeq2 approach at an uncorrected p < 0.05.

# The number of significant genes based on the Wilcoxon rank-sum test is: 
sum(wilcox_result$p_value < 0.05, na.rm = TRUE)

# The number of significant genes based on provided DESeq2 report is:
sum(DESeq$pvalue<0.05, na.rm = TRUE)

## c) What is the overlap between the set of genes deemed significant at an uncorrected p < 0.05 cutoff by the two different approaches (DESeq2 and Wilcoxon rank-sum statistics)? Report the total number of genes that overlap and a list of those overlapping genes.

# Significant genes based on the Wilcoxon rank-sum test:
sig_p_wilcoxon <- wilcox_result[wilcox_result$p_value < 0.05, , drop = FALSE]
sig_gene_wilcoxon <- rownames(sig_p_wilcoxon)

# Significant genes based on the provided DESeq2 file:
sig_p_deseq <- DESeq[DESeq$pvalue < 0.05, ]
sig_gene_DESeq <- sig_p_deseq$X

# Get the overlap genes:
overlap_gene <- intersect(sig_gene_wilcoxon, sig_gene_DESeq)
# The number of overlap genes is:
length(overlap_gene)
# Export the overlap genes into a txt file:
write.table(overlap_gene, "overlap_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

###############################################################################

## 5. Multiple hypothesis correction

## a) Use Bonferroni correction with the DESeq2 results to identify differentially expressed genes at a global significance level (corrected p < 0.05) and report the number of significant genes.

# Calculate adjusted p-value based on Bonferroni correction:
DESeq$adjusted_p_bonf <- p.adjust(DESeq$pvalue, method = "bonferroni")

# Find the significant genes based on adjusted p values:
sig_adjusted_p_bonf_deseq <- DESeq[DESeq$adjusted_p_bonf < 0.05, ]

# The number of significant genes after Bonferroni correction with DESeq2 results is:
nrow(sig_adjusted_p_bonf_deseq)

## b) Use the Benjamini-Hochberg step-up procedure to control the False Discovery Rate (FDR) with the DESeq2 statistics to identify differentially expressed genes at an FDR < 0.05.

# Calculate adjusted p-value based on Benjamini-Hochberg procedure:
DESeq$adjusted_p_BH <- p.adjust(DESeq$pvalue, method = "fdr")

# Find the significant genes based on adjusted p values:
sig_adjusted_p_BH_deseq <- DESeq[DESeq$adjusted_p_BH < 0.05, ]

# The number of significant genes after Benjamini-Hochberg procedure with DESeq2 results is:
nrow(sig_adjusted_p_BH_deseq)

## c) Rank the genes by their p-values in ascending order and plot the p-values and the threshold used for adjusted p-values as a function of the index of the ranked genes. 

# Order DESeq data based on the ascending order of original p-value:
ordered_deseq <- DESeq[order(DESeq$pvalue),]

# Calculate the adjusted p-value threshold based on the ordered p-value:
ordered_deseq$threshold <- 0.05 * c(1:nrow(ordered_deseq)) / nrow(ordered_deseq)

# Make the plots:
pdf("pvalue_vs_threshold.pdf", width = 8, height = 6)
plot(x=c(1:500), y=ordered_deseq$pvalue[1:500], xlab = "Indices", ylab = "p-value", col = "orange", type = "l", lwd = 2)
lines(x=c(1:500), y=ordered_deseq$threshold[1:500], type = "l", lty = 2, col = "purple", lwd = 2)
dev.off()

###############################################################################

## Supplementary analysis for 4a.

# Provide the summary statistics of KCTD3 gene expression level based on shor- and long-term patients:
library(psych)
describeBy(cleaned_data$KCTD3, cleaned_data$Group)
describeBy(cleaned_data$PRPF40A, cleaned_data$Group)
