#set working directory make sure it is not linked to onedrive
setwd("D:/New")
getwd()

# Set Bioconductor repositories
options(repos = BiocManager::repositories())
install.packages("httr2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
clinical_data <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")
colnames(clinical_data) # Displays the available clinical attributes
#Check the structure of the dataset to identify the column of interest
str(clinical_data)
#ajcc_pathologic_m ,  prior_malignancy

table(clinical_data$ajcc_pathologic_m)
table(clinical_data$prior_malignancy)
# counts data
query <- GDCquery(
  project = "TCGA-PAAD",                      # Target the TCGA-PAAD project
  data.category = "Transcriptome Profiling",   # Category for transcriptomic data
  experimental.strategy = "RNA-Seq",           # Experimental strategy for RNA-Seq
  access = "open"                             # Access level: open data (no restricted access)
)

GDCdownload(query)
expression_data <- GDCprepare(query)
paad_matrix <- assay(expression_data, 'unstranded')
write.csv(expression_data, "transcriptomic_data.csv", row.names = TRUE)

# i am now setting the barcodes as they should match in clinical and the counts df
# gonna compare the column in query results and clinical 
# Inspect the first few entries
head(query$results[[1]]$cases.submitter_id)
head(clinical_data$bcr_patient_barcode)
cases <- toupper(query$results[[1]]$cases.submitter_id)
barcodes <- toupper(clinical_data$bcr_patient_barcode)

# Trim any leading/trailing whitespaces
cases <- trimws(cases)
barcodes <- trimws(barcodes)

# Find matching barcodes
matching_barcodes <- cases[cases %in% barcodes]

# Show the matches
print(matching_barcodes)

# Check if there are unmatched cases
unmatched_cases <- cases[!cases %in% barcodes]
print(unmatched_cases)


# Extract cases and bcr_patient_barcode from the query results
query_cases <- data.frame(
  bcr_patient_barcode = toupper(query$results[[1]]$cases.submitter_id),
  cases = query$results[[1]]$cases # Adjust column reference if necessary
)

# Ensure both dataframes have consistent formats
clinical_data$bcr_patient_barcode <- toupper(clinical_data$bcr_patient_barcode)

# Merge the dataframes based on the bcr_patient_barcode
clinical_data <- merge(
  clinical_data, 
  query_cases, 
  by = "bcr_patient_barcode", 
  all.x = TRUE # Keep all rows from clinical_df
)

# View the updated clinical dataframe
head(clinical_data)
# Ensure the clinical data is ordered by sample (cases)
clinical_filtered <- clinical_data[match(colnames(paad_matrix), clinical_data$cases), ]

# Confirm alignment
all(colnames(paad_matrix) == clinical_filtered$cases)  # Should return TRUE
# the number of observations were decreased to 183 in clinical filtered
cases <- toupper(clinical_filtered$cases)
barcodes <- toupper(clinical_data$cases)
unmatched_cases <- cases[!cases %in% barcodes]
unique_to_clinical_data <- setdiff(barcodes, cases)

# so now i found that there much be repetition in clinical_data

# Convert the cases column to uppercase for consistency
clinical_data$cases <- toupper(clinical_data$cases)

# Identify duplicate entries in the cases column
duplicated_cases <- clinical_data$cases[duplicated(clinical_data$cases)]

# Print the repetitive entries
cat("Repetitive entries in the 'cases' column:\n")
print(duplicated_cases)
# i thnk there was a missing or NA entries in "cases" so i compared the metastatic and primary cases so now
# M0 M1 MX (Before) 
# 86  6 98 
# M0 M1 MX (After)
# 81  5 97
# as for prior malignancy (before) [no=170 , yes=20] | (after) [no=164 , yes=19]
table(clinical_filtered$ajcc_pathologic_m)
table(clinical_filtered$prior_malignancy)

clinical_combined_filtered <- as.data.frame(clinical_filtered)  # Ensure it's a data frame
rownames(clinical_combined_filtered) <- clinical_combined_filtered[[71]]  # Set the first column as row names

# trouble shooting clinical data
all(colnames(paad_matrix) %in% rownames(clinical_combined_filtered)) # should return true check if the names match
dim(clinical_combined_filtered)
is.factor(clinical_combined_filtered$prior_malignancy)  # Should return TRUE
levels(clinical_combined_filtered$prior_malignancy)
#If Condition is not a factor or has only one level, adjust it:
clinical_combined_filtered$prior_malignancy <- as.factor(clinical_combined_filtered$prior_malignancy)


# trouble shooting counts matrix
class(paad_matrix)  # Should return "matrix" and "integer"
dim(paad_matrix)
is.numeric(paad_matrix)  # Should return TRUE
sum(is.na(paad_matrix))  # Should return 0 if there are no NA values
range(paad_matrix, na.rm = TRUE)  # Should not include NA
typeof(paad_matrix) # Should be "integer"

#trouble shooting  if it was not numeric or had "NA" in the matrix (it didn't have that so skipped it)
paad_matrix <- matrix(as.integer(paad_matrix), nrow = nrow(paad_matrix), 
                        dimnames = dimnames(paad_matrix))  # Convert to integer


paad_matrix[is.na(paad_matrix)] <- 0  # Replace NAs with 0

#deseq2 analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = paad_matrix,
  colData = clinical_combined_filtered,
  design = ~ prior_malignancy
)
vsd <- vst(dds, blind = FALSE)  # Using blind=FALSE ensures that the transformation takes into account the experimental design
class(vsd)
plotPCA(vsd, intgroup = "prior_malignancy")

summary(dds)  # Overview of the data set
dds <- DESeq(dds)
install.packages("gplots")
install.packages("ggplot2")

library(gplots)
library(ggplot2)

heatmap.2(assay(vsd), trace="none", col=bluered(100))

res <- results(dds)
res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

volcano <- ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  theme_minimal() + labs(x="log2 Fold Change", y="-log10(p-value)")
print(volcano)

plotMA(res, main="DESeq2 MA Plot", ylim=c(-5,5))

write.csv(as.data.frame(res), "deseq2_resTcga.csv")
significant_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
selected_genes <- rownames(significant_genes)

norm_counts <- counts(dds, normalized = TRUE)
gene_data <- norm_counts[selected_genes, ]

svm_input <- t(gene_data)

labels <- factor(clinical_combined_filtered$prior_malignancy)  # Replace `metadata` with your sample metadata dataframe


scaled_data <- scale(svm_input)

set.seed(123)  # For reproducibility
train_indices <- sample(1:nrow(scaled_data), 0.8 * nrow(scaled_data))
train_data <- scaled_data[train_indices, ]
train_labels <- labels[train_indices]

test_data <- scaled_data[-train_indices, ]
test_labels <- labels[-train_indices]


install.packages("e1071")
library(e1071)

# Train the model with probability = TRUE
svm_model <- svm(train_data, train_labels, kernel = "linear", cost = 1, scale = FALSE, probability = TRUE)

predictions <- predict(svm_model, test_data)
# Get predicted probabilities for the test data
svm_probs <- predict(svm_model, test_data, probability = TRUE)

# svm_probs will now contain the predicted probabilities for the classes
# Access the predicted probabilities
probabilities <- attr(svm_probs, "probabilities")

library(caret)
confusionMatrix(predictions, test_labels)

important_genes <- colnames(train_data)[order(abs(svm_model$coefs), decreasing = TRUE)]
print(important_genes)


#visualize results
# Perform PCA on scaled data
pca <- prcomp(scaled_data)

# Plot the first two principal components
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = labels)
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data", x = "PC1", y = "PC2") +
  scale_color_manual(values = c("blue", "red"))


install.packages("Rtsne")
library(Rtsne)


# Perform t-SNE
tsne <- Rtsne(scaled_data, dims = 2, pca = TRUE)

# Create a data frame for the plot
tsne_df <- data.frame(TSNE1 = tsne$Y[, 1], TSNE2 = tsne$Y[, 2], Condition = labels)

# Plot the t-SNE results
ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = Condition)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE of Gene Expression Data", x = "t-SNE1", y = "t-SNE2") +
  scale_color_manual(values = c("yellow", "green"))

#roc curve
install.packages("pROC")
library(pROC)

# For binary classification, use the class probabilities for the positive class
roc_curve <- roc(test_labels, probabilities[, 2])  # Assuming the second column corresponds to the positive class
plot(roc_curve, main = "ROC Curve", col = "blue")


library(caret)

# Confusion Matrix
confusion_matrix <- confusionMatrix(predictions, test_labels)
print(confusion_matrix)

library(caret)
fourfoldplot(confusion_matrix$table, color = c("lightblue", "red"))

#feature importance
# Check the dimensions of the coefficients
length(svm_model$coefs)  # Number of coefficients
# Get the names of the genes used in the model
used_genes <- colnames(train_data)[svm_model$index]

# Create a data frame for feature importance
importance_df <- data.frame(Gene = used_genes, Importance = as.vector(svm_model$coefs))

# Sort by importance
importance_df <- importance_df[order(abs(importance_df$Importance), decreasing = TRUE), ]

# Plot the top 10 most important genes
library(ggplot2)
ggplot(importance_df[1:10, ], aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Important Genes", x = "Gene", y = "Importance")

print(used_genes)
