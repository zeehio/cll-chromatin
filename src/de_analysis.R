library(data.table)
library(DESeq2)


data_dir = "data_submission"
plots_dir = file.path("results", "plots")

# read in
df = read.table(file.path(data_dir, "cll_expression_matrix.log2.csv"), sep=",", header=TRUE)
# un log2
countData = (2 ** df[, 2:(ncol(df) - 1)]) - 1
# add genes as index
rownames(countData) <- df$ensembl_gene_id
# replace column names with patient id
colnames(countData) = t(as.data.table(lapply(
    colnames(countData),
    function(x) unlist(strsplit(gsub("_hg19", "", gsub("CLL_RNA.seq_", "", x)), "_")[[1]][1]))))
# coherce to integer
countData <- apply(countData, c(1, 2), function (x) {
    (as.integer(x))
})

# get IGHV mutation status of selected samples for RNA-seq
sel = read.table(file.path("metadata", "selected_samples.tsv"), header=TRUE)
colData = sel[sel$patient_id %in% colnames(countData), ]
rownames(colData) = colData$patient_id
colData <- colData[order(as.numeric(row.names(colData))), ]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ sample_cluster)
dds <- DESeq(dds)

alpha = 0.05

for (contrast in list(
    c("sample_cluster", "uCLL", "iCLL"),
    c("sample_cluster", "iCLL", "mCLL"),
    c("sample_cluster", "uCLL", "mCLL"))) {
    # name for condition
    condition_name = paste0(contrast[2:3], collapse="-")

    # get results
    res <- results(dds, contrast=contrast, alpha=alpha)
    res <- res[order(res$padj), ]

    # Plot MA
    pdf(file.path(plots_dir, paste0("gene_expression.differential.", condition_name, "_ma.pdf")))
    plotMA(res, main=paste0("MA plot: ", condition_name), ylim=c(-7, 7))
    dev.off()

    # Export all
    write.csv(as.data.frame(res), file=file.path(data_dir, paste0("gene_expression.differential.", condition_name, "_results.csv")))

    # Export DEGs
    resSig <- subset(res, padj < alpha)
    write.csv(as.data.frame(resSig), file=file.path(data_dir, paste0("gene_expression.differential.", condition_name, "_results.degs.csv")))
}

# transform values
rld <- rlog(dds)

# Plot PCA
pdf(file.path(plots_dir, paste0("gene_expression.differential.pca.pdf")))
plotPCA(rld, intgroup=c("sample_cluster"))
dev.off()

# distances between samples
sampleDists <- dist(t(assay(rld)))
sampleCor <- corr(t(assay(rld)))

# Plot heatmaps
library(gplots)
pdf(file.path(plots_dir, paste0("gene_expression.distance.pdf")))
heatmap.2(as.matrix(sampleDists, trace="none"))
dev.off()
pdf(file.path(plots_dir, paste0("gene_expression.correlation.pdf")))
heatmap.2(as.matrix(sampleCor, trace="none"))
dev.off()
