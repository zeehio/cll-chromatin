#!/usr/bin/env python

import yaml
import os
from pipelines import Project, ATACseqSample
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


sns.set_style("whitegrid")


def normalizeByIntervalLength(series):
    """
    Divides a pandas.Series with numerical values by the length of the interval
    (defined in "chr", "start", "end") in a bam file (RPK).
    """
    length = float(series["end"] - series["start"])

    rpK = (series.drop(["chrom", "start", "end"]) / length) * 1e3

    return series[["chrom", "start", "end"]].append(rpK)


def normalizeByLibrarySize(series, samples):
    """
    Divides a pandas.Series with numerical values by the total number of mapped
    reads in a bam file.
    """
    # get sample which this series belongs to
    sample = [sample for sample in samples if sample.name == series.name][0]

    # get number of aligned reads for that sample
    size = float(pysam.AlignmentFile(sample.filtered).mapped)

    return (series / size) * 1e6


def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
dataDir = os.path.join(config["paths"]["parent"], config["projectname"], "data")
resultsDir = os.path.join(config["paths"]["parent"], config["projectname"], "results")
plotsDir = os.path.join(resultsDir, "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("../metadata/sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]


# GET CONSENSUS SITES ACROSS SAMPLES
for i, sample in enumerate(samples):
    # Get summits of peaks and window around them
    peaks = pybedtools.BedTool(sample.peaks) # .slop(g=pybedtools.chromsizes('hg19'), b=250)
    # Merge overlaping peaks within a sample
    peaks = peaks.merge()
    if i == 0:
        sites = peaks
    else:
        # Concatenate all peaks
        sites = sites.cat(peaks)

# Merge overlaping peaks across samples
sites = sites.merge()
sites.saveas(os.path.join(dataDir, "all_sample_peaks.concatenated.bed"))

# Loop at summary statistics:
# interval lengths
sns.distplot([interval.length for interval in sites], bins=300, kde=False)
plt.xlabel("peak width (bp)")
plt.ylabel("frequency")
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.lengths.pdf"), bbox_inches="tight")

# calculate support (number of samples overlaping each merged peak)
for i, sample in enumerate(samples):
    if i ==0:
        support = sites.intersect(sample.peaks, wa=True, c=True)
    else:
        support = support.intersect(sample.peaks, wa=True, c=True)
support = support.to_dataframe()
support = support.reset_index()
support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
# divide sum (of unique overlaps) by total to get support value between 0 and 1
support["support"] = support[range(len(samples))].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(samples)), axis=1)
# save
support.to_csv(os.path.join(dataDir, "all_sample_peaks.support.csv"), index=False)
# plot
sns.distplot(support["support"], bins=10)
plt.ylabel("frequency")
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.support.pdf"), bbox_inches="tight")


# Remove some regions chrM, chr.*random
# After inspecting known genotypes, consider excluding more regions/chromosomes across all samples
# e.g. chr11q, chr12


# MEASURE CHROMATIN OPENNESS PER SITE PER PATIENT
# Get read counts for each sample in the union of the sites
coverage = sites.multi_bam_coverage(bams=[sample.filtered for sample in samples], p=True, q=30, D=True)
# make dataframe
coverage = coverage.to_dataframe().reset_index()
coverage.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
coverage.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.raw_coverage.bed"), sep="\t", index=False)

# Normalize by feature length (Reads per kilobase)
rpk = coverage.apply(normalizeByIntervalLength, axis=1)
# Normalize by library size - mapped reads (Reads per kilobase per million)
rpkm = rpk[[sample.name for sample in samples]].apply(normalizeByLibrarySize, args=(samples, ), axis=0)

# Save
rpkm = pd.concat([rpk[["chrom", "start", "end"]], rpkm], axis=1)
rpkm.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.rpkm.bed"), sep="\t", index=False)

# Log2 transform
rpkm[samples] = np.log2(1 + rpkm[[sample.name for sample in samples]])

# Plot
rpkmMelted = pd.melt(rpkm, id_vars=["chrom", "start", "end"], var_name="sample", value_name="rpkm")

# rpkm density
g = sns.FacetGrid(rpkmMelted, col="sample", aspect=2, col_wrap=4)
g.map(sns.distplot, "rpkm", hist=False);
plt.savefig(os.path.join(plotsDir, "fpkm_per_sample.distplot.pdf"), bbox_inches="tight")

# boxplot rpkm per sample
# Plot the orbital period with horizontal boxes
sns.boxplot(x="rpkm", y="sample", data=rpkmMelted)
plt.savefig(os.path.join(plotsDir, "fpkm_per_sample.boxplot.pdf"), bbox_inches="tight")

# pairwise rpkm scatter plot between samples
fig = plt.figure(figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')
g = sns.PairGrid(rpkm[[sample.name for sample in samples]])
g.map(hexbin)
g.fig.subplots_adjust(wspace=.02, hspace=.02);
plt.savefig(os.path.join(plotsDir, "fpkm_per_sample.pairwise_hexbin.pdf"), bbox_inches="tight")

# qv2 vs mean rpkm
rpkm['qv2'] = rpkm[[sample.name for sample in samples]].apply(lambda x: (np.std(x) / np.mean(x)) ** 2, axis=1)
rpkm['mean'] = rpkm[[sample.name for sample in samples]].apply(lambda x: np.mean(x), axis=1)

sns.jointplot('mean', "qv2", data=rpkm)
plt.xlabel("mean log2(1 + rpkm)")
plt.savefig(os.path.join(plotsDir, "fpkm_per_sample.qv2_vs_mean.pdf"), bbox_inches="tight")


# GLOBAL CHARACTERIZATION
# Overview:
# Correlation
sns.clustermap(rpkm[samples].corr(), square=True)
plt.savefig(os.path.join(dataDir, "all_sample_peaks.concatenated.correlation_clustering.pdf"), bbox_inches="tight")

# Heatmap sites vs patients
# hierarchical clustering

# PCA
pca = PCA()
X = pca.fit_transform(rpkm[[sample.name for sample in samples]])

# get unique colors
from matplotlib.pyplot import cm
color = iter(cm.rainbow(np.linspace(0,1,len(samples))))

# plot each point
for i, sample in enumerate(samples):
    plt.scatter(X[i, 0], X[X[i, 0] == i, 1], label=sample.name, color=color.next())
plt.legend()

# MDS


# INTER-SAMPLE VARIABILITY ANALYSIS
# Get most variable sites
# plot CV2 vs mean openness

# Enrichment of variable sites:
# LOLA
# get closest gene: GO, KEGG, OMIM, mSigDB
# de novo motif finding - enrichment

# Try calling differentially open sites (DOS)
# with DESeq2


# CLASSIFICATION
# Stratify patients on:
# - treated vs untreated
# - ...
# Train classifiers on groups
# Predict for all samples
# Assess: ROC, AUC
