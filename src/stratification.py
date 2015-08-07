#!/usr/bin/env python

"""
This script makes plots to ilustrate the stratification between CLL patients.
"""

import yaml
import os
from pipelines.models import Project, ATACseqSample
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing
import parmap
import pysam
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.decomposition import RandomizedPCA
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
import itertools


sns.set_style("whitegrid")
sns.set_context("paper")


def count_reads_in_intervals(bam, intervals):
    """
    """
    counts = dict()

    bam = pysam.Samfile(bam, 'rb')

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

    for interval in intervals:
        if interval.split(":")[0] not in chroms:
            continue
        counts[interval] = bam.count(region=interval)
    bam.close()

    return counts


def normalize_by_interval_length(series):
    """
    Divides a pandas.Series with numerical values by the length of the interval
    (defined in "chr", "start", "end") in a bam file (RPK).
    """
    length = float(series["end"] - series["start"])

    rpk = (series.drop(["chrom", "start", "end"]) / length) * 1e3

    return series[["chrom", "start", "end"]].append(rpk)


def normalize_by_library_size(series, samples, rm_mt=True):
    """
    Divides a pandas.Series with numerical values by the total number of mapped
    reads in a bam file.
    """
    # get sample which this series belongs to
    sample = [sample for sample in samples if sample.name == series.name][0]

    # get number of aligned reads for that sample
    size = float(pysam.AlignmentFile(sample.filtered).mapped)

    if rm_mt:
        mt = pysam.AlignmentFile(sample.filtered).count(reference="chrM")
        return (series / (size - mt)) * 1e6
    else:
        return (series / size) * 1e6


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_sample_id(name):
    return name.split("_")[3:4][0]


def annotate_igvh_mutations(samples):
    new_samples = list()

    for sample in samples:
        _id = name_to_sample_id(sample.name)
        if clinical.loc[clinical['SampleID'] == _id, 'IGVH mut/unmut'].tolist()[0] == 1:
            sample.mutated = True
        elif clinical.loc[clinical['SampleID'] == _id, 'IGVH mut/unmut'].tolist()[0] == 2:
            sample.mutated = False
        else:
            sample.mutated = None
        new_samples.append(sample)
    return new_samples


def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)


def get_cluster_classes(den, label='ivl'):
    from collections import defaultdict

    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][j] for j in l]
        cluster_classes[c] = i_l

    return cluster_classes


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
dataDir = os.path.join(config["paths"]["parent"], config["projectname"], "data")
resultsDir = os.path.join(config["paths"]["parent"], config["projectname"], "results")
plotsDir = os.path.join(resultsDir, "plots")

# Get clinical info
clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))
clinical = clinical[["PatientID", "SampleID", "timepoint", "IGVH mut/unmut", "gender", "Via (%)", "dob", "Date of Diagnosis"]].drop_duplicates()


# Start project
# prj = pickle.load(open("prj.pickle", 'rb'))
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample if s.cellLine == "CLL"]

# Annotate with igvh mutated status
# add "mutated" attribute to sample depending on IGVH mutation status
samples = annotate_igvh_mutations(samples)


# GET CONSENSUS SITES ACROSS SAMPLES
peak_count = dict()
for i, sample in enumerate(samples):
    print(sample.name)
    # Get summits of peaks and window around them
    peaks = pybedtools.BedTool(sample.peaks)  # .slop(g=pybedtools.chromsizes('hg19'), b=250)
    # Merge overlaping peaks within a sample
    peaks = peaks.merge()
    if i == 0:
        sites = peaks
    else:
        # Concatenate all peaks
        sites = sites.cat(peaks)
    # Let's keep track of the number of new peaks found with each new sample
    peak_count[i] = len(sites)

# Merge overlaping peaks across samples
sites = sites.merge()

# Remove blacklist regions
blacklist = pybedtools.BedTool(os.path.join(dataDir, "wgEncodeDacMapabilityConsensusExcludable.bed"))
# Remove chrM peaks and save
sites = sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(dataDir, "all_sample_peaks.concatenated.bed"))

sites = pybedtools.BedTool(os.path.join(dataDir, "all_sample_peaks.concatenated.bed"))

# Plot cumulative number of peaks
plt.plot(peak_count.keys(), peak_count.values(), 'o')
plt.ylim(0, max(peak_count.values()) + 5000)
plt.title("Cumulative peaks per sample")
plt.xlabel("Number of samples")
plt.ylabel("Total number of peaks")
plt.savefig(os.path.join(plotsDir, "total_peak_count.per_patient.pdf"), bbox_inches="tight")
plt.close('all')

# Loop at summary statistics:
# interval lengths
sns.distplot([interval.length for interval in sites], bins=300, kde=False)
plt.xlabel("peak width (bp)")
plt.ylabel("frequency")
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.lengths.pdf"), bbox_inches="tight")
plt.close('all')

# calculate support (number of samples overlaping each merged peak)
for i, sample in enumerate(samples):
    print(sample.name)
    if i == 0:
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
support = pd.read_csv(os.path.join(dataDir, "all_sample_peaks.support.csv"))
# plot
sns.distplot(support["support"], bins=40)
plt.ylabel("frequency")
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.support.pdf"), bbox_inches="tight")
plt.close('all')


# MEASURE CHROMATIN OPENNESS PER SITE PER PATIENT
# Get read counts for each sample in the union of the sites

# Old method with bedtools
# coverage = sites.multi_bam_coverage(bams=[sample.filtered for sample in samples], p=True, q=30, D=True)
# coverage = sites.multi_bam_coverage(bams=[sample.filtered for sample in samples], q=30, D=True)
# make dataframe
# coverage = coverage.to_dataframe().reset_index()

# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]

# Count reads with pysam
# make strings with intervals
sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites]
# count, create dataframe
coverage = pd.DataFrame(
    map(
        lambda x:
            pd.Series(x),
            parmap.map(
                count_reads_in_intervals,
                [sample.filtered for sample in samples],
                sites_str,
                parallel=True
            )
    ),
    index=[sample.name for sample in samples]
).T
coverage.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.raw_coverage.tsv"), sep="\t", index=True)
coverage = pd.read_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.raw_coverage.tsv"), sep="\t", index_col=0)

# Add interval description to df
ints = map(
    lambda x: (
        x.split(":")[0],
        x.split(":")[1].split("-")[0],
        x.split(":")[1].split("-")[1]
    ),
    coverage.index
)
coverage["chrom"] = [x[0] for x in ints]
coverage["start"] = [int(x[1]) for x in ints]
coverage["end"] = [int(x[2]) for x in ints]

# Normalize by feature length (Reads per kilobase)
rpk = coverage.apply(normalize_by_interval_length, axis=1)
rpk.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.rpk.tsv"), sep="\t", index=True)
rpk = pd.read_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.rpk.tsv"), sep="\t", index_col=0)


# Normalize by library size - mapped reads (Reads per kilobase per million)
rpkm = rpk[[sample.name for sample in samples]].apply(normalize_by_library_size, args=(samples, ), axis=0)

# Save
rpkm = pd.concat([rpk[["chrom", "start", "end"]], rpkm], axis=1)
rpkm.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t", index=False)
rpkm = pd.read_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t")

# Log2 transform
rpkm[[sample.name for sample in samples]] = np.log2(1 + rpkm[[sample.name for sample in samples]])

# Plot
rpkmMelted = pd.melt(rpkm, id_vars=["chrom", "start", "end"], var_name="sample", value_name="rpkm")

# rpkm density
# all in one plot
for sample in samples:
    sns.distplot(rpkm[[sample.name]], hist=False, label=sample.name)
# plt.legend()
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.distplot.all.pdf"), bbox_inches="tight")
plt.close()

# separately in one grid
g = sns.FacetGrid(rpkmMelted, col="sample", aspect=2, col_wrap=4)
g.map(sns.distplot, "rpkm", hist=False)
plt.xlim(0, 15)
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.distplot.pdf"), bbox_inches="tight")
plt.close()

# boxplot rpkm per sample
# Plot the orbital period with horizontal boxes
sns.boxplot(x="rpkm", y="sample", data=rpkmMelted)
plt.xlim(0, 15)
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.boxplot.pdf"), bbox_inches="tight")
plt.close()

# pairwise rpkm scatter plot between samples
# fig = plt.figure(figsize=(8, 6), dpi=300, facecolor='w', edgecolor='k')
# g = sns.PairGrid(rpkm[[sample.name for sample in samples]])
# g.map(hexbin)
# g.fig.subplots_adjust(wspace=.02, hspace=.02)
# plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.pairwise_hexbin.pdf"), bbox_inches="tight")
# plt.close()

# Variation:
rpkm = pd.merge(rpkm, support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
# mean rpkm
rpkm['mean'] = rpkm[[sample.name for sample in samples]].apply(lambda x: np.mean(x), axis=1)
# dispersion (variance / mean)
rpkm['dispersion'] = rpkm[[sample.name for sample in samples]].apply(lambda x: np.var(x) / np.mean(x), axis=1)
# qv2 vs mean rpkm
rpkm['qv2'] = rpkm[[sample.name for sample in samples]].apply(lambda x: (np.std(x) / np.mean(x)) ** 2, axis=1)

sns.jointplot('mean', "dispersion", data=rpkm)
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.dispersion.pdf"), bbox_inches="tight")
plt.close('all')

sns.jointplot('mean', "qv2", data=rpkm)
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.qv2_vs_mean.pdf"), bbox_inches="tight")
plt.close('all')

sns.jointplot('support', "qv2", data=rpkm)
plt.savefig(os.path.join(plotsDir, "rpkm_per_sample.support_vs_qv2.pdf"), bbox_inches="tight")
plt.close('all')


# After inspecting known genotypes, consider excluding more regions/chromosomes across all samples
# e.g. chr11q, chr12


# Filter out regions which the maximum across all samples is below a treshold
filtered = rpkm[rpkm[[sample.name for sample in samples]].apply(max, axis=1) > 3]

sns.jointplot('mean', "dispersion", data=filtered)
sns.jointplot('mean', "qv2", data=filtered)


# GLOBAL CHARACTERIZATION
# Overview:

# get colors depending on IGVH mut
df = pd.DataFrame([sample.asSeries() for sample in samples])


df = pd.merge(df, clinical, left_on="sampleID", right_on="SampleID")
mut_s = {"1.0": "red", "2.0": "black", "nan": "grey"}
sex_s = {"F": "red", "M": "blue", "nan": "grey"}
muts = [clinical.loc[clinical['SampleID'] == sample.sampleID, "IGVH mut/unmut"] for sample in samples]
sex = [clinical.loc[clinical['SampleID'] == sample.sampleID, "gender"] for sample in samples]
mut_colors = [mut_s[str(x.get(x.index[0]))] if len(x.tolist()) > 0 else "grey" for x in muts]
sex_colors = [sex_s[str(x.get(x.index[0]))] if len(x.tolist()) > 0 else "grey" for x in sex]


# Correlation
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

c = rpkm[[sample.name for sample in samples]].corr()
c.index = map(name_to_repr, c.index)
c.columns = map(name_to_repr, c.columns)

sns.clustermap(
    c,
    method="complete",
    cmap=cmap,
    row_colors=mut_colors,
    col_colors=sex_colors,
    square=True, annot=False,
    figsize=(20, 16)
)
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.correlation_clustering.pdf"), bbox_inches="tight")

# PCA
# normalize
X = normalize(rpkm[[sample.name for sample in samples]])

# random PCA
pca = RandomizedPCA()
X = pca.fit(X).transform(X)

variance = [np.round(i * 100, 0) for i in pca.explained_variance_ratio_]

# get unique colors
colors = cm.Paired(np.linspace(0, 1, len(samples)))

# 2 components
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pca.components_[i, 0], pca.components_[i, 1],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
fig.axes[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
fig.axes[0].set_ylabel("PC2 - {0}% variance".format(variance[1]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-2comp.1vs2.pdf"), bbox_inches="tight")

# 2vs3 components
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pca.components_[i, 1], pca.components_[i, 2],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
fig.axes[0].set_xlabel("PC2 - {0}% variance".format(variance[1]))
fig.axes[0].set_ylabel("PC3 - {0}% variance".format(variance[2]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-2comp.2vs3.pdf"), bbox_inches="tight")

# 3vs4 components
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pca.components_[i, 2], pca.components_[i, 3],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
fig.axes[0].set_xlabel("PC3 - {0}% variance".format(variance[2]))
fig.axes[0].set_ylabel("PC4 - {0}% variance".format(variance[3]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-2comp.3vs4.pdf"), bbox_inches="tight")

# 4vs5 components
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pca.components_[i, 3], pca.components_[i, 4],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
fig.axes[0].set_xlabel("PC4 - {0}% variance".format(variance[3]))
fig.axes[0].set_ylabel("PC5 - {0}% variance".format(variance[4]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-2comp.4vs5.pdf"), bbox_inches="tight")

# 1vs3 components
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pca.components_[i, 0], pca.components_[i, 2],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
fig.axes[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
fig.axes[0].set_ylabel("PC3 - {0}% variance".format(variance[2]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-2comp.1vs3.pdf"), bbox_inches="tight")

# 3 components
fig = plt.figure()
# plot each point
ax = fig.add_subplot(111, projection='3d')
for i, sample in enumerate(samples):
    ax.scatter(
        pca.components_[i, 0], pca.components_[i, 1], pca.components_[i, 2],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=100
    )
ax.set_xlabel("PC1 - {0}% variance".format(variance[0]))
ax.set_ylabel("PC2 - {0}% variance".format(variance[1]))
ax.set_zlabel("PC3 - {0}% variance".format(variance[2]))
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.PCA-3comp.pdf"), bbox_inches="tight")

# MDS
# on 1000 most variable sites
n = 2000
Xvar = normalize(rpkm.ix[rpkm[[sample.name for sample in samples]].apply(np.var, axis=1).order(ascending=False).index].head(n)[[sample.name for sample in samples]])

# convert two components as we're plotting points in a two-dimensional plane
# "precomputed" because we provide a distance matrix
# we will also specify `random_state` so the plot is reproducible.
mds = MDS()  # n_components=2, dissimilarity="precomputed", random_state=1)

pos = mds.fit_transform(Xvar)

# plot
fig = plt.figure()
for i, sample in enumerate(samples):
    plt.scatter(
        pos[i, 0], pos[i, 1],
        label=name_to_repr(sample.name),
        color=colors[i],
        s=50
    )
plt.title("MDS on %i most variable genes" % n)
plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
plt.savefig(os.path.join(plotsDir, "all_sample_peaks.concatenated.MDS.pdf"), bbox_inches="tight")


# COMPARISON mCLL - uCLL
# test if sites come from same population
pvalues = rpkm.apply(
    lambda x: mannwhitneyu(
        x.loc[[sample.name for sample in samples if sample.mutated]],
        x.loc[[sample.name for sample in samples if not sample.mutated]])[1],
    axis=1
)

# correct for multiple testing
qvalues = multipletests(pvalues)[1]

# get differential sites
sigs = rpkm.ix[pvalues[pvalues < 0.0001].index][[sample.name for sample in samples]]
sigs.to_csv("clust.tsc", sep="\t", index=False)

# correlate samples on significantly different sites
sns.clustermap(
    sigs.corr(),
    method="complete",
    cmap=cmap,
    row_colors=mut_colors,
    col_colors=sex_colors,
    square=True, annot=False,
    figsize=(20, 16)
)
plt.savefig(os.path.join(plotsDir, "dors.correlation_clustering.pdf"), bbox_inches="tight")

clustermap = sns.clustermap(
    sigs,
    cmap=plt.get_cmap('YlGn'),
    square=True, annot=False,
)
plt.savefig(os.path.join(plotsDir, "dors.clustering.pdf"), bbox_inches="tight")

# get cluster assignments from linkage matrix
lm = clustermap.dendrogram_col.linkage

# plot dendrogram
# determine height to separate clusters
dendr = dendrogram(lm, labels=sigs.columns, color_threshold=65)
# plt.show()

# assign each sample to one cluster
clusters = get_cluster_classes(dendr)

# concatenate clusters on the unmutated side
unmutated_cluster = ['b', 'c', 'm', 'y']

# annotate samples with cluster
new_samples = list()
for cluster, sample_names in clusters.items():
    for sample_name in sample_names:
        s = [sample for sample in samples if sample.name == sample_name][0]
        s.cluster = cluster if cluster not in unmutated_cluster else 'b'
        new_samples.append(s)
samples = new_samples


# Repeat again independence test and
# get all differential sites (from 3 comparisons)
all_sig = pd.DataFrame()
for g1, g2 in itertools.combinations(['r', 'g', 'b']):
    pvalues = rpkm.apply(
        lambda x: mannwhitneyu(
            x.loc[[sample.name for sample in samples if sample.cluster == "r"]],
            x.loc[[sample.name for sample in samples if sample.cluster == "g"]]),
        axis=1
    )
    sig = rpkm.ix[pvalues[pvalues < 0.0001].index][[sample.name for sample in samples]]
    sig['comparison'] = "-".join(g1, g2)
    all_sig = pd.concat([all_sig, sig])

all_sig = all_sig.drop_duplicates()

all_sig.to_csv(os.path.join("data", "dors.3_populations.tsv"), sep="\t", index=False)


# correlate samples on significantly different sites
sns.clustermap(
    all_sig.corr(),
    method="complete",
    cmap=cmap,
    square=True, annot=False,
    figsize=(20, 16)
)
plt.savefig(os.path.join(plotsDir, "dors.3_populations.correlation_clustering.pdf"), bbox_inches="tight")

# heatmap all significat sites
clustermap = sns.clustermap(
    all_sig,
    cmap=plt.get_cmap('YlGn'),
    square=True, annot=False
)
plt.savefig(os.path.join(plotsDir, "dors.3_populations.clustering.pdf"), bbox_inches="tight")


# INTER-SAMPLE VARIABILITY ANALYSIS

# Subsample peaks or reads and see the minimum required to form the clusters previously

# See which regions explain most of variability for each cluster


# Get statistically different peaks between groups
#
#

#


# other options:

# Get most variable sites
# plot CV2 vs mean openness

# Enrichment of variable sites:
# LOLA
# get closest gene: GO, KEGG, OMIM, mSigDB
# de novo motif finding - enrichment


# See relationship betgween openness at enhancers and promoters
# across samples, get co-ocurrence frequency


# CLASSIFICATION
# Stratify patients on:
# - treated vs untreated
# - ...
# Train classifiers on groups
# Predict for all samples
# Assess: ROC, AUC
