#!/usr/bin/env python

import yaml
import os
from pipelines import Project, ATACseqSample
import pybedtools
from collections import Counter
import pysam


def normalizeByIntervalLength(series):
    """
    Divides a pandas.Series with numerical values by the length of the interval
    (defined in "chr", "start", "end") in a bam file (RPK).
    """
    length = float(series["end"] - series["start"])

    rpK = (series.drop(["chr", "start", "end", "name", "score", "strand"]) / length) * 1e3

    return series[["chr", "start", "end", "name", "score", "strand"]].append(rpK)


def normalizeByLibrarySize(series, samples):
    """
    Divides a pandas.Series with numerical values by the total number of mapped
    reads in a bam file.
    """
    # get sample which this series belongs to
    sample = [sample for sample in samples if sample.name == series.name]

    # get number of aligned reads for that sample
    size = float(pysam.AlignmentFile(sample.filtered).mapped)

    return (series / size) * 1e6


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
dataDir = os.path.join(config["paths"]["parent"], config["projectname"], "data")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("../metadata/sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]


# GET CONSENSUS SITES ACROSS SAMPLES
for i, sample in enumerate(samples):
    # Get summits of peaks and window around them
    peaks = pybedtools.BedTool(sample.filteredPeaks).slop(g=pybedtools.chromsizes('hg19'), b=250)
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
# interval lengths, calculate support (number of samples overlaping each final one) etc...
lengths = Counter([interval.length for interval in sites])
support = sites.multi_intersect(i=[sample.filteredPeaks for sample in samples])
# plot these frequencies now

# Remove some regions chtMT, chr.*random

# After inspecting phenotypes, consider excluding more regions/chromosomes across all samples
# e.g. chr11q, chr12

# GET MEASURE OF CHROMATIN OPENNESS PER SITE PER PATIENT
# Get read counts for each sample in the union of the sites
coverage = sites.multi_bam_coverage(bams=[sample.filtered for sample in samples], p=True, q=30, D=True)
# make dataframe
sampleNames = [sample.name for sample in samples]
coverage = coverage.to_dataframe(names=["chr", "start", "end", "name", "score", "strand"] + sampleNames)
coverage.to_csv(os.path.join(dataDir, "all_sample_peaks.concatenated.raw_coverage.bed"), sep="\t", index=False)

# Divide by feature length
coverage.apply(normalizeByIntervalLength, axis=1)
# Divide by library size (mapped reads)
coverage.apply(normalizeByLibrarySize, args=(samples, ), axis=0)

coverage.saveas(os.path.join(dataDir, "all_sample_peaks.concatenated.fpkm.bed"))


# GLOBAL CHARACTERIZATION
# Overview:
# Heatmap sites vs patients
# hierarchical clustering

# heatmap with only sites near known CLL factors


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
