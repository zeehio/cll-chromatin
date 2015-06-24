#!/usr/bin/env python

import yaml
import os
from pipelines import Project, ATACseqSample
import pybedtools
from collections import Counter


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
# Get union of all sample's peaks
# with merged ovelaping intervals
for i, sample in enumerate(samples):
    if i == 0:
        sites = pybedtools.BedTool(sample.filteredPeaks)
    else:
        sites = sites.cat(pybedtools.BedTool(sample.filteredPeaks))

sites.saveas(os.path.join(dataDir, "peaks_concatenated.bed"))

# Loop at summary statistics:
# interval lengths, calculate support (number of samples overlaping each final one) etc...
lengths = Counter([interval.length for interval in sites])
support = sites.multi_intersect(i=[sample.filteredPeaks for sample in samples])
# plot these frequencies now

# Remove some regions chtMT, chr.*random

# After inspecting phenotypes, consider excluding more regions/chromosomes across all samples


# GET MEASURE OF CHROMATIN OPENNESS PER SITE PER PATIENT
# Get read counts for each sample in the union of the sites
# maybe consider enlarging a bit the sites (slop)

# e.g.:
# bedtools annotate -names [sample.name for sample in samples] -i consensus_peaks.bed -files [sample.filtered for sample in samples]
coverage = sites.multi_bam_coverage([sample.filtered for sample in samples], p=True, q=30)
coverage.saveas(os.path.join(dataDir, "peaks_concatenated.coverage.bed"))

# normalize counts by library size
# e.g. try using DESeq2 sizeFactors() to get resizing factors


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
