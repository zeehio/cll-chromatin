#!/usr/bin/env python

"""
This script aims to build patient- or patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import yaml
import os
from pipelines.models import Project, ATACseqSample
import numpy as np
import HTSeq
import pybedtools
import multiprocessing
import parmap
import cPickle as pickle
import pandas as pd


def prepare_motifs():
    """
    Prepare motifs for footprinting.
    """
    motifs_dir = "/data/groups/lab_bock/shared/resources/genomes/hg19/motifs"
    cmds = list()

    # Get ENCODE motif matches
    cmds.append("wget http://compbio.mit.edu/encode-motifs/matches.txt.gz")  # get motifs
    cmds.append("gzip -d matches.txt.gz")  # extract
    cmds.append("tr ' ' \\t < matches.txt > matches.tsv")  # replace spaces with tabs
    cmds.append("""perl -ane 'print "$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[0]\n"' matches.tsv > matches.bed""")  # make bed file

    # Run file preparation
    for cmd in cmds:
        os.system(cmd)

    # make output dir
    if not os.path.exists(os.path.join(motifs_dir, "TFs")):
        os.mkdir(os.path.join(motifs_dir, "TFs"))

    # read master file, split per TF
    with open("matches.bed", "r") as handle:
        for line in handle:
            tf = line.split("\t")[4].split("_")[0]
            with open(os.path.join("TFs", tf + ".bed"), 'a') as handle2:
                handle2.write(line)


def truncate_interval(interval, center=True, n=1):
    """
    Truncate genomic interval to a number of basepairs from its center (if center is True) or its start.

    :param interval: a pybedtools.Interval object.
    :type interval: pybedtools.Interval
    :param center: If interval should be truncated to its center position
    :type center: bool
    :param n: Number of basepairs to truncate interval to.
    :type n: int
    :returns: Truncated pybedtools.Interval object.
    :rtype: pybedtools.Interval
    """
    if center:
        center = ((interval.end - interval.start) / 2) + interval.start
        interval.start = center
    interval.end = interval.start + n
    return interval


def bedtools_interval_to_genomic_interval(interval):
    """
    Given a pybedtools.Interval object, returns a HTSeq.GenomicInterval object.

    :param interval: a pybedtools.Interval object.
    :type interval: pybedtools.Interval
    :returns: HTSeq.GenomicInterval object.
    :rtype: HTSeq.GenomicInterval
    """
    if interval.strand == "+" or interval.strand == 0 or interval.strand == str(0):
        return HTSeq.GenomicInterval(interval.chrom, interval.start, interval.end, "+")
    elif interval.strand == "-" or interval.strand == 0 or interval.strand == str(1):
        return HTSeq.GenomicInterval(interval.chrom, interval.start, interval.end, "-")
    else:
        return HTSeq.GenomicInterval(interval.chrom, interval.start, interval.end)


def coverage(bed_file, bam_file, sites, fragmentsize=1, orientation=True, duplicates=True, strand_specific=True):
    """
    Gets read coverage in genomic intervals.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.

    :param bed_file: Bed file.
    :type bed_file: str
    :param bam: HTSeq.BAM_Reader object, must be sorted and indexed with .bai file.
    :type bam: HTSeq.BAM_Reader
    :type fragmentsize: int
    :type stranded: bool
    :type duplicates: bool
    :returns: OrderedDict with regionName:numpy.array(coverage)
    :rtype: collections.OrderedDict
    """
    motifs = pybedtools.BedTool(bed_file)
    # add window around
    intervals = motifs.slop(b=100, genome="hg19")

    intervals = intervals.intersect(b=sites, u=True)

    # convert to HTSeq.GenomicInterval
    intervals = map(bedtools_interval_to_genomic_interval, intervals)

    # Bam file
    bam = HTSeq.BAM_Reader(bam_file)

    chroms_avoid = ['chrM', 'chrX', 'chrY']
    n = len(intervals)
    m = intervals[0].length

    # create empty array
    if not strand_specific:
        coverage = np.zeros((n, m), dtype=np.float64)
    else:
        coverage = np.zeros((n, m * 2), dtype=np.float64)

    # Loop through intervals, get coverage, append to dict
    for i, feature in enumerate(intervals):
        if i % 1000 == 0:
            print(n - i)

        # Check if feature is in bam index
        if feature.chrom in chroms_avoid:
            continue

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            # check it's aligned
            if not aln.aligned:
                continue

            aln.iv.length = fragmentsize  # adjust to size

            # get position in relative to window
            if orientation:
                if feature.strand == "+" or feature.strand == ".":
                    start_in_window = aln.iv.start - feature.start - 1
                    end_in_window = aln.iv.end - feature.start - 1
                else:
                    start_in_window = feature.length - abs(feature.start - aln.iv.end) - 1
                    end_in_window = feature.length - abs(feature.start - aln.iv.start) - 1
            else:
                start_in_window = aln.iv.start - feature.start - 1
                end_in_window = aln.iv.end - feature.start - 1

            # check fragment is within window; this is because of fragmentsize adjustment
            if start_in_window <= 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                coverage[i, start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    coverage[i, start_in_window: end_in_window] += 1
                else:
                    coverage[i, m + start_in_window: m + end_in_window] += 1
    return coverage


def call_footprints(cuts, annot, plot, motif_length):
    """
    Call footprints.
    Requires dataframe with cuts and dataframe with annotation (2> cols).
    """
    import rpy2.robjects as robj  # for ggplot in R
    import rpy2.robjects.pandas2ri  # for R dataframe conversion
    import rpy2.robjects.numpy2ri  # for R numpy objects conversion

    robj.pandas2ri.activate()

    # Plot with R
    footprint = robj.r("""
    library("CENTIPEDE")

    function(cuts, annot, plotFile, mLen) {
        centFit <- fitCentipede(
            Xlist = list(as.matrix(cuts)),
            Y = as.matrix(annot),
            sweeps = 1000
        )
        pdf(plotFile)
            plotProfile(centFit$LambdaParList[[1]],Mlen=mLen)
        dev.off()
        return(centFit$PostPr)
    }

    """)

    # run the plot function on the dataframe
    return footprint(cuts, annot, plot, motif_length)


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)

# Project dirs
data_dir = os.path.join(config["paths"]["parent"], config["projectname"], "data")
results_dir = os.path.join(config["paths"]["parent"], config["projectname"], "results")
plots_dir = os.path.join(results_dir, "plots")

# Directory for footprint plots
if not os.path.exists(os.path.join(plots_dir, "footprints")):
    os.mkdir(os.path.exists(os.path.join(plots_dir, "footprints")))

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]


# FOOTPRINTING
# Get all cll sites
sites = pybedtools.BedTool(os.path.join(data_dir, "all_sample_peaks.concatenated.bed"))
# Use motifs from here: http://compbio.mit.edu/encode-motifs/
# Split motifs by TF
motifs_dir = "/data/groups/lab_bock/shared/resources/genomes/hg19/motifs"
prepare_motifs()

# Get all TF motifs
motifs = os.listdir(os.path.join("/data/groups/lab_bock/shared/resources/genomes/hg19/motifs/TFs"))
motif_names = [motif.split(".")[0] for motif in motifs]
motif_files = [os.path.abspath(os.path.join("/data/groups/lab_bock/shared/resources/genomes/hg19/motifs/TFs", motif)) for motif in motifs]
motif_sizes = list()
for motif_file in motif_files:
    fields = open(motif_file, 'r').readline().split("\t")
    motif_sizes.append(int(fields[2]) - int(fields[1]))


# Get window around (e.g. +/- 100bp) each motif
# Get coverage there
# Footprint with centipede
# Save posterior probabilities for every site

# Loop through samples, get coverage, save, footprint
covs = dict()
foots = dict()
for sample in samples:
    # calculate coverage in parallel for several TFs
    covs[sample.name] = parmap.map(
        coverage,
        motif_files,
        sample.filtered,  # be more stringent here, ask for sample.filtered
        sites
    )
    # serialize
    pickle.dump(covs, open(os.path.join(data_dir, "all_samples.motif_coverage.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Call footprints
    for i, coverage in enumerate(covs[sample.name]):
        try:
            foots[sample.name] = call_footprints(
                coverage[:, 1:-1],
                np.ones([len(coverage), 1]),
                os.path.join(plots_dir, "footprints", sample.name + "." + motif_names[i] + ".pdf"),
                motif_sizes[i]
            )
        except:
            continue

    # serialize
    pickle.dump(foots, open(os.path.join(data_dir, "all_samples.motif_coverage.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)


# more:
# get factor occupancy score (FOS) for each factor in each site
# cluster again patients based on TF occupancy


# NETWORKS
# Network construction:

# Calculate a "regulation score" between EVERY TF and EVERY gene:
# - for each TF, for each chromossome, sum the weighted posterior probabilities of factor A regulating gene I.
# weigh these by dividing by the (log) of the distance of each binding site to the TSS of gene I.
# This will give a matrix of TS-gene "regulation scores".
# Investigate the distribution of scores.

# Compare regulation across patients:
# - for each TF, correlate the scores with the scores of all other transcription factors in the other patient.

# Additionally
# - build network of TF -> gene based on:
#   - binding at the promoter or at known enhancers for the gene
#   OR
#   - consider regulation only if score is above threshold

# Network types:
# patient-specific:
# - build network specific to each patient
# - compare networks
# - cluster patients based on network connections

# for groups of patients:
# - come up with a way of combining signal from several patients from one group
# - build networks specific to groups

# Compare Networks


#


# CIS-REGULATORY MODULE USAGE BETWEEN CLL TYPES
# Get CLL groups
# Select genome region of relevance
# Get reads for regions at e.g. 1kb resolution
# Correlate genome positions in each group

# Get differential between groups
