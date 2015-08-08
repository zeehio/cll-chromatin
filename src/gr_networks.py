#!/usr/bin/env python

"""
This script aims to build patient- or patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import yaml
import os
from pipelines.models import Project, ATACseqSample
from collections import OrderedDict
import numpy as np
import HTSeq
import pybedtools
import multiprocessing
import parmap
import cPickle as pickle
import pandas as pd


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
        center = interval.end - interval.start
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


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True, strand_specific=True):
    """
    Gets read coverage in genomic intervals.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.

    :param bam: HTSeq.BAM_Reader object, must be sorted and indexed with .bai file.
    :type bam: HTSeq.BAM_Reader
    :param intervals: dict with HTSeq.GenomicInterval objects as values.
    :type intervals: dict
    :type fragmentsize: int
    :type stranded: bool
    :type duplicates: bool
    :returns: OrderedDict with regionName:numpy.array(coverage)
    :rtype: collections.OrderedDict
    """
    chroms_avoid = ['chrM', 'chrX', 'chrY']
    cov = OrderedDict()
    n = len(intervals)
    i = 0

    # Loop through intervals, get coverage, append to dict
    for name, feature in intervals.iteritems():
        if i % 1000 == 0:
            print(n - i)
        # Initialize empty array for this feature
        if not strand_specific:
            profile = np.zeros(feature.length, dtype=np.float64)
        else:
            profile = np.zeros((2, feature.length), dtype=np.float64)

        # Check if feature is in bam index
        if feature.chrom in chroms_avoid:
            i += 1
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
                profile[start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window: end_in_window] += 1
                else:
                    profile[1][start_in_window: end_in_window] += 1

        # append feature profile to dict
        cov[name] = profile
        i += 1
    return cov


def call_footprints(cuts, annot, plot, mLen):
    """
    Call footprints.
    Requires dataframe with cuts and dataframe with annotation (2> cols).
    """
    import rpy2.robjects as robj  # for ggplot in R
    import rpy2.robjects.pandas2ri  # for R dataframe conversion

    # Plot with R
    footprint = robj.r("""
    library(CENTIPEDE)

    function(cuts, annot, plotFile, mLen) {
        centFit <- fitCentipede(
            Xlist = list(as.matrix(cuts)),
            Y = as.matrix(annot)
        )
        pdf(plotFile)
            plotProfile(centFit$LambdaParList[[1]],Mlen=mLen)
        dev.off()
        return(centFit$PostPr)
    }

    """)

    # convert the pandas dataframe to an R dataframe
    robj.pandas2ri.activate()
    cuts_R = robj.conversion.py2ri(cuts)
    annot_R = robj.conversion.py2ri(annot)

    # run the plot function on the dataframe
    return np.ndarray.flatten(robj.conversion.ri2py(footprint(cuts_R, annot_R, plot, mLen)))


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
data_dir = os.path.join(config["paths"]["parent"], config["projectname"], "data")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("../metadata/sequencing_sample_annotation.csv")

# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]


# FOOTPRINTING
# - Go to each TF motif site annotated: http://compbio.mit.edu/encode-motifs/
# - get window around (e.g. +/- 100bp), footprint with CENTIPEDE or PIQ
# - save posterior probabilities for every site

# Separate sites by TF, get TF motif length

# Build matrix of TF-site

# Polulate with footprint score


# Preprocess motifs data objects
# load motifs
motifs = pybedtools.BedTool(
    os.path.join(
        "/data/groups/lab_bock/shared/resources/genomes/hg19/motifs",
        "matches.merged.bed"
    )
)
# truncate to center base-pair
motifs = map(truncate_interval, motifs)
# add window around
motifs = motifs.slop(b=100)
# convert to HTSeq.GenomicInterval
motifs = map(bedtools_interval_to_genomic_interval, motifs)

# Generate some dummy motif scores
# this should be replaced with the actual scores later
motif_scores = pd.DataFrame(np.ones([10, 2]))

# Loop through samples, get coverage, save, footprint
covs = dict()
foots = dict()
for sample in samples:
    # Divide intervals by chromossome, do each in parallel
    # I need to find a way to split the intervals according to chromosome and create an iterator with
    # a further list of intervals of the same chromosome as elements

    # this needs to be reduced
    # perhaps a pd.DataFrame().append() od pd.concat() will do the trick
    covs[sample.name] = parmap.map(
        coverage, {motif.name: motif for motif in motifs},
        HTSeq.BAM_Reader(sample.mapped),  # be more stringent here, ask for sample.filtered
        fragmentsize=1
    )
    # serialize
    pickle.dump(covs, open(os.path.join(data_dir, "all_samples.motif_coverage.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Call footprints
    foots[sample.name] = call_footprints(pd.DataFrame(covs[sample.name]), motif_scores)
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
