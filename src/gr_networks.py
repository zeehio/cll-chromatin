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


def save_pandas(data, fname):
    '''Save DataFrame or Series
    Parameters
    ----------
    fname : str
        filename to use
    data: Pandas DataFrame or Series
    '''
    np.save(open(fname, 'w'), data)
    if len(data.shape) == 2:
        meta = data.index, data.columns
    elif len(data.shape) == 1:
        meta = (data.index,)
    else:
        raise ValueError('save_pandas: Cannot save this type')
    s = pickle.dumps(meta)
    s = s.encode('string_escape')
    with open(fname, 'a') as f:
        f.seek(0, 2)
        f.write(s)


def load_pandas(fname, mmap_mode='r'):
    '''Load DataFrame or Series
    Parameters
    ----------
    fname : str
        filename
    mmap_mode : str, optional
        Same as numpy.load option
    '''
    values = np.load(fname, mmap_mode=mmap_mode)
    with open(fname) as f:
        np.lib.format.read_magic(f)
        np.lib.format.read_array_header_1_0(f)
        f.seek(values.dtype.alignment * values.size, 1)
        meta = pickle.loads(f.readline().decode('string_escape'))
    if len(meta) == 2:
        return pd.DataFrame(values, index=meta[0], columns=meta[1])
    elif len(meta) == 1:
        return pd.Series(values, index=meta[0])


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


def prepare_intervals(sites):
    # Get all TF motifs
    motifs_files = os.listdir(os.path.join("/data/groups/lab_bock/shared/resources/genomes/hg19/motifs/TFs"))
    motif_files = [os.path.abspath(os.path.join("/data/groups/lab_bock/shared/resources/genomes/hg19/motifs/TFs", motif)) for motif in motifs_files]

    # intersect motifs with all cll sites
    for motif_file in motif_files:
        # read in bedfile
        motifs = pybedtools.BedTool(motif_file)
        # get motif name
        motif_name = os.path.basename(motif_file.split(".")[0])

        # keep only motifs overlaping the consensus set of sites
        motifs = motifs.intersect(b=sites, u=True)

        # add window around
        intervals = motifs.slop(b=100, genome="hg19")

        # save new files
        intervals.saveas(os.path.join(data_dir, 'motifs', motif_name + ".bed"))


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


def footprint(bed_file, bam_file, sites, sample_name, fragmentsize=1, orientation=True, duplicates=True, strand_specific=True):
    """
    Gets read coverage in genomic intervals. Passes coverage to call_footprints and returns posterior probabilities.

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
    # read in bedfile
    motifs = pybedtools.BedTool(bed_file)
    # get motif name
    motif_name = os.path.basename(bed_file.split(".")[0])
    # get motif length (length of first interval)
    motif_length = motifs[0].length

    # convert intervals to HTSeq.GenomicInterval
    intervals = map(bedtools_interval_to_genomic_interval, motifs)

    # Handle bam file
    bam = HTSeq.BAM_Reader(bam_file)

    # exclude bad chroms
    chroms_exclude = ['chrM', 'chrX', 'chrY']

    # get dimensions of matrix to store profiles of Tn5 transposition
    n = len(intervals)
    m = intervals[0].length

    # create empty matrix
    if not strand_specific:
        coverage = np.zeros((n, m), dtype=np.float64)
    else:
        # if "strand_specific", get signal for both strands independently, but concatenated
        coverage = np.zeros((n, m * 2), dtype=np.float64)

    # Loop through intervals, get coverage, increment matrix count
    for i, feature in enumerate(intervals):
        # counter just to track
        if i % 1000 == 0:
            print(n - i)

        # Check if feature is not in bad chromosomes
        if feature.chrom in chroms_exclude:
            continue

        # Fetch alignments in interval
        for aln in bam[feature]:
            # check it's aligned
            if not aln.aligned:
                continue

            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue

            aln.iv.length = fragmentsize  # adjust reads to specified size

            # get position relative to window if required (motif-oriented)
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
            if start_in_window < 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                coverage[i, start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    coverage[i, start_in_window: end_in_window] += 1
                else:
                    coverage[i, m + start_in_window: m + end_in_window] += 1
    # Call footprints, get posterior probabilities
    try:
        probs = call_footprints(coverage, np.ones([len(coverage), 1]), motif_length, os.path.join(plots_dir, "footprints", sample_name + "." + motif_name + ".pdf"))
        if len(probs) != len(coverage):
            probs = np.zeros(len(coverage))
    except:
        # if error, return zeros
        probs = np.zeros(len(coverage))
    return probs


def call_footprints(cuts, annot, motif_length, plot):
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
    return footprint(cuts, annot, plot, motif_length).flatten()


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

# make motif bed files for each TF from bulk file
prepare_motifs()

# make intervals out of motifs
prepare_intervals(sites)

# get filtered motif files
motifs_files = [os.path.join(data_dir, "motifs", x) for x in os.listdir(os.path.join(data_dir, 'motifs'))]

# Get coverage around motifs there
# Footprint with centipede
# Get posterior probabilities for every site in an array
# Reduce by concatenation all arrays for various TFs into one long array
# Create matrix sites vs samples with probabilities as values

# Loop through samples, get coverage, footprint, reduce
# Make dataframe
all_probs = pd.DataFrame(columns=[sample.name for sample in samples])

try:
    for sample in samples[1:]:
        # calculate posterior probabilities in parallel for several TFs
        probs = reduce(
            lambda x, y: np.concatenate([x, y]),
            parmap.map(
                footprint,
                motifs_files,
                sample.filteredshifted,  # be stringent here, ask for sample.filtered
                sites,
                sample.name
            )
        )
        all_probs[sample.name] = probs
        # serialize
        save_pandas(all_probs, os.path.join(data_dir, "all_samples.footprint_probabilities.pdy"))
except KeyboardInterrupt:
    pass


# "all_probs" can further be seen as a multiindex dataframe with indice levels of:
#   TF
#   position (chrom, start, end)

# connect each motif to a gene:
# get TSS annotation
# get nearest TSS from motif
# for each motif with posterior probability > 0.99, create relation TF -> gene


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


# Examples:
# Open chromatin defined by DNaseI and FAIREidentifies regulatory elements that shape cell-type identity
# Lingyun Song et al.
# Figure 6A

# Patterns of regulatory activity across diverse human cell types predict tissue identity, transcription factor binding, and long-range interactions
# Nathan C. Sheffield et al.
