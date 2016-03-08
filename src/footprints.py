
import os
import re
from collections import OrderedDict
import numpy as np
import HTSeq
import pybedtools
import multiprocessing
import parmap

import matplotlib.pyplot as plt
import seaborn as sns


def bedTools_interval_to_genomic_interval(bedtool):
    """
    Given a pybedtools.BedTool object, returns dictionary of HTSeq.GenomicInterval objects.
    """
    intervals = OrderedDict()
    for iv in bedtool:
        name = "{}:{}-{}".format(iv.chrom, iv.start, iv.end)
        if iv.strand == "+" or iv.strand == 0 or iv.strand == str(0):
            intervals[name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
        elif iv.strand == "-" or iv.strand == 0 or iv.strand == str(1):
            intervals[name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
        else:
            intervals[name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)
    return intervals


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    stranded - boolean.
    duplicates - boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']
    cov = OrderedDict()
    n = len(intervals)
    i = 0
    try:
        for name, feature in intervals.iteritems():
            if i % 1000 == 0:
                print(n - i)
            # Initialize empty array for this feature
            if not strand_specific:
                profile = np.zeros(feature.length, dtype=np.float64)
            else:
                profile = np.zeros((2, feature.length), dtype=np.float64)

            # Check if feature is in bam index
            if feature.chrom not in chroms or feature.chrom == "chrM":
                i += 1
                continue

            # Replace chromosome reference 1 -> chr1 if not chr
            if switchChromsNames:
                feature.chrom = re.sub("chr", "", feature.chrom)

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
    except KeyboardInterrupt:
        return cov

    return cov


def coverage_single(feature, bam, fragmentsize=1, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    stranded - boolean.
    duplicates - boolean.
    """
    # Initialize empty array for this feature
    if not strand_specific:
        profile = np.zeros(feature.length, dtype=np.float64)
    else:
        profile = np.zeros((2, feature.length), dtype=np.float64)

    # Replace chromosome reference 1 -> chr1 if not chr
    if switchChromsNames:
        feature.chrom = re.sub("chr", "", feature.chrom)

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
    return profile


def get_reads_in_intervals(bam, intervals, strand_specific=False):
    """
    Counts reads in a iterable holding strings
    representing genomic intervals of the type chrom:start-end.
    """
    # count, create dataframe
    coverage = parmap.map(
        coverage_single,
        intervals.values(),
        bam,
        strand_specific=strand_specific,
        parallel=True)

    if not strand_specific:
        coverage = np.vstack(coverage)
    else:
        coverage = (np.vstack([x[0] for x in coverage]), np.vstack([x[1] for x in coverage]))
    return coverage
    # plt.plot(coverage.mean(0))

    # Test
    # coverage_single(intervals.values()[0], bam, strand_specific=False)
    # coverage_single(intervals.values()[:50], bam, strand_specific=True)


# INIT
# Get path configuration
data_dir = os.path.join('.', "data_submission")
scratch_dir = os.path.join("/scratch/users/arendeiro/piq")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# motifs to footprint
motifs = [
    "AP1", "EBF1", "RUNX3", "NFKB", "IRF4", "POU2F2",
    "PAX5",
    "CTCF"]

# groups of samples to footprint
groups = {
    "CLL": "merged-samples_all_all/merged.bam",
    "CD19": "external/CD19_DNase.bam"}

# all CLL peaks
peaks = pybedtools.BedTool("data/cll_peaks.bed")

for motif in motifs:
    # Prepare motif files
    bed = pybedtools.BedTool("/home/arendeiro/resources/genomes/hg19/motifs/TFs/{}.true.bed".format(motif))
    # intersect with peaks
    bed = bed.intersect(peaks, wa=True)
    # get leftmost basepair
    bed = bed.flank(s=True, l=1, r=0, genome="hg19")
    # get window around
    bed = bed.slop(b=250, genome="hg19")
    # convert to HTSeq
    intervals = bedTools_interval_to_genomic_interval(bed)

    print(motif, len(intervals))

    for group, file in groups.items():
        print(motif, group)
        # get coverage, strand specfic
        bam = HTSeq.BAM_Reader(os.path.join("data", file))

        # get coverage
        cov = get_reads_in_intervals(bam, intervals)

        # save to disk
        np.savetxt(os.path.join(data_dir, "footprint.cuts_matrix.{}.{}.csv".format(motif, group)), cov, delimiter=",", fmt='%1.0f')
        # cov = np.loadtxt(os.path.join(data_dir, "footprint.cuts_matrix.{}.{}.csv".format(motif, group)), delimiter=",")

        # plot footprint
        fig, axis = plt.subplots(1)
        axis.plot(cov.mean(0)[1:-1])
        fig.savefig(os.path.join(plots_dir, "footprint.cuts_matrix.{}.{}.svg".format(motif, group)), bbox_inches="tight")

        # cluster regions
        # sns.clustermap(cov, col_cluster=False, row_cluster=True)


# Cuts in specific regions
mCLL = HTSeq.BAM_Reader("data/merged-samples_mutated_True/merged.sorted.bam")
uCLL = HTSeq.BAM_Reader("data/merged-samples_mutated_False/merged.sorted.bam")
allCLL = HTSeq.BAM_Reader("data/merged-samples_all_all/merged.bam")

import HTSeq
intervals = {
    "PAX9": HTSeq.GenomicInterval("chr14", 37069382, 37148611),
    "PAX9_intron": HTSeq.GenomicInterval("chr14", 37130494, 37132045)
}
fragmentsize = 1

covT = coverage(mCLL, intervals, fragmentsize, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False)
covF = coverage(uCLL, intervals, fragmentsize, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False)
covall = coverage(allCLL, intervals, fragmentsize, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False)

plt.plot(covF["PAX9_intron"])
plt.plot(covT["PAX9_intron"])
plt.plot(covall["PAX9_intron"])

# "chr14:37069382-37148611"
