#!/usr/bin/env python

"""
This script aims to build patient- or patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import os
from pipelines.models import Project, ATACseqSample
import re
import pandas as pd
import pybedtools
import multiprocessing
import parmap


def piq_to_network(results_dir, motif_numbers):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """

    def calculate_score(gene, chrom_footprints):
        """
        As defined in 10.1016/j.cels.2015.06.003
        """
        return sum(chrom_footprints.apply(lambda x: 2 * (x["purity"] - 0.5) * 10 ** (-abs(tsss.ix[gene]['start'] - x['start']) / 1000000), axis=1))

    # list results_dir
    files = os.listdir(results_dir)
    # get all cll peaks to filter data
    all_peaks = pybedtools.BedTool("data/cll_peaks.bed")
    # read in gene info
    tsss = pd.read_csv("data/hg19.refSeq.TSS.mRNA.bed", sep="\t", header=None)
    tsss.columns = ["chrom", "start", "end", "id"]

    # prepare TF vs Gene matrix
    scores = pd.DataFrame(index=tsss["id"], columns=motif_numbers)

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in motif_numbers:
        # get both forward and reverse complement PIQ output files
        result_files = list()
        for f in files:
            m = re.match(r'%i-.*\.csv$' % motif, f)
            if hasattr(m, "string"):
                result_files.append(m.string)

        # make bed file from it
        # concatenate files (forward and reverse complement are treated differently by PIQ)
        for i, result_file in enumerate(result_files):
            df = pd.read_csv(os.path.join(results_dir, result_file), index_col=0)
            df.rename(columns={"coord": "start"}, inplace=True)
            # fix coordinates
            if "RC-calls.csv" not in result_file:
                df["end"] = df["start"] + 1
            else:
                df["end"] = df["start"]
                df["start"] = df["start"] - 1
            # concatenate
            if i == 0:
                df2 = df
            else:
                df2 = pd.concat([df, df2])

        # Filter for purity
        footprints = df2[df2["purity"] > 0.7]

        # If empty give 0 to every gene for this TF
        if len(footprints) < 500:
            continue

        footprints[['chr', 'start', 'end', 'pwm', 'shape', 'score', 'purity']].to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # filter for motifs overlapping CLL peaks
        footprints = pybedtools.BedTool(os.path.join("tmp.bed")).intersect(all_peaks, wa=True).to_dataframe()
        footprints.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity"]

        # If empty give 0 to every gene for this TF
        if len(footprints) < 500:
            continue
        print("footprint number", len(footprints))

        # CONNECT
        # Now assign a score between this TF and every gene:
        # get distance to nearest gene TSS in the same chromosome as footprint

        for chrom in footprints["chrom"].unique():
            print(chrom)
            # calculate the distance between each footprint and every gene in the chromosome
            chrom_footprints = footprints[footprints["chrom"] == chrom]

            scores.loc[tsss[tsss["chrom"] == chrom]['id'], motif] = parmap.map(calculate_score, tsss[tsss["chrom"] == chrom].index, chrom_footprints)

    # everything else gets 0
    return scores.fillna(0)


# Get path configuration
data_dir = os.path.abspath(os.path.join('.', "data"))
scratch_dir = os.path.join("/scratch/users/arendeiro/piq")
results_dir = os.path.abspath(os.path.join('.', "results"))
plots_dir = os.path.join(results_dir, "plots")

# Get clinical info
clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))

# Get all CLL peaks
cll_peaks = os.path.join(data_dir, "cll_peaks.bed")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Select ATAC-seq samples
prj.samples = [s for s in prj.samples if type(s) == ATACseqSample]

to_exclude_sample_id = ['1-5-45960']

# FOOTPRINTING
motifs_file = "~/workspace/piq-single/pwms/jasparfix.txt"
n_motifs = 1316

df = pd.read_csv("data/tf_gene_matching.txt", sep="\t", header=None)
motif_numbers = df[0].tolist()

# prepare motifs for footprinting (done once)
# cmds = piq_prepare_motifs(motifs_file, n_motifs)
# for cmd in cmds:
#     os.system(cmd)

# parse output,
# connect each motif to a gene
for sample in prj.samples[1:5]:
    print(sample)
    scores = piq_to_network(os.path.join(sample.dirs.sampleRoot, "footprints"), motif_numbers)
    scores.to_csv(os.path.join(sample.dirs.sampleRoot, "footprints", "piq.TF-gene_scores.csv"))

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
