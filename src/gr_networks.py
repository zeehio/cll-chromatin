#!/usr/bin/env python

"""
This script aims to build patient- or patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import os
from pipelines.models import Project, ATACseqSample
from pipelines import toolkit as tk
import textwrap
import re
import pandas as pd
import pybedtools
import networkx as nx


def name_to_sample_id(name):
    return name.split("_")[3:4][0]


def annotate_igvh_mutations(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            if clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 1:
                sample.mutated = True
            elif clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 2:
                sample.mutated = False
            else:
                sample.mutated = None
        else:
            sample.mutated = None
        new_samples.append(sample)
    return new_samples


def piq_prepare_motifs(motif_file="~/workspace/piq-single/pwms/jasparfix.txt", n_motifs=1316):
    """
    """
    cmds = list()
    for motif in range(1, n_motifs + 1):
        cmd = """
    Rscript ~/workspace/piq-single/pwmmatch.exact.r"""
        cmd += " ~/workspace/piq-single/common.r"
        cmd += " {0}".format(motif_file)
        cmd += " " + str(motif)
        cmd += """ /scratch/users/arendeiro/piq/motif.matches/
    """
        cmds.append(cmd)
    return cmds


def piq_prepare_bams(bams, output_cache):
    """
    :type bams: list
    :type output_cache: str
    """
    # add before:
    # salloc -c 12 -p develop --mem-per-cpu 8000
    cmd = """
    Rscript ~/workspace/piq-single/bam2rdata.r ~/workspace/piq-single/common.r {0} """.format(output_cache)
    cmd += " ".join(bams)
    cmd += """
    """
    return cmd


def piq_footprint(bam_cache, motif_numbers, tmp_dir, results_dir):
    """
    Footprint using PIQ.
    """
    cmds = list()
    for motif in motif_numbers:
        cmd = """
    Rscript ~/workspace/piq-single/pertf.r"""
        cmd += " ~/workspace/piq-single/common.r"
        cmd += " /scratch/users/arendeiro/piq/motif.matches/"
        cmd += " " + os.path.join(tmp_dir, str(motif))
        cmd += " " + results_dir
        cmd += " " + bam_cache
        cmd += " " + str(motif)
        cmd += """
    """
        cmds.append(cmd)

    return cmds


def piq_footprint_single(bam_cache, motif_number, tmp_dir, results_dir):
    """
    Footprint using PIQ.
    """
    cmd = """
    Rscript ~/workspace/piq-single/pertf.r"""
    cmd += " ~/workspace/piq-single/common.r"
    cmd += " /scratch/users/arendeiro/piq/motif.matches/"
    cmd += " " + tmp_dir
    cmd += " " + results_dir
    cmd += " " + bam_cache
    cmd += " " + str(motif_number)
    cmd += """
"""
    return cmd


def piq_to_network(results_dir, motif_numbers):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """
    # list results_dir
    files = os.listdir(results_dir)
    # get all cll peaks to filter data
    all_peaks = pybedtools.BedTool("data/cll_peaks.bed")
    # read in gene info
    refseq_mrna_tss = pybedtools.BedTool("data/hg19.refSeq.TSS.mRNA.deduplicated.bed")

    # dict to store TF->gene interactions
    interactions = pd.DataFrame()

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in motif_numbers:
        # get both forward and reverse complement PIQ output files
        result_files = list()
        for f in files:
            m = re.match(r'%i-.*-calls\.csv$' % motif, f)
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
        footprints.to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # Get closest gene
        closest_tss = pybedtools.BedTool(os.path.join("tmp.bed")).closest(refseq_mrna_tss, d=True).to_dataframe()
        closest_tss.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity", "chrom_gene", "start_gene", "end_gene", "gene", "distance"]

        # Get weighted values
        # weigh with footprint purity and distance to tss
        scores = closest_tss
        scores['interaction_score'] = scores.apply(lambda x: 2 * (x['purity'] - 0.5) * 10 ** -(x['distance'] / 1000000.), axis=1)
        # sum scores for each gene
        scores = scores.groupby(['gene'])['interaction_score'].apply(sum).reset_index()

        scores["TF"] = motif

        interactions = pd.concat([interactions, scores])

    return interactions


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
samples_to_exclude = [
    'CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19',
    'CLL_ATAC-seq_5186_1-5-57350_ATAC17-4_hg19',
    'CLL_ATAC-seq_4784_1-5-52817_ATAC17-6_hg19',
    'CLL_ATAC-seq_981_1-5-42480_ATAC16-6_hg19',
    'CLL_ATAC-seq_5277_1-5-57269_ATAC17-8_hg19',
    'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19',
    'CLL_ATAC-seq_5147_1-5-48105_ATAC17-2_hg19',
    'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19']  # 'CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19']
samples = [s for s in prj.samples if type(s) == ATACseqSample and s.name not in samples_to_exclude]

samples = annotate_igvh_mutations(samples, clinical)

# FOOTPRINTING
motifs_file = "~/workspace/piq-single/pwms/jasparfix.txt"
n_motifs = 1316

# read list of tfs to do
df = pd.read_csv("data/tf_gene_matching.txt", sep="\t", header=None)
df[1] = [x.upper() for x in df[1]]
tfs = df[1]
number2tf = dict(zip(df[0], df[1]))
motif_numbers = df[0]

# get refseq -> gene symbol mapping
os.system("""mysql --user=genome -N --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e "select name,name2 from refGene" > data/Refseq2Gene.txt""")
refseq2gene = pd.read_csv("data/Refseq2Gene.txt", sep="\t", header=None)
refseq2gene = dict(zip(refseq2gene[0], refseq2gene[1]))

# stupid PIQ hard-coded links
os.chdir("/home/arendeiro/workspace/piq-single/")

# prepare motifs for footprinting (done once)
cmds = piq_prepare_motifs(motifs_file, n_motifs)
for cmd in cmds:
    cmd2 = tk.slurmHeader("PIQ_preparemotifs", os.path.join("/home/arendeiro/", "piq_preparemotifs.slurm.log"), cpusPerTask=1, queue="shortq")

    # stupid PIQ hard-coded links
    cmd2 += cmd

    # write job to file
    with open("/home/arendeiro/tmp.sh", 'w') as handle:
        handle.writelines(textwrap.dedent(cmd2))

    tk.slurmSubmitJob("/home/arendeiro/tmp.sh")


# for each sample create R cache with bam file
jobs = list()
for sample in samples:
    if sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(sample.dirs.sampleRoot, "footprints")
    r_data = os.path.join(foots_dir, sample.name + ".filteredshifted.RData")
    if os.path.isfile(r_data):
        continue
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)

    tmp_dir = os.path.join(scratch_dir, sample.name)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    job_file = os.path.join(foots_dir, "slurm_job.sh")

    # prepare slurm job header
    cmd = tk.slurmHeader(sample.name + "_PIQ_prepareBam", os.path.join(foots_dir, "slurm.log"), cpusPerTask=2, queue="shortq")

    # stupid PIQ hard-coded links
    cmd += """
    cd /home/arendeiro/workspace/piq-single/
    """

    # prepare bams
    cmd += piq_prepare_bams([sample.filteredshifted], r_data)

    # slurm footer
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    # append file to jobs
    jobs.append(job_file)

# submit jobs (to create bam file caches)
for job in jobs:
    tk.slurmSubmitJob(job)


# for each sample launch several jobs (>500) to footprint
jobs = list()
for sample in samples:
    if sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(sample.dirs.sampleRoot, "footprints")
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)
    r_data = os.path.join(foots_dir, sample.name + ".filteredshifted.RData")

    # if footprint files exist, skip sample:
    if os.path.exists(os.path.join(foots_dir, "936-PB00421Mafk1-diag.pdf")):  # this is an example
        continue
    else:
        print(sample.name)

    for motif in motif_numbers:
        if not os.path.exists("/scratch/users/arendeiro/piq/motif.matches/%i.pwmout.RData" % motif):
            continue

        t_dir = os.path.join(scratch_dir, sample.name)
        if not os.path.exists(t_dir):
            os.mkdir(t_dir)
        tmp_dir = os.path.join(scratch_dir, sample.name, str(motif))
        # if not os.path.exists(tmp_dir):
        #     os.mkdir(tmp_dir)
        job_file = os.path.join(foots_dir, "slurm_job_motif%i.sh" % motif)
        slurm_log = os.path.join(foots_dir, "slurm_motif%i.log" % motif)

        # prepare slurm job header
        cmd = tk.slurmHeader(sample.name + "_PIQ_footprinting_motif%i" % motif, slurm_log, cpusPerTask=2, queue="shortq", memPerCpu=8000)

        # stupid PIQ hard-coded links
        cmd += """
    cd /home/arendeiro/workspace/piq-single/
        """

        # footprint
        cmd += piq_footprint_single(r_data, motif, tmp_dir, results_dir=foots_dir)

        # slurm footer
        cmd += tk.slurmFooter()

        # write job to file
        with open(job_file, 'w') as handle:
            handle.writelines(textwrap.dedent(cmd))

        # append file to jobs
        jobs.append(job_file)


# submit jobs (to footprint)
for job in jobs:
    tk.slurmSubmitJob(job)


# parse output,
# connect each motif to a gene
for sample in samples:
    print sample
    interactions = piq_to_network(os.path.join(sample.dirs.sampleRoot, "footprints"), motif_numbers)

    # Drop TFBS with no gene in the same chromosome (random chroms) - very rare cases
    interactions = interactions[interactions['gene'] != '.']

    # Get original TF name and gene symbol
    interactions['TF'] = [number2tf[tf] for tf in interactions['TF']]
    interactions['gene'] = [refseq2gene[gene] for gene in interactions['gene']]

    interactions['interaction_type'] = "pd"
    interactions.to_csv(os.path.join(data_dir, "footprints", sample.name + ".piq.TF-gene_interactions.tsv"), sep="\t", index=False)

    # Filter for TF-> TF interactions
    interactions_TF = interactions[interactions['gene'].isin(tfs)]
    interactions_TF.to_csv(os.path.join(data_dir, "footprints", sample.name + ".piq.TF-TF_interactions.tsv"), sep="\t", index=False)

    # Filter for nodes with more than 2 edges
    interactions_TF_filtered = interactions_TF[interactions_TF['interaction_score'] >= 1]
    interactions_TF_filtered.to_csv(os.path.join(data_dir, "footprints", sample.name + ".piq.TF-TF_interactions.filtered.tsv"), sep="\t", index=False)

# Network types:
# patient-specific:
# - build network specific to each patient
# - compare networks
# - cluster patients based on network connections

# for groups of patients:
# - come up with a way of combining signal from several patients from one group
# - build networks specific to groups


# Describe networks
for sample in samples:
    df = pd.read_csv(os.path.join(data_dir, "footprints", sample.name + ".piq.TF-gene_interactions.tsv"), sep="\t")

    G = nx.Graph()
    for i in df.index:
        G.add_edge(df.ix[i]['TF'], df.ix[i]['gene'], weight=df.ix[i]['interaction_score'])

    nx.shortest_path(G, 'PAX5', 'NFKB1', weight='weight')


# classify nodes into regulator/regulated
# color in visualization


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
