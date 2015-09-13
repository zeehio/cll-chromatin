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
import numpy as np
import re
import cPickle as pickle
import pandas as pd
import pybedtools


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
    for motif in range(n_motifs):
        cmd = "Rscript ~/workspace/piq-single/pwmmatch.exact.r"
        cmd += " ~/workspace/piq-single/common.r"
        cmd += " {0}".format(motif_file)
        cmd += " " + str(motif)
        cmd += " /scratch/users/arendeiro/piq/motif.matches/"
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


def piq_footprint(bam_cache, n_motifs, tmp_dir, results_dir):
    """
    """
    cmds = list()
    for motif in range(1, n_motifs + 1):
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


def piq_to_network(results_dir, n_motifs):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """
    # list results_dir
    files = os.listdir(results_dir)
    # get all cll peaks to filter data
    all_peaks = pybedtools.BedTool("data/cll_peaks.bed")
    # read in gene info
    tsss = pd.read_csv("data/hg19.refSeq.TSS.mRNA.bed", sep="\t", header=None)
    tsss.columns = ["chrom", "start", "end", "id"]

    # prepare TF vs Gene matrix
    scores = pd.DataFrame(index=tsss["id"], columns=range(1, n_motifs + 1))

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in range(1, n_motifs + 1):
        print(motif)
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
        if len(footprints) < 1:
            print "continuing 1"
            continue

        footprints[['chr', 'start', 'end', 'pwm', 'shape', 'score', 'purity']].to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # filter for motifs overlapping CLL peaks
        footprints = pybedtools.BedTool(os.path.join("tmp.bed")).intersect(all_peaks, wa=True).to_dataframe()
        footprints.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity"]

        # If empty give 0 to every gene for this TF
        if len(footprints) < 1:
            print "continuing 2"
            continue

        # CONNECT
        # Now assign a score between this TF and every gene:
        # get distance to nearest gene TSS in the same chromosome as footprint
        for chrom in footprints["chrom"].unique():
            # calculate the distance between each footprint and every gene in the chromosome
            for i in tsss[tsss["chrom"] == chrom].index:
                gene_scores = list()
                for j in footprints[footprints["chrom"] == chrom].index:
                    dist = abs(footprints.ix[j]["start"] - tsss.ix[i]["start"])
                    gene_scores.append(2 * (footprints.ix[j]["purity"] - 0.5) * 10 ** -(dist / 1e6))
                scores.loc[tsss.ix[i]["id"], motif] = sum(gene_scores)

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

prj.samples = annotate_igvh_mutations(prj.samples, clinical)

to_exclude_sample_id = ['1-5-45960']

# FOOTPRINTING
motifs_file = "~/workspace/piq-single/pwms/jasparfix.txt"
n_motifs = 1316

# prepare motifs for footprinting (done once)
# cmds = piq_prepare_motifs(motifs_file, n_motifs)
# for cmd in cmds:
#     os.system(cmd)

# stupid PIQ hard-coded links
os.chdir("/home/arendeiro/workspace/piq-single/")

# for each sample create R cache with bam file
jobs = list()
for sample in prj.samples:
    if sample.sampleID in to_exclude_sample_id or sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(sample.dirs.sampleRoot, "footprints")
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)
    r_data = os.path.join(foots_dir, sample.name + ".filteredshifted.RData")
    tmp_dir = os.path.join(scratch_dir, sample.name)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    job_file = os.path.join(foots_dir, "slurm_job.sh")

    # prepare slurm job header
    cmd = tk.slurmHeader(sample.name + "_PIQ_footprinting", os.path.join(foots_dir, "slurm.log"), cpusPerTask=2, queue="shortq")

    # stupid PIQ hard-coded links
    cmd += """
    cd /home/arendeiro/workspace/piq-single/
    """

    # prepare bams
    cmd += piq_prepare_bams([sample.filteredshifted], r_data)

    # footprint
    cmd_list = piq_footprint(r_data, n_motifs, tmp_dir, results_dir=foots_dir)
    cmd += "\n".join(cmd_list)

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


# for each sample launch several jobs (>1000) to footprint
jobs = list()
for sample in prj.samples[1:5]:
    if sample.sampleID in to_exclude_sample_id or sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(sample.dirs.sampleRoot, "footprints")
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)
    r_data = os.path.join(foots_dir, sample.name + ".filteredshifted.RData")
    tmp_dir = os.path.join(scratch_dir, sample.name)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    cmd_list = piq_footprint(r_data, n_motifs, tmp_dir, results_dir=foots_dir)

    for motif, foot_cmd in enumerate(cmd_list):
        motif += 1
        job_file = os.path.join(foots_dir, "slurm_job_motif%i.sh" % motif)
        slurm_log = os.path.join(foots_dir, "slurm_motif%i.log" % motif)

        # prepare slurm job header
        cmd = tk.slurmHeader(sample.name + "_PIQ_footprinting_motif%i" % motif, slurm_log, cpusPerTask=2, queue="shortq")

        # stupid PIQ hard-coded links
        cmd += """
        cd /home/arendeiro/workspace/piq-single/
        """

        # footprint
        cmd += """
        {0}
        """.format(foot_cmd)

        # delete its own slurm files
        cmd += """
        rm {0}
        rm {1}
        """.format(job_file, slurm_log)

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


# FOOTPRINT GROUPS UNMUTATED/MUTATED

# get unmutated/mutated samples
muts = list()
unmuts = list()
for sample in prj.samples:
    if sample.sampleID in to_exclude_sample_id or sample.technique != "ATAC-seq" or sample.cellLine != "CLL" or sample.readType == "SE":
        continue
    if sample.mutated is True:
        muts.append(sample.filteredshifted)
    elif sample.mutated is False:
        unmuts.append(sample.filteredshifted)

# merge bam files manually
muts_bam = os.path.join(data_dir, "CLL_all_igvhmutated_samples.filteredshifted.bam")
muts_cache = os.path.join(data_dir, "CLL_all_igvhmutated_samples.filteredshifted.RData")
unmuts_bam = os.path.join(data_dir, "CLL_all_igvhunmutated_samples.filteredshifted.bam")
unmuts_cache = os.path.join(data_dir, "CLL_all_igvhunmutated_samples.filteredshifted.RData")

os.system("sambamba merge -t 12 -p %s " % muts_bam + " ".join(muts))
os.system("sambamba merge -t 12 -p %s " % unmuts_bam + " ".join(unmuts))

# prepare merged bam files from IGVH mutated and unmutated samples
cmd = piq_prepare_bams(muts, muts_cache)
os.system(cmd)
cmd = piq_prepare_bams(unmuts, unmuts_cache)
os.system(cmd)

# prepare also Nextera background
nextera_dir = "/scratch/users/arendeiro/PGA1Nextera"
nextera_bam = os.path.join(nextera_dir, "PGA_0001_Nextera-2.bam")
nextera_cache = os.path.join(nextera_dir, "PGA_0001_Nextera-2.RData")

cmd = piq_prepare_bams([nextera_bam], nextera_cache)
os.system(cmd)

# run footprinting
os.mkdir(os.path.join(scratch_dir, "mutated"))
cmds = piq_footprint(
    os.path.join(data_dir, "CLL_all_igvhmutated_samples.filteredshifted.RData"), n_motifs, scratch_dir,
    results_dir=os.path.join(scratch_dir, "mutated")
)
for cmd in cmds:
    os.system(cmd)
os.mkdir(os.path.join(scratch_dir, "unmutated"))
cmds = piq_footprint(os.path.join(data_dir, "CLL_all_igvhunmutated_samples.filteredshifted.RData"), n_motifs, scratch_dir, os.path.join(scratch_dir, "unmutated"))
for cmd in cmds:
    os.system(cmd)


# parse output,
# connect each motif to a gene
for sample in prj.samples[1:5]:
    print sample
    scores = piq_to_network(os.path.join(sample.dirs.sampleRoot, "footprints"), n_motifs)
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
