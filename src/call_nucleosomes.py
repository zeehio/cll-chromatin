#!/usr/bin/env python

"""
Call nucleosomes with NucleoATAC
"""

import os
from pipelines.models import Project, ATACseqSample
from pipelines import toolkit as tk
import textwrap


def nucleoatac(bed_file, bam_file, output_dir, fasta):
    cmd = """
    nucleoatac run --bed {0} --bam {1} --out {2} --fasta {3} --cores 8
    """.format(bed_file, bam_file, output_dir, fasta)

    return cmd


# Get path configuration
data_dir = os.path.abspath(os.path.join('.', "data"))
scratch_dir = os.path.join("/scratch/users/arendeiro/piq")
results_dir = os.path.abspath(os.path.join('.', "results"))
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Select ATAC-seq samples
prj.samples = [s for s in prj.samples if type(s) == ATACseqSample]

to_exclude_sample_id = ['1-5-45960']

# static info
fasta = "/data/groups/lab_bock/shared/resources/genomes/hg19/hg19.fa"
peaks = "/home/arendeiro/cll-patients/data/cll_peaks.bed"

# for each sample create R cache with bam file
jobs = list()
for sample in prj.samples:
    if sample.sampleID in to_exclude_sample_id or sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    nucleoatac_dir = os.path.join(sample.dirs.sampleRoot, "nucleoatac")
    if not os.path.exists(nucleoatac_dir):
        os.mkdir(nucleoatac_dir)
    r_data = os.path.join(nucleoatac_dir, sample.name + ".filteredshifted.RData")
    tmp_dir = os.path.join(scratch_dir, sample.name)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    job_file = os.path.join(nucleoatac_dir, "slurm_job.sh")

    # prepare slurm job header
    cmd = tk.slurmHeader(sample.name + "_NucleoATAC_nucleosome_calling", os.path.join(nucleoatac_dir, "slurm.log"), cpusPerTask=8, queue="longq")

    # footprint
    cmd += nucleoatac(peaks, sample.filtered, nucleoatac_dir, fasta) + "\n"

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    # append file to jobs
    jobs.append(job_file)


# submit jobs (to create bam file caches)
for job in jobs:
    tk.slurmSubmitJob(job)
