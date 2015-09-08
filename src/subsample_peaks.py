#!/usr/bin/env python

"""
"""

import os
from pipelines.models import Project
from pipelines import toolkit as tk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def subsample_bam(input_bam, n, output_bam):
    """
    Subsample n reads from bam file.
    """
    pass


def count_peaks_thresholds(peak_file, thresholds, output_file):
    """
    Count number of peaks depending on various q-value thresholds.

    :param thresholds: List with floats indicating various q-value thresholds.
    :type thresholds: list
    """
    pass


# Get path configuration
data_dir = os.path.join('.', "data")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# scratch dir (bam files and peaks)
scratch_dir = os.path.join("/scratch/users/arendeiro/", "subsample_peaks")
# dir for storing counts
output_dir = os.path.join(results_dir, "subsample_peaks")


# q-value thresholds to count at
thresholds = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.00001, 0.000005]

# M reads to subsample
millions = [1, 2, 4, 8, 16, 32, 64, 128]

# Loop through samples, make job
jobs = list()
for sample in prj.samples:
    # subsample bam file in log scale
    for reads in millions:
        # run name
        run_name = "_".join([sample.name, str(millions) + "M"])
        # make tmp bam file
        tmp_bam = os.path.join(scratch_dir, run_name + ".bam")
        # make tmp folder for peaks
        tmp_peaks_dir = os.path.join(scratch_dir, run_name)
        os.mkdir(tmp_peaks_dir)
        # peaks
        peak_file = os.path.join(tmp_peaks_dir, run_name + "_peaks.narrowPeak")
        # count output
        output_file = os.path.join(output_dir, run_name + ".tsv")

        # slurm header
        cmd = tk.slurmHeader(
            jobName="_".join([sample.name, str(millions)]),
            output=os.path.join(scratch_dir, run_name + ".log")
        )

        # subsample
        subsample_bam(sample.mapped, reads, tmp_bam)

        # call peaks
        tk.macs2CallPeaksATACSeq(tmp_bam, tmp_peaks_dir, run_name, sample.genome)

        # count peaks
        count_peaks_thresholds(peak_file, thresholds, output_file)

        # slurm footer
        cmd += tk.slurmFooter()

        jobs.append(cmd)

# submit jobs
for job in jobs:
    tk.slurmSubmitJob(job)

# collect

# plot
