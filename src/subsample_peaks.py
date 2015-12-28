#!/usr/bin/env python

"""
Plot number of peaks per reads subsampled at various q-value thresholds.
"""

import os
from pipelines.models import Project
from pipelines import toolkit as tk
import numpy as np
import textwrap
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

sns.set_style("whitegrid")
sns.set_context("paper")


def subsample_bam(input_bam, n, output_bam):
    """
    Subsample n reads from bam file.
    """
    cmd = """
    macs2 randsample -t {0} -n {1} -o {2}
    """.format(input_bam, n, output_bam)
    return cmd


def count_peaks_thresholds(peak_file, thresholds, output_file):
    """
    Count number of peaks depending on various q-value thresholds.

    :param thresholds: List with floats indicating various q-value thresholds.
    :type thresholds: list
    """
    cmds = list()
    for threshold in thresholds:
        threshold = -np.log10(threshold)
        cmd = """
    COUNT=`awk '{ if ($8 >= %f) print }' %s | wc -l`; echo "%f\t$COUNT" >> %s
        """ % (threshold, peak_file, threshold, output_file)
        cmds.append(cmd)

    return "".join(cmds)


# Get path configuration
data_dir = os.path.join('.', "data")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# scratch dir (bam files and peaks)
scratch_dir = os.path.join("/scratch/users/arendeiro/", "subsample_peaks")
if not os.path.exists(scratch_dir):
    os.mkdir(scratch_dir)
# dir for storing counts
output_dir = os.path.join(results_dir, "subsample_peaks")
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# q-value thresholds to count at
thresholds = [1e-2]  # , 1e-4, 1e-8, 1e-16, 1e-32]

# M reads to subsample
millions = range(1, 33)

# Loop through samples, make job
jobs = list()
for sample in prj.samples:
    if sample.technique == "ATAC-seq" and sample.cellLine == "CLL":
        # subsample bam file in log scale
        for reads in millions:
            # run name
            run_name = "_".join([sample.name, str(reads) + "M"])
            # job file
            job_file = os.path.join(scratch_dir, run_name + ".slurmjob.sh")
            # make tmp bam file
            tmp_bam = os.path.join(scratch_dir, run_name + ".bam")
            # make tmp folder for peaks
            tmp_peaks_dir = os.path.join(scratch_dir, run_name)
            if not os.path.exists(tmp_peaks_dir):
                os.mkdir(tmp_peaks_dir)
            # peaks
            peak_file = os.path.join(tmp_peaks_dir, run_name + "_peaks.narrowPeak")
            # count output
            output_file = os.path.join(output_dir, run_name + ".tsv")

            # slurm header
            cmd = tk.slurmHeader(
                jobName=run_name,
                output=os.path.join(scratch_dir, run_name + ".log"),
                queue="shortq",
                cpusPerTask=4
            )

            # subsample
            cmd += subsample_bam(sample.mapped, reads * 1e6, tmp_bam)

            # call peaks
            cmd += tk.macs2CallPeaksATACSeq(tmp_bam, tmp_peaks_dir, run_name, sample.genome)

            # count peaks
            cmd += count_peaks_thresholds(peak_file, thresholds, output_file)

            # slurm footer
            cmd += "\n"
            cmd += tk.slurmFooter()

            # write job to file
            with open(job_file, 'w') as handle:
                handle.writelines(textwrap.dedent(cmd))

            # append file to jobs
            jobs.append(job_file)

# submit jobs
for job in jobs:
    tk.slurmSubmitJob(job)

# collect
counts = pd.DataFrame()

for sample in prj.samples:
    if sample.technique == "ATAC-seq" and sample.cellLine == "CLL":
        # subsample bam file in log scale
        for reads in millions:
            run_name = "_".join([sample.name, str(reads) + "M"])
            output_file = os.path.join(output_dir, run_name + ".tsv")

            # read outfile
            try:
                with open(output_file, "r") as handle:
                    lines = handle.readlines()
            except:
                continue

            # parse
            for line in lines:
                threshold, count = line.strip().split("\t")
                threshold = int(float(threshold))
                count = int(count)
                series = pd.Series([reads, threshold, count], index=['reads', 'threshold', 'peak_count'])
                counts = counts.append(series, ignore_index=True)

cc = counts[
    (counts['peak_count'] > 100) &
    (counts['threshold'] == 2)
]

# save
cc.to_csv(os.path.join(results_dir, "read_subsampling_peak_number.csv"), index=False)

# average over samples
ccc = cc.groupby(['reads']).apply(np.mean)
plt.plot(range(1 + len(ccc['peak_count'])), [0] + ccc['peak_count'].tolist(), 'o-')

ccc = cc.groupby(['reads']).apply(lambda x: pd.Series(stats.norm.interval(0.95, loc=x.mean(), scale=x.std() / np.sqrt(len(x)))))

plt.plot(range(1 + len(ccc[0])), [0] + [x[0] for x in ccc[0]], 'o-')
plt.plot(range(1 + len(ccc[1])), [0] + [x[0] for x in ccc[1]], 'o-')
plt.savefig(os.path.join(plots_dir, "read_subsampling_peak_number.95CI.svg"))
