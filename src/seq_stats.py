#!/usr/bin/env python

"""
Determine number of mitochondrial reads.
Determine what percentage of the total number of duplicate reads is explained by the mitochondrial reads.
"""

import os
from pipelines.models import Project, ATACseqSample
from pipelines import toolkit as tk
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

sns.set(style="whitegrid")
pd.set_option("date_dayfirst", True)
sns.set_context("poster")
sns.set_palette(sns.color_palette("colorblind"))

import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('text', usetex='false')


# Get path configuration
data_dir = os.path.join(".", "data")
results_dir = os.path.join(".", "results")
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")


# Select ATAC-seq samples
samples_to_exclude = [
    "CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19",
    "CLL_ATAC-seq_5186_1-5-57350_ATAC17-4_hg19",
    "CLL_ATAC-seq_4784_1-5-52817_ATAC17-6_hg19",
    "CLL_ATAC-seq_981_1-5-42480_ATAC16-6_hg19",
    "CLL_ATAC-seq_5277_1-5-57269_ATAC17-8_hg19",
    "CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19",
    "CLL_ATAC-seq_5147_1-5-48105_ATAC17-2_hg19",
    "CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19"]

samples = [s for s in prj.samples if type(s) == ATACseqSample and s.cellLine == "CLL" and s.name not in samples_to_exclude]


# GET DUPLICATES IN MITOCHONDRIAL GENOME ONLY
# Submit job marking duplicates for each sample
for sample in samples:
    dupsLog = os.path.join(sample.dirs.sampleRoot, sample.name + ".dupLog.txt")

    cmd = tk.slurmHeader("_".join([sample.name, "mitochondria_duplicates"]), "/home/arendeiro/scratch/mitoLog.txt", cpusPerTask=4)
    cmd += """sambamba slice {0} chrM | sambamba markdup -t 4 /dev/stdin /home/arendeiro/scratch/{1}.dups.rmMe 2>  {2}\n""".format(sample.mapped, sample.name, dupsLog)
    cmd += tk.slurmFooter()

    jobFile = "/home/arendeiro/scratch/" + sample.name + ".slurm.sh"

    with open(jobFile, "w") as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(jobFile)

# Collect duplicate counts, calculate other variables
df = pd.DataFrame(columns=["name", "total_reads", "%duplicates", "total MT", "%dups nuclear", "%dups MT"])
for sample in samples:
    s = pd.Series(index=["name", "total_reads", "%duplicates", "total MT", "%dups MT", "%dups nuclear"])

    dupsLog = os.path.join(sample.dirs.sampleRoot, sample.name + ".dupLog.txt")

    allDups = tk.parseDuplicateStats(sample.dupsMetrics)
    mtDups = tk.parseDuplicateStats(dupsLog)

    # is it empty?
    e = False if len(allDups) != 0 else True
    mtE = False if len(mtDups) != 0 else True

    s["name"] = sample.name if not e else None
    s["total_reads"] = float(allDups["single-ends"]) + (float(allDups["paired-ends"]) * 2) if not e else None
    s["duplicates"] = float(allDups["duplicates"])
    s["%duplicates"] = (s["duplicates"] / s["total_reads"]) * 100 if not e else None
    s["total MT"] = float(mtDups["single-ends"]) + (float(mtDups["paired-ends"]) * 2) if not mtE else None
    s["% MT"] = (s["total MT"] / s["total_reads"]) * 100 if not mtE else None
    s["%dups MT"] = (float(mtDups["duplicates"]) / float(allDups["duplicates"])) * 100 if not mtE else None
    s["%dups nuclear"] = ((float(allDups["duplicates"]) - float(mtDups["duplicates"])) / s["total_reads"]) * 100 if not mtE else None
    s["usable"] = s["total_reads"] - s["duplicates"]
    s["%usable"] = (s["usable"] / s["total_reads"]) * 100

    df = df.append(s, ignore_index=True)

#
df = df.replace(to_replace="None", value=np.nan)

# order columns
df = df[["name", "total_reads", "duplicates", "%duplicates", "total MT", "% MT", "%dups MT", "%dups nuclear", "usable", "%usable"]]

# save table
df.to_csv(os.path.join(plots_dir, "mitochondria_duplicates.csv"), index=False)

# PLOT
# stripplot on a grid
g = sns.PairGrid(df.sort(["total_reads"], ascending=False), x_vars=df.columns[1:], y_vars=["name"], size=15, aspect=.15)
g.map(sns.stripplot, size=10, orient="h", palette="Reds_r", edgecolor="gray")

for i, ax in enumerate(g.axes.flat):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # fix scales
    if i in range(1, 7):
        ax.set_xlim(0, 100)

sns.despine(left=True, bottom=True)

# save plot
plt.savefig(os.path.join(plots_dir, "mitochondria_duplicates.svg"), bbox_inches="tight")
