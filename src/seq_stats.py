#!/usr/bin/env python

"""
Determine number of mitochondrial reads.
Determine what percentage of the total number of duplicate reads is explained by the mitochondrial reads.
"""

import yaml
import os
from pipelines.models import Project, ATACseqSample
from pipelines import toolkit as tk
import textwrap
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

sns.set(style="whitegrid")


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
dataDir = os.path.join(config["paths"]["parent"], config["projectname"], "data")
resultsDir = os.path.join(config["paths"]["parent"], config["projectname"], "results", "plots")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample and s.cellLine == "CLL"]


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
df = pd.DataFrame(columns=["name", "total", "%duplicates", "total MT", "%dups nuclear", "%dups MT"])
for sample in samples:
    s = pd.Series(index=["name", "total", "%duplicates", "total MT", "%dups MT", "%dups nuclear"])

    dupsLog = os.path.join(sample.dirs.sampleRoot, sample.name + ".dupLog.txt")

    allDups = tk.parseDuplicateStats(sample.dupsMetrics)
    mtDups = tk.parseDuplicateStats(dupsLog)

    # is it empty?
    e = False if len(allDups) != 0 else True
    mtE = False if len(mtDups) != 0 else True

    s["name"] = sample.name if not e else None
    s["total"] = float(allDups["single-ends"]) + (float(allDups["paired-ends"]) * 2) if not e else None
    s["%duplicates"] = (float(allDups["duplicates"]) / s["total"]) * 100 if not e else None
    s["total MT"] = float(mtDups["single-ends"]) + (float(mtDups["paired-ends"]) * 2) if not mtE else None
    s["% MT"] = (s["total MT"] / s["total"]) * 100 if not mtE else None
    s["%dups MT"] = (float(mtDups["duplicates"]) / float(allDups["duplicates"])) * 100 if not mtE else None
    s["%dups nuclear"] = ((float(allDups["duplicates"]) - float(mtDups["duplicates"])) / s["total"]) * 100 if not mtE else None
    s['usable'] = s['total'] * (s['%duplicates'] / 100)
    s['%usable'] = (s['usable'] / s["total"]) * 100

    df = df.append(s, ignore_index=True)

#
df = df.replace(to_replace="None", value=np.nan)

# PLOT
# stripplot on a grid
g = sns.PairGrid(df.sort(["total"], ascending=False), x_vars=df.columns[1:], y_vars=["name"], size=15, aspect=.15)
g.map(sns.stripplot, size=10, orient="h", palette="Reds_r", edgecolor="gray")

for i, ax in enumerate(g.axes.flat):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # fix scales
    if i in [1, 3, 4, 5, 6]:
        ax.set_xlim(0, 100)

sns.despine(left=True, bottom=True)

# save plot
plt.savefig(os.path.join(resultsDir, "mitochondria_duplicates.pdf"), bbox_inches="tight")
