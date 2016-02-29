#!/usr/bin/env python

"""
Determine number of mitochondrial reads.
Determine what percentage of the total number of duplicate reads is explained by the mitochondrial reads.
"""

import os
import re
from pipelines.models import Project
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


def parse_duplicate_stats(stats_file):
    """
    Parses sambamba markdup output, returns series with values.

    :param stats_file: sambamba output file with duplicate statistics.
    :type stats_file: str
    """
    import pandas as pd
    import re

    series = pd.Series()

    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return series

    try:
        line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
        series["single-ends"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
        series["paired-ends"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
        line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
        series["duplicates"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
    except IndexError:
        pass
    return series


def parse_FRiP(frip_file, series):
    """
    Calculates the fraction of reads in peaks for a given sample.

    :param sample: A Sample object with the "peaks" attribute.
    :type sampleName: pipelines.Sample
    """
    import re
    import pandas as pd

    try:
        with open(frip_file, "r") as handle:
            content = handle.readlines()
    except:
        return pd.np.nan

    if content[0].strip() == "":
        return pd.np.nan

    readsInPeaks = int(re.sub("\D", "", content[0]))
    mappedReads = float(series["total_reads"])  # - series["unaligned"]

    return readsInPeaks / mappedReads


def calculateFRiP(inputBam, inputBed, output):
    cmd = """
    cut -f 1,2,3 {0} |""".format(inputBed)
    cmd += " bedtools coverage -counts -abam {0} -b - |".format(inputBam)
    cmd += """ awk '{{sum+=$4}} END {{print sum}}' > {0}
    """.format(output)

    return cmd


def parse_rna_stats(stats_file):
    """
    Parses sambamba markdup output, returns series with values.

    :param stats_file: sambamba output file with duplicate statistics.
    :type stats_file: str
    """
    import pandas as pd

    try:
        return pd.read_csv(stats_file, sep="\t", header=None).set_index(0)[1]
    except:
        return pd.Series()

# Get path configuration
data_dir = os.path.join(".", "data")
results_dir = os.path.join(".", "results")
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("metadata/project_config.yaml")
prj.add_sample_sheet()

samples = [s for s in prj.samples if s.library != "RNA-seq" and s.cell_line == "CLL"]
samples = [s for s in prj.samples if s.library == "ChIPmentation" and s.cell_line == "CLL"]


# GET DUPLICATES IN MITOCHONDRIAL GENOME ONLY
# Submit job marking duplicates for each sample
for sample in samples:
    dups_log = os.path.join(sample.paths.sample_root, sample.name + ".dupLog.txt")

    cmd = tk.slurmHeader("_".join([sample.name, "mitochondria_duplicates"]), "/home/arendeiro/scratch/mitoLog.txt", cpusPerTask=4)
    cmd += """
    sambamba index -t 4 {0}
    """.format(sample.mapped)
    cmd += """
    sambamba slice {0} chrM | sambamba markdup -t 4 /dev/stdin /home/arendeiro/scratch/{1}.dups.rmMe 2>  {2}
    """.format(sample.mapped, sample.name, dups_log)
    cmd += tk.slurmFooter()

    jobFile = "/home/arendeiro/scratch/" + sample.name + ".slurm.sh"

    with open(jobFile, "w") as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(jobFile)

r = list()
# Collect duplicate counts, calculate other variables
df = pd.DataFrame(columns=["name", "total_reads", "%duplicates", "total MT", "%dups nuclear", "%dups MT"])
for sample in prj.samples:
    s = pd.Series(index=["name", "total_reads", "%duplicates", "total MT", "%dups MT", "%dups nuclear"])

    if sample.library != "RNA-seq":
        if sample.library == "ATAC-seq":
            all_dups_log = re.sub("dups_metrics", "duplicates", sample.dups_metrics)
        else:
            all_dups_log = sample.dups_metrics
        mt_dups_log = os.path.join(sample.paths.sample_root, sample.name + ".dupLog.txt")

        allDups = parse_duplicate_stats(all_dups_log)
        mtDups = parse_duplicate_stats(mt_dups_log)

        # is it empty?
        e = False if len(allDups) != 0 else True
        mtE = False if len(mtDups) != 0 else True

        if e or mtE:
            print(sample.name, e, mtE)
            continue

        s["name"] = sample.name if not e else None
        s["total_reads"] = float(allDups["single-ends"]) + (float(allDups["paired-ends"]) * 2) if not e else None
        s["duplicates"] = float(allDups["duplicates"]) if not e else None
        s["%duplicates"] = ((s["duplicates"] / s["total_reads"]) * 100) if not e else None
        s["total MT"] = (float(mtDups["single-ends"]) + (float(mtDups["paired-ends"]) * 2)) if not mtE else None
        s["% MT"] = (s["total MT"] / s["total_reads"]) * 100 if not mtE else None
        s["%dups MT"] = (float(mtDups["duplicates"]) / float(allDups["duplicates"]) * 100) if not mtE else None
        s["%dups nuclear"] = ((float(allDups["duplicates"]) - float(mtDups["duplicates"])) / s["total_reads"]) * 100 if not mtE else None
        s["usable"] = s["total_reads"] - s["duplicates"] if not e else None
        s["%usable"] = (s["usable"] / s["total_reads"]) * 100 if not e else None

        s["frip"] = parse_FRiP(sample.frip, s)

        if pd.isnull(s["frip"]):
            r.append(sample)

    else:
        s2 = parse_rna_stats(os.path.join(sample.paths.sample_root, sample.name + ".rnaBitSeq_stats.tsv"))
        s["name"] = sample.name
        s["total_reads"] = float(s2["Fastq_reads"])
        s["usable"] = float(s2["Deduplicated_reads"])
        s["duplicates"] = s["total_reads"] - s["usable"]
        s["%duplicates"] = ((s["duplicates"] / s["total_reads"]) * 100)
        s["%usable"] = (s["usable"] / s["total_reads"]) * 100

    # add to dataframe
    df = df.append(s, ignore_index=True)

# replace None with np.nan
df = df.replace(to_replace="None", value=np.nan)

# order columns
df = df[["name", "total_reads", "duplicates", "%duplicates", "total MT", "% MT", "%dups MT", "%dups nuclear", "usable", "%usable", "frip"]]

# save table
df.to_csv(os.path.join(plots_dir, "seq_stats.mitochondria_duplicates.csv"), index=False)

# PLOT
df2 = df[df.name.str.contains("ATAC-seq")].sort(["total_reads"], ascending=False)
# stripplot on a grid
g = sns.PairGrid(df2, x_vars=df2.columns[1:], y_vars=["name"], size=30, aspect=.15)
g.map(sns.stripplot, size=10, orient="h", palette="Reds_r", edgecolor="gray")

for i, ax in enumerate(g.axes.flat):
    # Make the grid horizontal instead of vertical
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    # fix scales
    if i in [2, 4, 5, 6, 8]:
        ax.set_xlim(0, 100)

    if i == 9:
        ax.set_xlim(0, 0.4)

sns.despine(left=True, bottom=True)

# save plot
plt.savefig(os.path.join(plots_dir, "seq_stats.mitochondria_duplicates.svg"), bbox_inches="tight")
