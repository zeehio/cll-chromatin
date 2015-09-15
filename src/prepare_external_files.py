#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import pybedtools


def get_tss(entry):
    if entry['strand'] == "+":
        entry['end'] = entry['start'] + 1
    elif entry['strand'] == "-":
        entry['start'] = entry['end'] - 1
    return entry


def sra2bam(sra_acession, output_bam):
    from pipelines import toolkit as tk
    import textwrap
    # Slurm header
    cmd = tk.slurmHeader("_".join(["sra2bam", sra_acession]), "/scratch/users/arendeiro/%s_sra2bam.log" % sra_acession, cpusPerTask=4)

    # SRA to BAM
    cmd += "\n    sam-dump {0} | sambamba view -t 4 -S -f bam /dev/stdin > /home/arendeiro/cll-patients/data/external/{0}.bam\n".format(sra_acession)
    # Slurm footer
    cmd += tk.slurmFooter() + "\n"

    # Write job to file
    job_file = "/scratch/users/arendeiro/%s_sra2bam.sh" % sra_acession
    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurmSubmitJob(job_file)


# Get encode blacklisted regions
blacklist = "wgEncodeDacMapabilityConsensusExcludable.bed.gz"

os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/{0}".format(blacklist))
os.system("gzip -d {0}".format(blacklist))
os.system("mv wgEncodeDacMapabilityConsensusExcludable.bed ../data/wgEncodeDacMapabilityConsensusExcludable.bed")


# # Get ensembl genes and transcripts from grch37 and hg19 chromosomes (to get TSSs)
# ensembl_genes = "ensGene.txt.gz"
# os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/{0}".format(ensembl_genes))
# os.system("gzip -d {0}".format(ensembl_genes))
# os.system("mv ensGene.txt ../data/ensGene.txt")
# # get names
# ensembl_names = "ensemblToGeneName.txt.gz"
# os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/{0}".format(ensembl_names))
# os.system("gzip -d {0}".format(ensembl_names))
# os.system("mv ensemblToGeneName.txt ../data/ensemblToGeneName.txt")

# Get ensembl genes from grch37 and hg19 chromosomes from UCSC (you get introns this way)
regions = [
    "./data/ensembl_tss.bed",
    "./data/ensembl_tss2kb.bed",
    "./data/ensembl_utr5.bed",
    "./data/ensembl_exons.bed",
    "./data/ensembl_introns.bed",
    "./data/ensembl_utr3.bed",
]

# remove annotation after Ensembl transcriptID
for region in regions:
    r = pybedtools.BedTool(region).sort()
    r = pd.read_csv(region, sep="\t", header=None)
    r[3] = r[3].apply(lambda x: x.split("_")[0])
    r.to_csv(region, sep="\t", header=None, index=False)

for i, region in enumerate(regions[1:]):
    r = pybedtools.BedTool(region)
    if i == 0:
        genes = r
    else:
        genes.cat(r)
genes.sort().saveas("./data/ensembl_genes.bed")

# Make bed file
genes = pd.read_csv("../data/ensGene.txt", sep="\t", header=None)
genes = genes[[2, 4, 5, 12, 1, 3]]
genes.columns = ['chrom', 'start', 'end', 'gene', 'transcript', 'strand']

# Annotate with gene names
names = pd.read_csv("../data/ensemblToGeneName.txt", sep="\t", header=None)
names.columns = ['transcript', 'name']
annotation = pd.merge(genes, names)
annotation.to_csv("../data/GRCh37_hg19_ensembl_genes.bed", sep="\t", index=False, header=False)

# Get TSSs
tsss = annotation.apply(get_tss, axis=1)
tsss.to_csv("../data/GRCh37_hg19_ensembl_genes.tss.bed", sep="\t", index=False, header=False)
# sort bed file
tsss = pybedtools.BedTool("../data/GRCh37_hg19_ensembl_genes.tss.bed")
tsss.sort().saveas("../data/GRCh37_hg19_ensembl_genes.tss.bed")


# Get refseq bed file from UCSC, add 1 bp upstream, name as hg19.refSeq.TSS.bed
"sed 's/_up_1_.*//' hg19.refSeq.TSS.bed > t"
# Filter out ncRNAs
"grep NM t > hg19.refSeq.TSS.mRNA.bed"


# Get roadmap CD19 perypheral blood HMM state annotation
# read more about it here http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
roadmap_url = "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"
roadmap_15statesHMM = "all.mnemonics.bedFiles.tgz"
roadmap_15statesHMM_CD19 = "E032_15_coreMarks_mnemonics.bed.gz"

os.system("wget {0}{1}".format(roadmap_url, roadmap_15statesHMM))
os.system("tar zxvf {0} {1}".format(roadmap_15statesHMM, roadmap_15statesHMM_CD19))
os.system("gzip -d {0}".format(roadmap_15statesHMM_CD19))
os.system("mv E032_15_coreMarks_mnemonics.bed ../data/E032_15_coreMarks_mnemonics.bed")


# Get GOtermID - GOtermDescription mapping
# download data from Biomart (for some reason the GOterm name cannot be get automatically)
# http://www.ensembl.org/biomart/martview/6451fcd5296302994382deee7bd9c8eb?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.name_1006|hsapiens_gene_ensembl.default.feature_page.go_id&FILTERS=&VISIBLEPANEL=resultspanel
"mv x data/goID_goName.csv"


# Puente 2015
# Get WGS variants
cmds = list()
cmds.append("wget http://www.nature.com/nature/journal/vaop/ncurrent/extref/nature14666-s3.zip")
cmds.append("unzip nature14666-s3.zip")
supFile = "data/Puente2015.supplementaryTable2.csv"
cmds.append("mv supplementaryTable2.tsv %s" % supFile)
for cmd in cmds:
    os.system(cmd)

# Read in data
df = pd.read_csv(supFile)
df.rename(columns={"#CHROM": "CHROM"}, inplace=True)
df.drop(["WG", "ID"], inplace=True, axis=1)  # this are uninformative columns -> unique(col) == 1

# Filter for SNPs
df2 = df[
    (np.array([len(x) for x in df['REF']]) == 1) &
    (np.array([len(x) for x in df['ALT']]) == 1)
]

# Collapse across patients
colapsed = df2.groupby(['CHROM', 'POS', 'REF', 'ALT']).aggregate(lambda x: ",".join([str(i) for i in x]))

# Filter for recurrent genes
# get extended data table 2
# interesect


# External samples, SRA to unaligned bam

samples = {
    "GM12878": "SRR891268",  # rep1
    "CD4_T_day1": "SRR891275",  # 1hour - rep1
    "CD4_T_day2": "SRR891277",  # 2hour - rep1
    "CD4_T_day3": "SRR891279",  # 3hour - rep1
    "Raji": "SRR1787814",
    "RJ2.2.5": "SRR1787816",
    "SKNMC": "SRR1594026",
    "Monocyte-derived_dendritic_cells": "SRR1725732",
    "Monocyte-derived_dendritic_cells_infected": "SRR1725731",
    # "IMR90": [
    #     "SRR1448792",
    #     "SRR1448793"
    # ],
    # "IMR90_Nutlin-3a": [
    #     "SRR1448794",
    #     "SRR1448795"
    # ]
}

for accession in samples.values():
    sra2bam(accession, os.path.join("data/external/%s.bam" % accession))
