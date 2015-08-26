#!/usr/bin/env python

import os
import pandas as pd
import pybedtools


def get_tss(entry):
    if entry['strand'] == "+":
        entry['end'] = entry['start'] + 1
    elif entry['strand'] == "-":
        entry['start'] = entry['end'] - 1
    return entry


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


# Get roadmap CD19 perypheral blood HMM state annotation
# read more about it here http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
roadmap_url = "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"
roadmap_15statesHMM = "all.mnemonics.bedFiles.tgz"
roadmap_15statesHMM_CD19 = "E032_15_coreMarks_mnemonics.bed.gz"

os.system("wget {0}{1}".format(roadmap_url, roadmap_15statesHMM))
os.system("tar zxvf {0} {1}".format(roadmap_15statesHMM, roadmap_15statesHMM_CD19))
os.system("gzip -d {0}".format(roadmap_15statesHMM_CD19))
os.system("mv E032_15_coreMarks_mnemonics.bed ../data/E032_15_coreMarks_mnemonics.bed")
