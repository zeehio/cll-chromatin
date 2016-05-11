#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import pybedtools
import re
from Bio import motifs


def get_tss(x):
    if x['strand'] == "+":
        x['end'] = x['start'] + 1
    elif x['strand'] == "-":
        x['start'] = x['end'] - 1
    return x


def get_promoter(x, radius=2500):
    if x['strand'] == "+":
        x['start'] = x['start'] - radius if x['start'] - radius > 0 else 0
        x['end'] = x['start'] + (radius * 2)
        return x
    elif x['strand'] == "-":
        x['end'] += radius
        x['start'] = x['end'] - (radius * 2) if x['end'] - (radius * 2) > 0 else 0
        return x


def get_promoter_and_genebody(x, radius=2500):
    if x['strand'] == "+":
        x['start'] = x['start'] - radius if x['start'] - radius > 0 else 0
        return x
    elif x['strand'] == "-":
        x['end'] += radius
        return x


def bowtie2Map(inputFastq1, outputBam, log, metrics, genomeIndex, maxInsert, cpus, inputFastq2=None):
    import re

    outputBam = re.sub("\.bam$", "", outputBam)
    # Admits 2000bp-long fragments (--maxins option)
    cmd = "bowtie2 --very-sensitive -p {0}".format(cpus)
    cmd += " -x {0}".format(genomeIndex)
    cmd += " --met-file {0}".format(metrics)
    if inputFastq2 is None:
        cmd += " {0} ".format(inputFastq1)
    else:
        cmd += " --maxins {0}".format(maxInsert)
        cmd += " -1 {0}".format(inputFastq1)
        cmd += " -2 {0}".format(inputFastq2)
    cmd += " 2> {0} | samtools view -S -b - | samtools sort - {1}".format(log, outputBam)

    return cmd


# Get encode blacklisted regions
blacklist = "wgEncodeDacMapabilityConsensusExcludable.bed.gz"

os.system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/{0}".format(blacklist))
os.system("gzip -d {0}".format(blacklist))
os.system("mv {} ../data/{}".format(re.sub(".gz", "", blacklist)))


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
    "data/external/ensembl_tss.bed",
    "data/external/ensembl_tss2kb.bed",
    "data/external/ensembl_utr5.bed",
    "data/external/ensembl_exons.bed",
    "data/external/ensembl_introns.bed",
    "data/external/ensembl_utr3.bed",
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
genes.sort().saveas("data/external/ensembl_genes.bed")

# Make bed file
genes = pd.read_csv("data/ensGene.txt", sep="\t", header=None)
genes = genes[[2, 4, 5, 12, 1, 3]]
genes.columns = ['chrom', 'start', 'end', 'gene', 'transcript', 'strand']

# Annotate with gene names
names = pd.read_csv("data/external/ensemblToGeneName.txt", sep="\t", header=None)
names.columns = ['transcript', 'name']
annotation = pd.merge(genes, names)
annotation.to_csv("data/external/GRCh37_hg19_ensembl_genes.bed", sep="\t", index=False, header=False)

# Get TSSs
tsss = annotation.apply(get_tss, axis=1)
tsss.to_csv(".data/external/GRCh37_hg19_ensembl_genes.tss.bed", sep="\t", index=False, header=False)
# sort bed file
tsss = pybedtools.BedTool("../data/GRCh37_hg19_ensembl_genes.tss.bed")
tsss.sort().saveas("../data/GRCh37_hg19_ensembl_genes.tss.bed")


# Get refseq bed file from UCSC, add 1 bp upstream, name as hg19.refSeq.TSS.bed
"sed 's/_up_1_.*//' hg19.refSeq.TSS.bed > t"
# Filter out ncRNAs
"grep NM t > hg19.refSeq.TSS.mRNA.bed"


# Get roadmap HMM state annotation
# read more about it here http://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state
roadmap_url = "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/"
roadmap_15statesHMM = "all.mnemonics.bedFiles.tgz"
os.system("wget {0}{1}".format(roadmap_url, roadmap_15statesHMM))

# Get all cell's HMM state annotation
os.system("gzip -d {0}".format(roadmap_15statesHMM))

# concatenate all files
all_states = "all_states_all_lines.bed"
os.system("cat *.bed > {0}".format(all_states))

# Get CD19 perypheral blood HMM state annotation
roadmap_15statesHMM_CD19 = "E032_15_coreMarks_mnemonics.bed.gz"
os.system("tar zxvf {0} {1}".format(roadmap_15statesHMM, roadmap_15statesHMM_CD19))
os.system("gzip -d {0}".format(roadmap_15statesHMM_CD19))
os.system("mv E032_15_coreMarks_mnemonics.bed ../data/E032_15_coreMarks_mnemonics.bed")


# Footprinting
# get all jaspar motifs
"wget http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt"
jaspar = motifs.parse(open("data/external/pfm_all.txt", 'r'), "jaspar")
# motif annotation
"wget http://jaspar.genereg.net/html/DOWNLOAD/database/MATRIX.txt"
annot = pd.read_table("data/external/MATRIX.txt", names=["index", "db", "id", 0, "TF"])
# get species annotation
"wget http://jaspar.genereg.net/html/DOWNLOAD/database/MATRIX_SPECIES.txt"
spec = pd.read_table("data/external/MATRIX_SPECIES.txt", names=["index", "species_id"])
# merge both
annot = pd.merge(annot, spec, on=['index'])

# get ids of only human motifs
human_annot = annot[annot['species_id'] == "9606"]

# filter out any not uniquely mappable gene name
human_annot = human_annot[
    (~human_annot['TF'].str.contains("\(")) &
    (~human_annot['TF'].str.contains(":")) &
    (~human_annot['TF'].str.contains("LM\d")) &  # filter out LM* motifs
    (human_annot['TF'] != "EWSR1-FLI1")
]

# get these from all jaspar motifs
human_motifs = [m for m in jaspar if m.base_id in human_annot['id'].tolist()]

# write back
with open("data/external/jaspar_human_motifs.txt", "w") as handle:
    handle.write(motifs.jaspar.write(human_motifs, format='jaspar'))

# write mapping of TF index and name
with open("data/external/jaspar_human_motifs.id_mapping.txt", "w") as handle:
    handle.write("\n".join(["\t".join([str(i), human_motifs[i - 1].base_id, human_motifs[i - 1].name.split(";")[1].upper()]) for i in range(1, 1 + len(human_motifs))]))


# TFBS-gene assignment:
# Get ensembl annotation (ensGenes) for mm10
names = [
    "name",  # Name of gene (usually transcript_id from GTF)
    "chrom",  # Chromosome name
    "strand",  # + or - for strand
    "txStart",  # Transcription start position
    "txEnd",  # Transcription end position
    "cdsStart",  # Coding region start
    "cdsEnd",  # Coding region end
    "exonCount",  # Number of exons
    "exonStarts",  # Exon start positions
    "exonEnds",  # Exon end positions
    "id",  # Unique identifier
    "name2",  # Alternate name (e.g. gene_id from GTF)
    "cdsStartStat",  # enum('none','unk','incmpl','cmpl')
    "cdsEndStat",  # enum('none','unk','incmpl','cmpl')
    "exonFra"  # Exon frame offsets {0,1,2}
]
genes = pd.read_table("data/external/ensGene.txt", names=names).reset_index(drop=True)
genes = genes[['chrom', 'txStart', 'txEnd', 'name2', 'name', 'strand']]
genes.columns = ['chrom', 'start', 'end', 'ensembl_gene', 'ensembl_transcript', 'strand']

# get longest transcript of each gene
indexes = genes.groupby("ensembl_gene").apply(lambda x: np.argmax(x['end'] - x['start']))
genes = genes.ix[indexes.unique().tolist()]

# make region with TSS
tss = genes.apply(get_tss, axis=1)
tss.to_csv("data/external/ensembl.tss.bed", sep="\t", index=False, header=False)

# make region with promoter + gene body
promoter_and_genesbody = genes.apply(get_promoter_and_genebody, axis=1)
promoter_and_genesbody.to_csv("data/external/ensembl.promoter_and_genesbody.bed", sep="\t", index=False, header=False)


# Get GOtermID - GOtermDescription mapping
# download data from Biomart (for some reason the GOterm name cannot be get automatically - probably I cannot figure out how :s)
# http://www.ensembl.org/biomart/martview/6451fcd5296302994382deee7bd9c8eb?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.name_1006|hsapiens_gene_ensembl.default.feature_page.go_id&FILTERS=&VISIBLEPANEL=resultspanel
"mv x data/external/goID_goName.csv"


# Puente 2015
# Get WGS variants
cmds = list()
cmds.append("wget http://www.nature.com/nature/journal/vaop/ncurrent/extref/nature14666-s3.zip")
cmds.append("unzip nature14666-s3.zip")
supFile = "data/external/Puente2015.supplementaryTable2.csv"
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


# GWAS studies
# get GWAS catalog
os.system("wget -O gwas_catalog.tsv http://www.ebi.ac.uk/gwas/api/search/downloads/alternative")  # gwas db dump
os.system("wget http://www.ebi.ac.uk/fgpt/gwas/ontology/GWAS-EFO-Mappings201405.xlsx")  # gwas mapping/ontology

# read in catalog and mappings
df = pd.read_csv("gwas_catalog.tsv", sep="\t")  # hg38!
mapping = pd.read_excel("GWAS-EFO-Mappings201405.xlsx")

# merge both
df2 = pd.merge(df, mapping, on=['PUBMEDID'])

# subset columns
df3 = df2[['CHR_ID', 'CHR_POS', 'PUBMEDID', 'DISEASE/TRAIT', 'PARENT', 'SNPS', 'STRONGEST SNP-RISK ALLELE', 'P-VALUE', 'OR or BETA']]
df3.columns = ['chr', 'pos', 'pubmed_id', 'trait', 'ontology_group', 'snp', 'snp_strongest_allele', 'p_value', 'beta']
df3.to_csv("data/external/gwas_catalog.csv", index=False)

# filter on p-value
df4 = df3[df3['p_value'] < 5e-8]

if not os.path.exists("regions"):
    os.makedirs("regions")

# export bed file for each ontology group
regionset_index = pd.DataFrame()
for group in df4['ontology_group'].unique():
    df5 = df4[df4['ontology_group'] == group]
    df5 = df5[['chr', 'pos']]
    df5.columns = ['chr', 'start']
    # drop entries without a position
    df5.dropna(how='any', subset=['chr', 'start'], inplace=True)
    df5['chr'] = ['chr' + str(int(i)) for i in df5['chr']]
    df5['end'] = df5['start'] + 1

    df5['start'] = df5['start'].astype(int)
    df5['end'] = df5['end'].astype(int)

    # write bed file
    group = re.sub(" ", "_", group).lower()
    df5.drop_duplicates().to_csv("regions/gwas_catalog.%s.bed" % group, sep="\t", header=False, index=False)

    # save in regionset index
    regionset_index = regionset_index.append(pd.Series(["Human", "SNPs in gwas catalog - %s" % group, "GWAS catalog", "gwas_catalog.%s.bed" % group]), ignore_index=True)

# save regionset index
regionset_index.columns = ["species", "description", "dataSource", "filename"]
regionset_index.to_csv("data/external/index.txt", sep="\t", index=False)


# CD19 DNase, Roadmap
accs = ["SRR171503", "SRR171504", "SRR171505", "SRR171506", "SRR171527", "SRR171528"]
os.system("cd %s " % os.path.join("data", "external"))
# get data
for acc in accs:
    os.system("fastq-dump %s " % acc)
os.system("cat %s > CD19_DNase.fastq" % " ".join([acc + ".fastq" for acc in accs]))

# map
cmd = bowtie2Map(
    inputFastq1="CD19_DNase.fastq",
    outputBam=os.path.join("data", "external", "CD19_DNase.bam"),
    log=os.path.join("data", "external", "CD19_DNase.log"),
    metrics=os.path.join("data", "external", "CD19_DNase.alnmetrics.txt"),
    genomeIndex="/data/groups/lab_bock/shared/resources/genomes/hg19/indexed_bowtie2/hg19",
    maxInsert=2000, cpus=24)
os.system(cmd)


# CpG islands
os.system("wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz")
os.system("cut -f 2,3,4 cpgIslandExt.txt > data/external/cpgIsland.hg19.bed")

#

# DNAme Kulis 2012

os.system("wget http://www.nature.com/ng/journal/v44/n11/extref/ng.2443-S2.xls")
df = pd.read_excel("ng.2443-S2.xls", skiprows=1)
df.rename(columns={"CHR": "chrom", "MAPINFO": "start"}, inplace=True)

# add chr
df["chrom"] = "chr" + df.chrom.astype(str)

# replace strs
df["Group"] = df["Group"].str.replace(" / ", "-").str.replace("Met ", "").str.replace("unmet ", "").str.replace("-CLL", "")

for group in df["Group"].unique():
    df2 = df[df["Group"] == group][["chrom", "start"]]
    df2["end"] = df2["start"] + 1
    df2[['chrom', 'start', 'end']].to_csv("data/external/Kulis_DNAme.%s.hg19.bed" % group, sep="\t", index=False, header=None)

# Save whole thing
df.to_csv("data/external/ng.2443-S2.csv", index=False)


# DNAme windows from Oakes 2016
os.system("wget http://www.nature.com/ng/journal/v48/n3/extref/ng.3488-S3.xlsx")
df = pd.read_excel("ng.3488-S3.xlsx", skiprows=1)

# get only windows with at least one TF motif
tfs = ["AP-1", "EBF1", "RUNX3", "IRF4", "OCT2", "NFkB"]
df2 = df[tfs]
df2 = df.ix[df2[df2.any(1)].index]

# Export Bed file
df2['chrom'] = df2["Window position (hg19)"].apply(lambda x: x.split(":")[0])
df2['start'] = df2["Window position (hg19)"].apply(lambda x: x.split(":")[1].split("-")[0])
df2['end'] = df2["Window position (hg19)"].apply(lambda x: x.split(":")[1].split("-")[1])

df2[['chrom', 'start', 'end']].to_csv("data/external/TF_DNAme_windows.hg19.bed", sep="\t", index=False, header=None)

# Split by TF
for TF in tfs:
    df2[df2[TF] == 1][['chrom', 'start', 'end']].to_csv("data/external/TF_DNAme_windows.%s.hg19.bed" % TF, sep="\t", index=False, header=None)

# Save whole thing
df2.to_csv("data/external/ng.3488-S3.csv", index=False)
