#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import pybedtools
import re
from Bio import motifs


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
genes = pd.read_csv("data/ensGene.txt", sep="\t", header=None)
genes = genes[[2, 4, 5, 12, 1, 3]]
genes.columns = ['chrom', 'start', 'end', 'gene', 'transcript', 'strand']

# Annotate with gene names
names = pd.read_csv("data/ensemblToGeneName.txt", sep="\t", header=None)
names.columns = ['transcript', 'name']
annotation = pd.merge(genes, names)
annotation.to_csv("data/GRCh37_hg19_ensembl_genes.bed", sep="\t", index=False, header=False)

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
jaspar = motifs.parse(open("data/pfm_all.txt", 'r'), "jaspar")
# motif annotation
"wget http://jaspar.genereg.net/html/DOWNLOAD/database/MATRIX.txt"
annot = pd.read_table("data/MATRIX.txt", names=["index", "db", "id", 0, "TF"])
# get species annotation
"wget http://jaspar.genereg.net/html/DOWNLOAD/database/MATRIX_SPECIES.txt"
spec = pd.read_table("data/MATRIX_SPECIES.txt", names=["index", "species_id"])
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
with open("data/jaspar_human_motifs.txt", "w") as handle:
    handle.write(motifs.jaspar.write(human_motifs, format='jaspar'))

# write mapping of TF index and name
with open("data/jaspar_human_motifs.id_mapping.txt", "w") as handle:
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
genes = pd.read_table("data/ensGene.txt", names=names).reset_index(drop=True)
genes = genes[['chrom', 'txStart', 'txEnd', 'name2', 'name', 'strand']]
genes.columns = ['chrom', 'start', 'end', 'ensembl_gene', 'ensembl_transcript', 'strand']

# get longest transcript of each gene
indexes = genes.groupby("ensembl_gene").apply(lambda x: np.argmax(x['end'] - x['start']))
genes = genes.ix[indexes.unique().tolist()]

# make region with TSS
tss = genes.apply(get_tss, axis=1)
tss.to_csv("data/ensembl.tss.bed", sep="\t", index=False, header=False)

# make region with promoter + gene body
promoter_and_genesbody = genes.apply(get_promoter_and_genebody, axis=1)
promoter_and_genesbody.to_csv("data/ensembl.promoter_and_genesbody.bed", sep="\t", index=False, header=False)


# Get GOtermID - GOtermDescription mapping
# download data from Biomart (for some reason the GOterm name cannot be get automatically - probably I cannot figure out how :s)
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
df3.to_csv("gwas_catalog.csv", index=False)

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
regionset_index.to_csv("index.txt", sep="\t", index=False)
