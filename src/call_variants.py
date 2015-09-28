#!/usr/bin/env python

"""
Variant-calling from ATAC-seq data.
"""

import os
from pipelines.models import Project

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if s.technique == "ATAC-seq" and s.cellLine == "CLL" and s.readType == "PE"]

cmds = list()
cmds.append("module load samtools/0.1.19")
cmds.append("samtools mpileup -uf ~/resources/genomes/hg19/hg19.fa {0} | bcftools view -bvcg - > data/variants.raw.bcf".format(" ".join([s.filtered for s in samples])))
cmds.append("bcftools view data/variants.raw.bcf | vcfutils.pl varFilter -D100 > data/variants.flt.vcf")

for cmd in cmds:
    os.system(cmd)


# Make bed file, intersect with CLL peaks
# cut header, save
# grep "#" data/variants.flt.vcf > data/variants.flt.header
# cut header out, duplicate position column -> bed file
# grep -v "#" data/variants.flt.vcf | awk -v OFS='\t' '{$2=$2 "\t" $2}{print}' > data/variants.flt.bed

# intersect with peaks
# bedtools intersect -wa -a data/variants.flt.bed -b data/cll-peaks.bed > data/variants.flt.in-peaks.bed
