#!/usr/bin/env python

"""
Variant-calling from ATAC-seq data.
"""

import os
from pipelines.models import Project, ATACseqSample

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]

cmds = list()
cmds.append("module load samtools/0.1.19")
cmds.append("samtools mpileup -uf ~/resources/genomes/hg19/hg19.fa {0} | bcftools view -bvcg - > data/variants.raw.bcf".format(" ".join([s.mapped for s in samples])))
cmds.append("bcftools view data/variants.raw.bcf | vcfutils.pl varFilter -D100 > data/variants.flt.vcf")

for cmd in cmds:
    os.system(cmd)
