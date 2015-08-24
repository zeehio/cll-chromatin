#!/usr/bin/env python

import yaml
import os
from pipelines import Project, ATACseqSample


# Read configuration file
with open("config.yaml", 'r') as handle:
    config = yaml.load(handle)
dataDir = os.path.join(config["paths"]["parent"], config["projectname"], "data")

# Start project
prj = Project("cll-patients")
prj.addSampleSheet("../metadata/sample_annotation.csv")


# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample]


# COMPARE CLL WITH NORMAL CELLS
# we'll have matched B-cells from the patients as well


# DE NOVO/CLL-SPECIFIC ENHANCERS
# Find unique enhancers across CLL samples compared with normal B-cells
# Search other cell types for overlaping enhancers:
# - if possitive -> enhancer activation
# - if negative -> de novo enhancer -> explore mechanism
# validate with H3K27ac ChIP-seq
# validate with RNA expression
