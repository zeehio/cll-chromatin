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


# FOOTPRINTING
# ideally with PIQ (seems to be the best performing method)

# get factor occupancy score (FOS) for each factor in each site
# cluster again patients based on TF occupancy


# NETWORKS
# Network construction:
# - after footprinting
# - assign a confident TF footprint to its target gene:
#     - is it in a promoter? is it in a known enhancer of a gene?
# - build network of A -> B based on this

# Network types:
# patient-specific:
# - build network specific to each patient
# - compare networks
# - cluster patients based on network connections

# for groups of patients:
# - come up with a way of combining signal from several patients from one group
# - build networks specific to groups
# - compare networks
