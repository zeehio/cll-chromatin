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
# - Go to each TF motif site annotated: http://compbio.mit.edu/encode-motifs/
# - get window around (e.g. +/- 100bp), footprint with CENTIPEDE or PIQ
# - save posterior probabilities for every site

# ideally with PIQ (seems to be the best performing method)

# more:
# get factor occupancy score (FOS) for each factor in each site
# cluster again patients based on TF occupancy


# NETWORKS
# Network construction:

# Calculate a "regulation score" between EVERY TF and EVERY gene:
# - for each TF, for each chromossome, sum the weighted posterior probabilities of factor A regulating gene I.
# weigh these by dividing by the (log) of the distance of each binding site to the TSS of gene I.
# This will give a matrix of TS-gene "regulation scores".
# Investigate the distribution of scores.

# Compare regulation across patients:
# - for each TF, correlate the scores with the scores of all other transcription factors in the other patient.

# Additionally
# - build network of TF -> gene based on:
# 	- binding at the promoter or at known enhancers for the gene
#   OR
#   - consider regulation only if score is above threshold

# Network types:
# patient-specific:
# - build network specific to each patient
# - compare networks
# - cluster patients based on network connections

# for groups of patients:
# - come up with a way of combining signal from several patients from one group
# - build networks specific to groups
# - compare networks
