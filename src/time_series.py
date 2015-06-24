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


# CHROMATIN OPENNESS AS A FUNCTION OF DISEASE PROGRESSION AND TREATMENT
# Cluster sites based on pattern of oppenness across time-points
# test enrichment of site clusters

# Think how to analyse data from few patients tretated with several drugs at different times
