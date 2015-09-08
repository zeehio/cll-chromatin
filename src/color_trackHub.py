#!/usr/bin/env python

"""
Add info to track hub.
Color code according to sample features.
"""

import os
from pipelines.models import Project
import pandas as pd
import re


def name_to_sample_id(name):
    return name.split("_")[3:4][0]


def annotate_igvh_mutations(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            if clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 1:
                sample.mutated = True
            elif clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 2:
                sample.mutated = False
            else:
                sample.mutated = None
        else:
            sample.mutated = None
        new_samples.append(sample)
    return new_samples


def annotate_gender(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            if clinical.loc[clinical['sample_id'] == _id, 'patient_gender'].tolist()[0] == "F":
                sample.patient_gender = "F"
            elif clinical.loc[clinical['sample_id'] == _id, 'patient_gender'].tolist()[0] == "M":
                sample.patient_gender = "M"
            else:
                sample.patient_gender = None
        else:
            sample.patient_gender = None
        new_samples.append(sample)
    return new_samples


# Get path configuration
data_dir = os.path.join('.', "data")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# Get clinical info
clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))

# Start project
# prj = pickle.load(open("prj.pickle", 'rb'))
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Annotate with igvh mutated status
# add "mutated" attribute to sample depending on IGVH mutation status
# Select ATAC-seq samples
prj.samples = annotate_igvh_mutations(prj.samples, clinical)
prj.samples = annotate_gender(prj.samples, clinical)


# read tracks
with open(os.path.join(prj.dirs.html, "trackHub_hg19.txt"), 'r') as handle:
    lines = handle.readlines()

mut = list()
unmut = list()
rest = list()

for i, line in enumerate(lines):
    if i == 0:
        continue
    fields = re.search("'(.*)'", line).group().split("\'")[1].strip().split("_")
    if fields[2] == "ChIP-seq" and len(fields) > 6:
        patiend_id = fields[4]
        sample_id = fields[5]
    elif fields[0] != "CLL":
        rest.append(line + "\n")
        continue
    else:
        patiend_id = fields[2]
        sample_id = fields[3]

    sample = [s for s in prj.samples if s.sampleID == sample_id][0]

    # add gender and mutation status to title
    if sample.mutated:
        mutated = "Y"
    elif sample.mutated is None:
        mutated = "?"
    elif not sample.mutated:
        mutated = "N"
    line = re.sub("_hg19 ", "; gender={0}; igvh mutated={1}".format(str(sample.patient_gender), mutated), line)

    # change color according to IGVH
    if sample.mutated:
        mut.append(re.sub("color=.*\n", "color=0,127,102", line) + "\n")
    if not sample.mutated:
        unmut.append(re.sub("color=.*\n", "color=196,43,89", line) + "\n")
    elif sample.mutated is None:
        rest.append(re.sub("color=.*\n", "color=196,43,89", line) + "\n")

# write back
with open(os.path.join(prj.dirs.html, "trackHub_hg19.color_coded.txt"), 'w') as handle:
    handle.write(lines[0])
    handle.writelines(unmut)
    handle.writelines(mut)
    handle.writelines(rest)

# make executable
os.chmod(os.path.join(prj.dirs.html, "trackHub_hg19.color_coded.txt"), 0655)
