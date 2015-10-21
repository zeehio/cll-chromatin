#!/usr/bin/env python

"""
This script merges bam files from groups of samples
and generated bigWig files from these.
"""

import os
import sys
import numpy as np
import pandas as pd
from pipelines.models import Project
from pipelines import toolkit as tk
import textwrap


def name_to_sample_id(name):
    return name.split("_")[3]


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


def annotate_treatments(samples, clinical):
    """
    Annotate samples with timepoint, treatment_status, treatment_type
    """
    def string_to_date(string):
        if type(string) is str:
            if len(string) == 10:
                return pd.to_datetime(string, format="%d/%m/%Y")
            if len(string) == 7:
                return pd.to_datetime(string, format="%m/%Y")
            if len(string) == 4:
                return pd.to_datetime(string, format="%Y")
        return pd.NaT

    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL":
            # get sample id
            sample_id = name_to_sample_id(sample.name)

            # get corresponding series from "clinical"
            sample_c = clinical[clinical['sample_id'] == sample_id].squeeze()

            # Get sample collection date
            sample.collection_date = string_to_date(sample_c['sample_collection_date'])
            # Get diagnosis date
            sample.diagnosis_date = string_to_date(sample_c['diagnosis_date'])
            # Get diagnosis disease
            sample.diagnosis_disease = sample_c['diagnosis_disease'] if type(sample_c['diagnosis_disease']) is str else None
            # Get time since diagnosis
            sample.time_since_diagnosis = sample.collection_date - sample.diagnosis_date

            # Get all treatment dates
            treatment_dates = [string_to_date(date) for date in sample_c[["treatment_%i_date" % (x) for x in range(1, 5)]].squeeze()]
            # Get treatment end date
            treatment_end_date = string_to_date(clinical[(clinical['sample_id'] == sample_id)][["treatment_end_date"]])
            # Check if there are earlier "timepoints"
            earlier_dates = [treatment_date for treatment_date in treatment_dates if treatment_date < sample.collection_date]

            # Annotate samples with active treatment
            for treatment_date in treatment_dates:
                # if one of the treatment dates is earlier
                if treatment_date < sample.collection_date:
                    # this sample was not collected at diagnosis time
                    sample.diagnosis_collection = False
                    # and no treatment end date in between, mark as under treatment
                    if treatment_end_date is pd.NaT:
                        sample.treatment_active = True
                    else:
                        if treatment_date < treatment_end_date < sample.collection_date:
                            sample.treatment_active = False
                        elif treatment_date < sample.collection_date < treatment_end_date:
                            sample.treatment_active = True
            # if there were no treatments before collection, consider untreated
            if not hasattr(sample, "treatment_active"):
                sample.treatment_active = False
                # if there were no treatments before collection, and collection was within 30 days of diagnosis, tag as collected at diagnosis
                if sample.time_since_diagnosis is not pd.NaT:
                    if abs(sample.time_since_diagnosis) < pd.to_timedelta(30, unit="days"):
                        sample.diagnosis_collection = True
            if not hasattr(sample, "diagnosis_collection"):
                sample.diagnosis_collection = False

            # Annotate treatment type, time since treatment
            if sample.treatment_active:
                if len(earlier_dates) > 0:
                    # Find out which earlier "timepoint" is closest and annotate treatment and response
                    previous_dates = [date for date in clinical[(clinical['sample_id'] == sample_id)][["treatment_%i_date" % (x) for x in range(1, 5)]].squeeze()]
                    closest_date = previous_dates[np.argmin([abs(date - sample.collection_date) for date in earlier_dates])]

                    # Annotate previous treatment date
                    sample.previous_treatment_date = string_to_date(closest_date)
                    # Annotate time since treatment
                    sample.time_since_treatment = sample.collection_date - string_to_date(closest_date)

                    # Get closest clinical "timepoint", annotate response
                    closest_timepoint = [tp for tp in range(1, 5) if closest_date == sample_c["treatment_%i_date" % tp]][0]

                    sample.treatment_type = sample_c['treatment_%i_regimen' % closest_timepoint]
                    sample.treatment_response = sample_c['treatment_%i_response' % closest_timepoint]

            # Annotate relapses
            # are there previous timepoints with good response?
            # Get previous clinical "timepoints", annotate response
            if len(earlier_dates) > 0:
                closest_timepoint = [tp for tp in range(1, 5) if closest_date == sample_c["treatment_%i_date" % tp]][0]

                # Annotate with previous known response
                sample.previous_response = sample_c['treatment_%i_response' % closest_timepoint]

                # if prior had bad response, mark current as relapse
                if sample_c['treatment_%i_response' % closest_timepoint] in ["CR", "GR"]:
                    sample.relapse = True
                else:
                    sample.relapse = False
            else:
                sample.relapse = False

        # If any attribute is not set, set to None
        for attr in ['diagnosis_collection', 'diagnosis_date', 'diagnosis_disease', 'time_since_treatment', 'treatment_type',
                     'treatment_response', "treatment_active", "previous_treatment_date", "previous_response", 'relapse']:
            if not hasattr(sample, attr):
                setattr(sample, attr, None)

        # Append sample
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


def annotate_mutations(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            sample.mutations = clinical[clinical['sample_id'] == _id]['mutations'].tolist()[0]
        else:
            sample.mutations = None
        new_samples.append(sample)
    return new_samples


def merge_bams(bams, output_bam):
    """
    Decorator for some methods of Analysis class.
    """
    job_file = "/scratch/users/arendeiro/tmp.sh"
    cmd = tk.slurmHeader("merge_bams", os.path.join("/scratch/users/arendeiro/", "merge_bams.slurm.log"), cpusPerTask=8, time='6-10:00:00', queue="longq", memPerCpu=8000)

    cmd += """
    samtools merge {0} {1}
    """.format(output_bam, " ".join(bams))
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


def bamToBigWig(inputBam, outputBigWig, tagmented=False, normalize=False):
    import os
    import re

    genomeSizes = "/data/groups/lab_bock/shared/resources/genomes/hg19/hg19.chromSizes"
    genome = "hg19"

    cmd = tk.slurmHeader("bam_to_bigwig", os.path.join("/scratch/users/arendeiro/", "merge_bams.slurm.log"), cpusPerTask=8, time='6-10:00:00', queue="longq", memPerCpu=8000)

    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))

    cmd1 = """
    bedtools bamtobed -i {0} |""".format(inputBam)
    if not tagmented:
        cmd1 += " bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
    cmd1 += " genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
            "-5 " if tagmented else "",
            genomeSizes,
            transientFile
    )
    cmd += cmd1

    if normalize:
        cmd += """
    awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov
    """.format(transientFile)

    cmd += """
    bedGraphToBigWig {0}{1}.cov {2} {3}
    """.format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig)

    # remove tmp files
    cmd += """
    if [[ -s {0}.cov ]]; then rm {0}.cov; fi
    """.format(transientFile)
    if normalize:
        cmd += """
    if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi
    """.format(transientFile)

    cmd += """
    chmod 755 {0}
    """.format(outputBigWig)

    cmd += tk.slurmFooter()

    # write job to file
    job_file = "/scratch/users/arendeiro/tmp.sh"
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


def main():
    # Start project
    prj = Project("cll-patients")
    prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

    # Annotate with clinical data
    clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))
    prj.samples = annotate_igvh_mutations(prj.samples, clinical)
    prj.samples = annotate_treatments(prj.samples, clinical)
    prj.samples = annotate_mutations(prj.samples, clinical)
    prj.samples = annotate_gender(prj.samples, clinical)

    # Start analysis object
    # only with ATAC-seq samples that passed QC
    samples_to_exclude = [
        'CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19',
        'CLL_ATAC-seq_5186_1-5-57350_ATAC17-4_hg19',
        'CLL_ATAC-seq_4784_1-5-52817_ATAC17-6_hg19',
        'CLL_ATAC-seq_981_1-5-42480_ATAC16-6_hg19',
        'CLL_ATAC-seq_5277_1-5-57269_ATAC17-8_hg19',
        'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19',
        'CLL_ATAC-seq_5147_1-5-48105_ATAC17-2_hg19',
        'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19']
    # Use these samples only
    samples = [sample for sample in prj.samples if sample.cellLine == "CLL" and sample.name not in samples_to_exclude and sample.technique == "ATAC-seq"]

    # ALL CLL SAMPLES
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "merged-samples", "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "merged-samples", "merged.bigwig"))
    merge_bams([sample.filtered for sample in samples], merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    # TRAIT-SPECIFIC
    # "gender" and "mutated"
    features = {
        "patient_gender": ("F", "M"),  # gender
        "mutated": (True, False),  # ighv mutation
    }
    for i, (feature, (group1, group2)) in enumerate(features.items()):
        # example : i, (feature, (group1, group2)) = (0, (features.items()[0]))
        # get dataframe subset with groups
        g1 = [[sample.filtered for sample in samples if getattr(sample, feature) == group1]]
        merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", feature, str(group1)]), "merged.bam"))
        merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", feature, str(group1) + ".bigwig"])))
        if not os.path.exists(merged_bam):
            merge_bams(g1, merged_bam)
        if not os.path.exists(merged_bigwig):
            bamToBigWig(merged_bam, merged_bigwig)

        g2 = [[sample.filtered for sample in samples if getattr(sample, feature) == group2]]
        merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", feature, str(group2)]), "merged.bam"))
        merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", feature, str(group2) + ".bigwig"])))
        if not os.path.exists(merged_bam):
            merge_bams(g2, merged_bam)
        if not os.path.exists(merged_bigwig):
            bamToBigWig(merged_bam, merged_bigwig)

    # untreated vs treated
    g1 = [sample.filtered for sample in samples if not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "untreated_vs_treated", "untreated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "untreated_vs_treated", "untreated" + ".bigwig"])))
    merge_bams(g1, merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    g2 = [sample.filtered for sample in samples if sample.treatment_active]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "untreated_vs_treated", "treated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "untreated_vs_treated", "treated" + ".bigwig"])))
    merge_bams(g2, merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    # untreated vs 1st line chemotherapy +~ B cell antibodies
    g1 = [sample.filtered for sample in samples if not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "untreated_vs_1stline", "untreated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "untreated_vs_1stline", "untreated" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g1, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR', 'Alemtuz']
    g2 = [sample.filtered for sample in samples if sample.treatment_active and sample.treatment_type in drugs]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "untreated_vs_1stline", "1stline"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "untreated_vs_1stline", "1stline" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g2, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    # Disease at Diagnosis - comparison in untreated samples
    # CLL vs MBL
    g1 = [sample.filtered for sample in samples if sample.diagnosis_disease == "CLL" and not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "CLL_vs_MBL", "CLL"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "CLL_vs_MBL", "CLL" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g1, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    g2 = [sample.filtered for sample in samples if sample.diagnosis_disease == "MBL" and not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.dirs.data, "_".join(["merged-samples", "CLL_vs_MBL", "MBL"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.dirs.html, "_".join(["merged-samples", "CLL_vs_MBL", "MBL" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g2, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
