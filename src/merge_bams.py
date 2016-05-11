#!/usr/bin/env python

"""
This script merges bam files from groups of samples
and generated bigWig files from these.
"""

import os
import sys
import pandas as pd
from pipelines.models import Project
from pipelines import toolkit as tk
import textwrap


def name_to_sample_id(name):
    return name.split("_")[3]


def annotate_clinical_traits(samples):
    # Annotate traits
    chemo_drugs = ["Chlor", "Chlor R", "B Of", "BR", "CHOPR"]  # Chemotherapy
    target_drugs = ["Alemtuz", "Ibrutinib"]  # targeted treatments
    muts = ["del13", "del11", "tri12", "del17"]  # chrom abnorms
    muts += ["SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "TP53", "MYD88", "CHD2", "NFKIE"]  # mutations
    for s in samples:
        # Gender
        s.gender = 1 if s.patient_gender == "M" else 0 if s.patient_gender == "F" else pd.np.nan
        # IGHV mutation status
        s.IGHV = s.ighv_mutation_status

    # Annotate samples which are under treament but with different types
    for sample in samples:
        if not sample.under_treatment:
            sample.chemo_treated = pd.np.nan
            sample.target_treated = pd.np.nan
        else:
            sample.chemo_treated = 1 if sample.treatment_regimen in chemo_drugs else 0
            sample.target_treated = 1 if sample.treatment_regimen in target_drugs else 0
        for mut in muts:
            setattr(sample, mut, 1 if sample.mutations is not pd.np.nan and mut in sample.mutations else 0)

    return samples


def annotate_disease_treatments(samples):
    """
    Annotate samples with timepoint, treatment_status, treatment_regimen
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
        if sample.cell_line == "CLL":
            # Get sample collection date
            sample.collection_date = string_to_date(sample.sample_collection_date)
            # Get diagnosis date
            sample.diagnosis_date = string_to_date(sample.diagnosis_date)
            # Get diagnosis disease
            sample.primary_CLL = 1 if sample.diagnosis_disease == "CLL" else 0  # binary label useful for later

            # Get time since diagnosis
            sample.time_since_diagnosis = sample.collection_date - sample.diagnosis_date

            # Annotate treatment type, time since treatment
            if sample.under_treatment:
                sample.time_since_treatment = sample.collection_date - string_to_date(sample.treatment_date)

        # Append sample
        new_samples.append(sample)
    return new_samples


def annotate_samples(samples, attrs):
    new_samples = list()
    for sample in samples:
        # If any attribute is not set, set to NaN
        for attr in attrs:
            if not hasattr(sample, attr):
                setattr(sample, attr, pd.np.nan)
        new_samples.append(sample)

    # read in file with IGHV group of samples selected for ChIPmentation
    selected = pd.read_csv(os.path.join("metadata", "selected_samples.tsv"), sep="\t")  # .astype(str)
    # annotate samples with the respective IGHV group
    for sample in samples:
        group = selected[
            (selected["patient_id"] == sample.patient_id) &
            (selected["sample_id"] == sample.sample_id)
        ]["sample_cluster"]
        if len(group) == 1:
            sample.ighv_group = group.squeeze()
        else:
            sample.ighv_group = pd.np.nan

    return annotate_clinical_traits(annotate_disease_treatments(new_samples))


def merge_bams(bams, output_bam):
    """
    Decorator for some methods of Analysis class.
    """
    job_file = "/scratch/users/arendeiro/tmp.sh"
    cmd = tk.slurmHeader("merge_bams", os.path.join("/scratch/users/arendeiro/", "merge_bams.slurm.log"), cpusPerTask=8, time='10:00:00', queue="develop", memPerCpu=8000)

    cmd += """
    samtools merge {0} {1}
    """.format(output_bam, " ".join(bams))
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


def sort_index_bam(bam):
    """
    Decorator for some methods of Analysis class.
    """
    job_file = "/scratch/users/arendeiro/tmp.sh"
    cmd = tk.slurmHeader("sort_merge_bam", os.path.join("/scratch/users/arendeiro/", "sort_merge_bam.slurm.log"), cpusPerTask=8, time='10:00:00', queue="shortq", memPerCpu=8000)

    cmd += """
    sambamba sort -t 8 {}

    """.format(bam)
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


def bamToBigWig(inputBam, outputBigWig, tagmented=False, normalize=True):
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


def bamToBigWig_RNA(inputBam, outputBigWig, tagmented=False, normalize=True):
    import os
    import re

    genomeSizes = "/data/groups/lab_bock/shared/resources/genomes/hg19/hg19.chromSizes"

    cmd = tk.slurmHeader("bam_to_bigwig", os.path.join("/scratch/users/arendeiro/", "merge_bams.slurm.log"), cpusPerTask=8, time='6-10:00:00', queue="longq", memPerCpu=8000)

    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))

    cmd1 = """
    bedtools bamtobed -i {0} |""".format(inputBam)
    cmd1 += " genomeCoverageBed -bg -g {1} -i stdin > {2}.cov".format(genomeSizes, transientFile)
    cmd += cmd1

    if normalize:
        cmd += """
    awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov
    """.format(transientFile)

    cmd += """
    bedGraphToBigWig {0}{1}.cov {2} {3}
    """.format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig)

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
    prj = Project("metadata/project_config.yaml")
    prj.add_sample_sheet()

    # annotated samples with a few more things:
    prj.samples = annotate_samples(prj.samples, prj.sheet.df.columns.tolist())

    samples = [
        sample for sample in prj.samples if
        (sample.library == "ATAC-seq")]

    # ALL CLL SAMPLES
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "merged-samples", "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "merged-samples", "merged.bigwig"))
    merge_bams([sample.filtered for sample in samples], merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    # TRAIT-SPECIFIC
    # "gender" and "mutated"
    features = {
        "patient_gender": ("F", "M"),  # gender
        "ighv_mutation_status": (True, False),  # ighv mutation
    }
    for i, (feature, (group1, group2)) in enumerate(features.items()):
        # example : i, (feature, (group1, group2)) = (0, (features.items()[0]))
        # get dataframe subset with groups
        g1 = [sample.filtered for sample in samples if getattr(sample, feature) == group1]
        merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", feature, str(group1)]), "merged.sorted.bam"))
        merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", feature, str(group1) + ".bigwig"])))
        if not os.path.exists(merged_bam):
            merge_bams(g1, merged_bam)
        if not os.path.exists(merged_bigwig):
            bamToBigWig(merged_bam, merged_bigwig)

        g2 = [sample.filtered for sample in samples if getattr(sample, feature) == group2]
        merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", feature, str(group2)]), "merged.sorted.bam"))
        merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", feature, str(group2) + ".bigwig"])))
        if not os.path.exists(merged_bam):
            merge_bams(g2, merged_bam)
        if not os.path.exists(merged_bigwig):
            bamToBigWig(merged_bam, merged_bigwig)

    # untreated vs treated
    g1 = [sample.filtered for sample in samples if not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "untreated_vs_treated", "untreated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "untreated_vs_treated", "untreated" + ".bigwig"])))
    merge_bams(g1, merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    g2 = [sample.filtered for sample in samples if sample.treatment_active]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "untreated_vs_treated", "treated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "untreated_vs_treated", "treated" + ".bigwig"])))
    merge_bams(g2, merged_bam)
    bamToBigWig(merged_bam, merged_bigwig)

    # untreated vs 1st line chemotherapy +~ B cell antibodies
    g1 = [sample.filtered for sample in samples if not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "untreated_vs_1stline", "untreated"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "untreated_vs_1stline", "untreated" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g1, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR', 'Alemtuz']
    g2 = [sample.filtered for sample in samples if sample.treatment_active and sample.treatment_type in drugs]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "untreated_vs_1stline", "1stline"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "untreated_vs_1stline", "1stline" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g2, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    # Disease at Diagnosis - comparison in untreated samples
    # CLL vs MBL
    g1 = [sample.filtered for sample in samples if sample.diagnosis_disease == "CLL" and not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "CLL_vs_MBL", "CLL"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "CLL_vs_MBL", "CLL" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g1, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    g2 = [sample.filtered for sample in samples if sample.diagnosis_disease == "MBL" and not sample.treatment_active and not sample.relapse]
    merged_bam = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged-samples", "CLL_vs_MBL", "MBL"]), "merged.bam"))
    merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged-samples", "CLL_vs_MBL", "MBL" + ".bigwig"])))
    if not os.path.exists(merged_bam):
        merge_bams(g2, merged_bam)
    if not os.path.exists(merged_bigwig):
        bamToBigWig(merged_bam, merged_bigwig)

    # Merge ChIPmentation samples of same group
    samples = [
        sample for sample in prj.samples if
        (sample.library == "ChIPmentation")]

    for group in ["uCLL", "iCLL", "mCLL"]:
        out_dir = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged_chipmentation_samples", group])))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for ip in ["H3K27ac", "H3K4me1", "H3K27me3", "IgG"]:
            # get samples from group
            bams = [sample.filtered for sample in samples if sample.ighv_group == group and sample.ip == ip]
            merged_bam = os.path.join(out_dir, "merged_{}.bam".format(ip))
            merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged_chipmentation_samples", group + "_{}.bigwig".format(ip)])))

            # Launch jobs
            if not os.path.exists(merged_bam):
                merge_bams(bams, merged_bam)
            print merged_bigwig
            if not os.path.exists(merged_bigwig):
                bamToBigWig(merged_bam, merged_bigwig)

    # Merge RNA-seq samples of same group
    samples = [
        sample for sample in prj.samples if
        (sample.library == "RNA-seq")]

    for group in ["uCLL", "iCLL", "mCLL"]:
        out_dir = os.path.abspath(os.path.join(prj.paths.results_subdir, "_".join(["merged_RNA-seq_samples", group])))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # get samples from group
        bams = [os.path.join(sample.paths.sample_root, sample.name + ".genome_aligned.bam") for sample in samples if sample.ighv_group == group]
        if not all([os.path.exists(bam) for bam in bams]):
            print "{} don't exist!".format(" ".join(bams))
            continue

        merged_bam = os.path.join(out_dir, "merged_RNA-seq_{}.bam".format(group))
        merged_bigwig = os.path.abspath(os.path.join(prj.config['trackhubs']['trackhub_dir'], "_".join(["merged_RNA-seq_samples", group + "_{}.bigwig".format(group)])))

        # Launch jobs
        # merge_bams(bams, merged_bam)
        # sort_index_bam(merged_bam)

        merged_bam = os.path.join(out_dir, "merged_RNA-seq_{}.sorted.bam".format(group))
        if not os.path.exists(merged_bigwig):
            bamToBigWig_RNA(merged_bam, merged_bigwig)

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
