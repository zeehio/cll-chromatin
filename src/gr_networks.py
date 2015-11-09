#!/usr/bin/env python

"""
This script builds patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import os
from pipelines.models import Project
from pipelines import toolkit as tk
import pandas as pd
import numpy as np
import textwrap
import re
import pybedtools
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt
import seaborn as sns


# Set settings
pd.set_option("date_dayfirst", True)
sns.set_style("whitegrid")
sns.set_context("paper")


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    """This returns joined patient and sample IDs"""
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_patient_id(name):
    return name.split("_")[2]


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


def piq_prepare_bams(bams, output_cache):
    """
    :type bams: list
    :type output_cache: str
    """
    # add before:
    # salloc -c 12 -p develop --mem-per-cpu 8000
    cmd = """
    Rscript ~/workspace/piq-single/bam2rdata.r ~/workspace/piq-single/common.r {0} """.format(output_cache)
    cmd += " ".join(bams)
    cmd += """
    """
    return cmd


def piq_footprint_single(bam_cache, motif_number, tmp_dir, results_dir):
    """
    Footprint using PIQ.
    """
    cmd = """
    Rscript ~/workspace/piq-single/pertf.r"""
    cmd += " ~/workspace/piq-single/common.r"
    cmd += " /scratch/users/arendeiro/piq/motif.matches/"
    cmd += " " + tmp_dir
    cmd += " " + results_dir
    cmd += " " + bam_cache
    cmd += " " + str(motif_number)
    cmd += """
"""
    return cmd


def run_merged(feature, bam_files, group_label):
    merge_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, group_label])))
    if not os.path.exists(merge_dir):
        os.mkdir(merge_dir)
    foots_dir = os.path.join(merge_dir, "footprints")
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)

    merged_bam = os.path.join(merge_dir, "merged.bam")
    merged_cache = os.path.join(merge_dir, "merged.RData")
    job_file = os.path.join(merge_dir, "slurm.sh")

    # Build job
    cmd = tk.slurmHeader(
        "_".join(["CLL_merged-samples", feature, group_label]), os.path.join(merge_dir, "slurm.log"),
        cpusPerTask=12, queue="longq", time="7-12:00:00", memPerCpu=4000
    )

    # merge all bam files
    cmd += """
    sambamba merge -t 12 {0} {1}
    """.format(merged_bam, " ".join(bam_files))

    # stupid PIQ hard-coded links
    cmd += """
    cd /home/arendeiro/workspace/piq-single/
    """
    # prepare cache
    cmd += piq_prepare_bams([merged_bam], merged_cache)

    # slurm footer
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    # append file to jobs
    return job_file


def footprint(feature, group_label, motif_numbers):
    scratch_dir = os.path.join("/scratch/users/arendeiro/piq")
    merge_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, group_label])))
    foots_dir = os.path.join(merge_dir, "footprints")
    merged_cache = os.path.join(merge_dir, "merged.RData")
    label = "_".join(["merged-samples", feature, group_label])

    jobs = list()
    for motif in motif_numbers:
        if not os.path.exists("/scratch/users/arendeiro/piq/motif.matches/%i.pwmout.RData" % motif):
            print("piq file for motif %i does not exist" % motif)
            continue

        t_dir = os.path.join(scratch_dir, label)
        if not os.path.exists(t_dir):
            os.mkdir(t_dir)
        tmp_dir = os.path.join(scratch_dir, label, str(motif))
        # if not os.path.exists(tmp_dir):
        #     os.mkdir(tmp_dir)
        job_file = os.path.join(foots_dir, "slurm_job_motif%i.sh" % motif)
        slurm_log = os.path.join(foots_dir, "slurm_motif%i.log" % motif)

        # prepare slurm job header
        cmd = tk.slurmHeader(label + "_PIQ_footprinting_motif%i" % motif, slurm_log, cpusPerTask=2, queue="shortq", memPerCpu=8000)

        # stupid PIQ hard-coded links
        cmd += """
    cd /home/arendeiro/workspace/piq-single/
        """

        # footprint
        cmd += piq_footprint_single(merged_cache, motif, tmp_dir, results_dir=foots_dir)

        # slurm footer
        cmd += tk.slurmFooter()

        # write job to file
        with open(job_file, 'w') as handle:
            handle.writelines(textwrap.dedent(cmd))

        # append file to jobs
        jobs.append(job_file)
    return jobs


def tfbs_to_gene(bed_file):
    # read in gene body + promoter info
    promoter_and_genesbody = pybedtools.BedTool("data/ensembl.promoter_and_genesbody.bed")
    # read in TSS info
    tss = pybedtools.BedTool("data/ensembl.tss.bed")
    # columns
    columns = ["chrom", "start", "end", "pwm", "shape", "strand", "score", "purity",
               "chrom_gene", "start_gene", "end_gene", "gene", "transcript", "strand_gene"]

    # Assign TFBS to gene if they overlap with gene body or promoter (5kb around TSS -> 2.5kb upstream)
    gene_assignments = pybedtools.BedTool(os.path.join(bed_file)).intersect(promoter_and_genesbody, wa=True, wb=True).to_dataframe().reset_index()
    gene_assignments.columns = columns

    # For the remaining TFBSs, assign TFBS to closest TSS regardless of distance
    # (distance is not so important for assignment because distance is a penalyzing effect during TF-gene interaction score calculation)
    # 1. get genes not assigned previously
    all_ = pybedtools.BedTool(os.path.join(bed_file)).to_dataframe().reset_index()
    merged = pd.merge(all_, gene_assignments, how="left", on=['chrom', 'start', 'end'])
    remaining = merged[merged['gene'].isnull()]
    remaining.icol(range(1, 9)).to_csv(os.path.join("tmp_rest.bed"), sep="\t", index=False, header=False)

    # 2. assign to nearest
    closest_tss = pybedtools.BedTool(os.path.join("tmp_rest.bed")).closest(tss, d=True).to_dataframe().reset_index()
    closest_tss.columns = columns + ['distance']

    # put the two together
    gene_assignments = pd.concat([gene_assignments, closest_tss])

    # set overlapping distance to 0
    gene_assignments.loc[gene_assignments['distance'].isnull(), 'distance'] = 0

    return gene_assignments


def piq_to_network(results_dir, motif_numbers):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """
    # list results_dir
    files = os.listdir(results_dir)
    # get all cll peaks to filter data
    all_peaks = pybedtools.BedTool("data/cll_peaks.bed")

    # dataframe to store TFBS assignment to genes
    assignments = pd.DataFrame()

    # dataframe to store TF->gene interactions
    interactions = pd.DataFrame()

    # dataframe to store stats about the TFBS and the interactions
    stats = pd.DataFrame()

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in motif_numbers:
        # get both forward and reverse complement PIQ output files
        result_files = list()
        for f in files:
            m = re.match(r'%i-.*-calls\.csv$' % motif, f)
            if hasattr(m, "string"):
                result_files.append(m.string)

        # make bed file from it
        # concatenate files (forward and reverse complement are treated differently by PIQ)
        for i, result_file in enumerate(result_files):
            df = pd.read_csv(os.path.join(results_dir, result_file), index_col=0)
            df.rename(columns={"coord": "start"}, inplace=True)
            # fix coordinates
            if "RC-calls.csv" not in result_file:
                df["end"] = df["start"] + 1
                df['strand'] = "+"
            else:
                df["end"] = df["start"]
                df["start"] = df["start"] - 1
                df['strand'] = "-"
            # concatenate
            if i == 0:
                df2 = df
            else:
                df2 = pd.concat([df, df2])

        # add total TFBS to stats
        stats.loc[motif, "TFBS"] = len(df2)
        stats.loc[motif, "TFBS_+"] = len(df2[df2['strand'] == "+"])
        stats.loc[motif, "TFBS_-"] = len(df2[df2['strand'] == "-"])

        # Filter for purity
        footprints = df2[df2["purity"] > 0.7]
        stats.loc[motif, "pur0.7"] = len(footprints)

        # If less than 500 significant interactions, ignore TF
        if len(footprints) < 500:
            continue

        footprints[['chr', 'start', 'end', 'pwm', 'shape', 'strand', 'score', 'purity']].to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # filter for motifs overlapping CLL peaks
        footprints = pybedtools.BedTool(os.path.join("tmp.bed")).intersect(all_peaks, wa=True).to_dataframe()
        footprints.columns = ['chr', 'start', 'end', 'pwm', 'shape', 'strand', 'score', 'purity']
        footprints.to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)
        stats.loc[motif, "overlap_cll"] = len(footprints)

        # assign TFBS to gene
        gene_assignments = tfbs_to_gene(os.path.join("tmp.bed"))
        stats.loc[motif, "gene_overlap_count"] = len(gene_assignments[gene_assignments['distance'] == 0])
        stats.loc[motif, "dist_gene_median"] = gene_assignments['distance'].median()
        stats.loc[motif, "dist_gene_std"] = gene_assignments['distance'].std()

        # Get weighted values
        # weigh with footprint purity and distance to tss
        gene_assignments['interaction_score'] = gene_assignments.apply(lambda x: 2 * (x['purity'] - 0.5) * 10 ** -(x['distance'] / 1000000.), axis=1)
        # sum scores for each gene
        scores = gene_assignments.groupby(['gene'])['interaction_score'].apply(sum).reset_index()
        scores["TF"] = motif

        # add mean score for each gene
        stats.loc[motif, "score_gene_mean"] = scores['interaction_score'].mean()
        stats.loc[motif, "score_gene_std"] = scores['interaction_score'].std()

        # add to dataframe with all TF-gene interactions
        interactions = pd.concat([interactions, scores])
        assignments = pd.concat([assignments, gene_assignments])

    return (assignments, interactions, stats)


def piq_to_change(results_dir1, results_dir2, motif_numbers):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """

    def get_purity(motif, results_dir, files):
        # get all cll peaks to filter data
        all_peaks = pybedtools.BedTool("data/cll_peaks.bed")
        # get both forward and reverse complement PIQ output files
        result_files = list()
        for f in files:
            m = re.match(r'%i-.*-calls\.all\.csv$' % motif, f)
            if hasattr(m, "string"):
                result_files.append(m.string)

        if len(result_files) != 2:
            return None
        # make bed file from it
        # concatenate files (forward and reverse complement are treated differently by PIQ)
        for i, result_file in enumerate(result_files):
            df = pd.read_csv(os.path.join(results_dir, result_file), index_col=0)
            df.rename(columns={"coord": "start"}, inplace=True)
            # fix coordinates
            if "RC-calls.csv" not in result_file:
                df["end"] = df["start"] + 1
            else:
                df["end"] = df["start"]
                df["start"] = df["start"] - 1
            # concatenate
            if i == 0:
                df2 = df
            else:
                df2 = pd.concat([df, df2])

        # Filter for purity
        df2[['chr', 'start', 'end', 'pwm', 'shape', 'score', 'purity']].to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # filter for motifs overlapping CLL peaks
        footprints = pybedtools.BedTool(os.path.join("tmp.bed")).intersect(all_peaks, wa=True).to_dataframe()
        footprints.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity"]

        return footprints["purity"]

    # list results_dir
    files1 = os.listdir(results_dir1)
    files2 = os.listdir(results_dir2)

    # dict to store TF->gene interactions
    fold_changes = list()

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in motif_numbers:
        print motif

        pur1 = get_purity(motif, results_dir1, files1)
        pur2 = get_purity(motif, results_dir2, files2)

        if pur1 is None or pur2 is None:
            continue

        pur1 = pur1.apply(lambda x: x if x > 0.7 else 0)
        pur2 = pur2.apply(lambda x: x if x > 0.7 else 0)

        fold_changes.append((motif, pur1.mean(), pur2.mean(), mannwhitneyu(pur1, pur2)[1]))

    return fold_changes


def collect_networks(foots_dir, motif_numbers, label):
    gene_assignments, interactions, stats = piq_to_network(foots_dir, motif_numbers)

    # Write gene assignments to disk
    gene_assignments.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.gene_assignments.tsv"), sep="\t", index=False)
    # Write stats to disk
    stats.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.stats.tsv"), sep="\t", index=False)

    # Drop TFBS with no gene in the same chromosome (random chroms) - very rare cases
    interactions = interactions[interactions['gene'] != '.']

    # Get original TF name and gene symbol
    interactions['TF'] = [number2tf[tf] for tf in interactions['TF']]
    interactions['gene'] = [ensembl2gene[gene] for gene in interactions['transcript']]

    # Save all TF-gene interactions
    interactions['interaction_type'] = "pd"
    interactions.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.tsv"), sep="\t", index=False)

    # Filter for TF-gene interactions stronger than 1
    interactions_filtered = interactions[interactions['interaction_score'] >= 1]
    interactions_filtered.to_csv(os.path.join(foots_dir, label + ".piq.TF-TF_interactions.filtered.tsv"), sep="\t", index=False)

    # Filter for TF-> TF interactions
    interactions_TF = interactions[interactions['gene'].isin(number2tf.values())]
    interactions_TF.to_csv(os.path.join(foots_dir, label + ".piq.TF-TF_interactions.tsv"), sep="\t", index=False)

    # Filter for TF-TF interactions stronger than 1
    interactions_TF_filtered = interactions_TF[interactions_TF['interaction_score'] >= 1]
    interactions_TF_filtered.to_csv(os.path.join(foots_dir, label + ".piq.TF-TF_interactions.filtered.tsv"), sep="\t", index=False)


# INIT
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

# Annotate with clinical data
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
    'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19']  # 'CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19']
samples = [sample for sample in prj.samples if sample.cellLine == "CLL" and sample.technique == "ATAC-seq" and sample.name not in samples_to_exclude]


# MERGE BAM FILES FROM SAMPLES, PREPARE R CACHE FOR PIQ
jobs = list()

# all CLL samples
jobs.append(run_merged("all", [sample.filtered for sample in samples], "all"))


# "gender" and "mutated"
features = {
    "patient_gender": ("F", "M"),  # gender
    "mutated": (True, False),  # ighv mutation
}
for i, (feature, (group1, group2)) in enumerate(features.items()):
    # example : i, (feature, (group1, group2)) = (0, (features.items()[0]))
    # get dataframe subset with groups
    g1 = [sample.filtered for sample in samples if getattr(sample, feature) == group1]
    g2 = [sample.filtered for sample in samples if getattr(sample, feature) == group2]

    # append file to jobs
    jobs.append(run_merged(feature, g1, str(group1)))
    jobs.append(run_merged(feature, g2, str(group2)))


# untreated vs 1st line chemotherapy +~ B cell antibodies
g1 = [sample.filtered for sample in samples if not sample.treatment_active and not sample.relapse]
drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR', 'Alemtuz']
g2 = [sample.filtered for sample in samples if sample.treatment_active and sample.treatment_type in drugs]
jobs.append(run_merged("untreated_vs_1stline", g1, "untreated"))
jobs.append(run_merged("untreated_vs_1stline", g2, "1stlinetreatment"))

# Disease at Diagnosis - comparison in untreated samples
# CLL vs MBL
g1 = [sample.filtered for sample in samples if sample.diagnosis_disease == "CLL" and not sample.treatment_active and not sample.relapse]
g2 = [sample.filtered for sample in samples if sample.diagnosis_disease == "MBL" and not sample.treatment_active and not sample.relapse]
jobs.append(run_merged("CLL_vs_MBL", g1, "CLL"))
jobs.append(run_merged("CLL_vs_MBL", g2, "MBL"))

for job in jobs:
    tk.slurmSubmitJob(job)


# FOOTPRINT
# get motifs
n_motifs = 366
motif_numbers = range(1, n_motifs + 1)
motifs_mapping = "data/jaspar_human_motifs.id_mapping.txt"
df = pd.read_table(motifs_mapping, header=None)
number2tf = dict(zip(df[0], df[2]))


# send out jobs
jobs = list()

# all samples
jobs += footprint("all", "all", motif_numbers)

# "gender" and "mutated"
for i, (feature, (group1, group2)) in enumerate(features.items()):
    # append file to jobs
    jobs += footprint(feature, str(group1), motif_numbers)
    jobs += footprint(feature, str(group2), motif_numbers)

# untreated vs 1st line chemotherapy +~ B cell antibodies
jobs += footprint("untreated_vs_1stline", "untreated", motif_numbers)
jobs += footprint("untreated_vs_1stline", "1stlinetreatment", motif_numbers)

# Disease at Diagnosis - comparison in untreated samples
# CLL vs MBL
jobs += footprint("CLL_vs_MBL", "CLL", motif_numbers)
jobs += footprint("CLL_vs_MBL", "MBL", motif_numbers)

for job in jobs:
    tk.slurmSubmitJob(job)


# BUILD NETWORKS
ensembl2gene = pd.read_csv("data/ensemblToGeneName.txt", sep="\t", header=None)
ensembl2gene = dict(zip(ensembl2gene[0], ensembl2gene[1]))

# parse PIQ output,
# connect each motif to a gene

# all samples
foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", "all", "all"]), "footprints"))
collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", "all", "all"]))


# "gender" and "mutated"
for i, (feature, (group1, group2)) in enumerate(features.items()):
    # i, (feature, (group1, group2)) = (0, (features.items()[0]))
    print feature
    # append file to jobs
    foots_dir1 = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group1)]), "footprints"))
    collect_networks(foots_dir1, motif_numbers, "_".join(["merged-samples", feature, str(group1)]))
    foots_dir2 = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group2)]), "footprints"))
    collect_networks(foots_dir2, motif_numbers, "_".join(["merged-samples", feature, str(group2)]))

# untreated vs 1st line chemotherapy +~ B cell antibodies
foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", "untreated_vs_1stline", "untreated"]), "footprints"))
collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", "untreated_vs_1stline", "1stlinetreatment"]))

# Disease at Diagnosis - comparison in untreated samples
# CLL vs MBL
foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", "CLL_vs_MBL", "CLL"]), "footprints"))
collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", "CLL_vs_MBL", "MBL"]))


# FOLD-CHANGE PLOT
# u/mCLL
i, (feature, (group1, group2)) = (0, (features.items()[0]))

foots_dir1 = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group1)]), "footprints"))
foots_dir2 = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group2)]), "footprints"))

df = pd.DataFrame(piq_to_change(foots_dir1, foots_dir2, motif_numbers))
df.columns = ['TF', 'mean1', 'mean2', 'pvalues']
df['TF'] = [number2tf[tf] for tf in df['TF']]

# sns.jointplot(df['mean1'], df['mean2'], xlim=(0.5, 1), ylim=(0.5, 1))

# # zscore transform
# df['mean1z'] = (df['mean1'] - df['mean1'].mean()) / df['mean1'].std(ddof=0)
# df['mean2z'] = (df['mean2'] - df['mean2'].mean()) / df['mean2'].std(ddof=0)
# sns.jointplot(df['mean1z'], df['mean2z'])

# # volcano plot
# df = pd.DataFrame(results)
# df.columns = ['TF', 'mean1', 'mean2', 'pvalues']
# df['TF'] = [number2tf[tf] for tf in df['TF']]
# sns.jointplot(df['fold'], -np.log10(1 + df['pvalues']))

# scatter
df['fold'] = df['mean1'] - df['mean2']
df['pvalues'] = multipletests(df.pvalues, method="fdr_bh", alpha=0.05)[1]

for i in range(len(df)):
    if df.ix[i]['fold'] < 0:
        col = "#DEDB43"
    elif df.ix[i]['fold'] > 0:
        col = "#31CC7E"
    else:
        col = "gray"

    plt.plot(
        df.ix[i]['mean1'], df.ix[i]['mean2'], 'o',
        alpha=1.0 if df.ix[i]['pvalues'] < 1e-80 else 0.2,
        color=col,
        label=df.ix[i]['TF'] if df.ix[i]['pvalues'] < 1e-80 else None
    )
    if df.ix[i]['pvalues'] < 1e-80:
        plt.annotate(
            df.ix[i]['TF'],
            xy=(df.ix[i]['mean1'], df.ix[i]['mean2']), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.3),
            arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0')
        )
plt.xlim((0.5, 1))
plt.ylim((0.5, 1))
plt.plot([0, 1], [0, 1], '--', color='black')

plt.savefig(os.path.join(plots_dir, "_".join(["merged-samples", feature, str(group1), str(group2), "fold_change.svg"])))
