
import os
from pipelines.models import Project
from pipelines import toolkit as tk
import pandas as pd
import numpy as np
import textwrap
import re
import pybedtools


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
        cpusPerTask=12, queue="longq", time="7-12:00:00", memPerCpu=8000
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


def piq_to_network(results_dir, motif_numbers):
    """
    Parse PIQ output, filter footprints.
    Returns matrix with likelyhood score of each TF regulating each gene.
    """
    # list results_dir
    files = os.listdir(results_dir)
    # get all cll peaks to filter data
    all_peaks = pybedtools.BedTool("data/cll_peaks.bed")
    # read in gene info
    refseq_mrna_tss = pybedtools.BedTool("data/hg19.refSeq.TSS.mRNA.deduplicated.bed")

    # dict to store TF->gene interactions
    interactions = pd.DataFrame()

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
            else:
                df["end"] = df["start"]
                df["start"] = df["start"] - 1
            # concatenate
            if i == 0:
                df2 = df
            else:
                df2 = pd.concat([df, df2])

        # Filter for purity
        footprints = df2[df2["purity"] > 0.7]

        # If empty give 0 to every gene for this TF
        if len(footprints) < 500:
            continue

        footprints[['chr', 'start', 'end', 'pwm', 'shape', 'score', 'purity']].to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # filter for motifs overlapping CLL peaks
        footprints = pybedtools.BedTool(os.path.join("tmp.bed")).intersect(all_peaks, wa=True).to_dataframe()
        footprints.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity"]
        footprints.to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)

        # Get closest gene
        closest_tss = pybedtools.BedTool(os.path.join("tmp.bed")).closest(refseq_mrna_tss, d=True).to_dataframe()
        closest_tss.columns = ["chrom", "start", "end", "pwm", "shape", "score", "purity", "chrom_gene", "start_gene", "end_gene", "gene", "distance"]

        # Get weighted values
        # weigh with footprint purity and distance to tss
        scores = closest_tss
        scores['interaction_score'] = scores.apply(lambda x: 2 * (x['purity'] - 0.5) * 10 ** -(x['distance'] / 1000000.), axis=1)
        # sum scores for each gene
        scores = scores.groupby(['gene'])['interaction_score'].apply(sum).reset_index()

        scores["TF"] = motif

        interactions = pd.concat([interactions, scores])

    return interactions


def collect_networks(foots_dir, motif_numbers, label):
    interactions = piq_to_network(foots_dir, motif_numbers)

    # Drop TFBS with no gene in the same chromosome (random chroms) - very rare cases
    interactions = interactions[interactions['gene'] != '.']

    # Get original TF name and gene symbol
    interactions['TF'] = [number2tf[tf] for tf in interactions['TF']]
    interactions['gene'] = [refseq2gene[gene] for gene in interactions['gene']]

    interactions['interaction_type'] = "pd"
    interactions.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.tsv"), sep="\t", index=False)

    # Filter for TF-> TF interactions
    interactions_TF = interactions[interactions['gene'].isin(tfs)]
    interactions_TF.to_csv(os.path.join(foots_dir, label + ".piq.TF-TF_interactions.tsv"), sep="\t", index=False)

    # Filter for nodes with more than 2 edges
    interactions_TF_filtered = interactions_TF[interactions_TF['interaction_score'] >= 1]
    interactions_TF_filtered.to_csv(os.path.join(foots_dir, label + ".piq.TF-TF_interactions.filtered.tsv"), sep="\t", index=False)


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

jobs = list()

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
motifs_file = "~/workspace/piq-single/pwms/jasparfix.txt"
n_motifs = 1316

# read list of tfs to do
df = pd.read_csv("data/tf_gene_matching.txt", sep="\t", header=None)
df[1] = [x.upper() for x in df[1]]
tfs = df[1]
number2tf = dict(zip(df[0], df[1]))
motif_numbers = df[0]

# send out jobs
jobs = list()
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
refseq2gene = pd.read_csv("data/Refseq2Gene.txt", sep="\t", header=None)
refseq2gene = dict(zip(refseq2gene[0], refseq2gene[1]))

# parse PIQ output,
# connect each motif to a gene

# "gender" and "mutated"
for i, (feature, (group1, group2)) in enumerate(features.items()):
    # append file to jobs
    foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group1)]), "footprints"))
    collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", feature, str(group1)]))
    foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", feature, str(group2)]), "footprints"))
    collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", feature, str(group2)]))

# untreated vs 1st line chemotherapy +~ B cell antibodies
foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", "untreated_vs_1stline", "untreated"]), "footprints"))
collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", "untreated_vs_1stline", "1stlinetreatment"]))

# Disease at Diagnosis - comparison in untreated samples
# CLL vs MBL
foots_dir = os.path.abspath(os.path.join(data_dir, "_".join(["merged-samples", "CLL_vs_MBL", "CLL"]), "footprints"))
collect_networks(foots_dir, motif_numbers, "_".join(["merged-samples", "CLL_vs_MBL", "MBL"]))
