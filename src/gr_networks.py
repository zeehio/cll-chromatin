#!/usr/bin/env python

"""
This script builds patient- or patient group-specific
gene regulatory networks infered from transcription-factor footprints
in ATAC-seq data.
"""

import os
from looper.models import Project
from pypiper import ngstk as tk
import pandas as pd
import textwrap
import re
import pybedtools

# Set settings
pd.set_option("date_dayfirst", True)


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    """This returns joined patient and sample IDs"""
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_patient_id(name):
    return name.split("_")[2]


def name_to_sample_id(name):
    return name.split("_")[3]


def annotate_clinical_traits(samples):
    """
    Annotate samples with clinical traits.
    """
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
            setattr(sample, mut, 1 if sample.mutations is not pd.np.nan and mut in str(sample.mutations) else 0)

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
    """
    Annotate samples with all available traits.
    """
    new_samples = list()
    for sample in samples:
        # If any attribute is not set, set to NaN
        for attr in attrs:
            if not hasattr(sample, attr):
                setattr(sample, attr, pd.np.nan)
        new_samples.append(sample)

    # read in file with IGHV group of samples selected for ChIPmentation
    selected = pd.read_csv(os.path.join("metadata", "selected_samples.tsv"), sep="\t").astype(str)
    # annotate samples with the respective IGHV group
    for sample in samples:
        group = selected[
            (selected["patient_id"].astype(str) == str(sample.patient_id)) &
            (selected["sample_id"].astype(str) == str(sample.sample_id))
        ]["sample_cluster"]
        if len(group) == 1:
            sample.ighv_group = group.squeeze()
        else:
            sample.ighv_group = pd.np.nan

    return annotate_clinical_traits(annotate_disease_treatments(new_samples))


def piq_prepare_motifs(motif_file="data/external/jaspar_human_motifs.txt", n_motifs=366):
    """
    """
    cmds = list()
    for motif in range(1, n_motifs + 1):
        a = os.path.exists("/scratch/users/arendeiro/piq/motif.matches/%i.pwmout.RData" % motif)
        b = os.path.exists("/scratch/users/arendeiro/piq/motif.matches/%i.pwmout.rc.RData" % motif)
        if a and b:
            continue

        cmd = """
    Rscript ~/workspace/piq-single/pwmmatch.exact.r"""
        cmd += " ~/workspace/piq-single/common.r"
        cmd += " {0}".format(os.path.abspath(motif_file))
        cmd += " " + str(motif)
        cmd += """ /scratch/users/arendeiro/piq/motif.matches/
    """
        cmds.append(cmd)
    return cmds


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


def piq_footprint_background_single(bam_cache, background_cache, motif_number, tmp_dir, results_dir):
    """
    Footprint using PIQ.
    """
    cmd = """
    Rscript ~/workspace/piq-single/pertf.bg.r"""
    cmd += " ~/workspace/piq-single/common.r"
    cmd += " /scratch/users/arendeiro/piq/motif.matches/"
    cmd += " " + tmp_dir
    cmd += " " + results_dir
    cmd += " " + bam_cache
    cmd += " " + background_cache
    cmd += " " + str(motif_number)
    cmd += """
"""
    return cmd


def piq_footprint_differential_single(bam_cache1, bam_cache2, motif_number, tmp_dir, results_dir):
    """
    Footprint using PIQ.
    """
    cmd = """
    Rscript ~/workspace/piq-single/pertf.bg.r"""
    cmd += " ~/workspace/piq-single/common.r"
    cmd += " /scratch/users/arendeiro/piq/motif.matches/"
    cmd += " " + tmp_dir
    cmd += " " + results_dir
    cmd += " " + bam_cache1
    cmd += " " + bam_cache2
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


def prepare_background(bam_file):
    out_dir = os.path.abspath(os.path.join(data_dir, "footprinting_background"))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_cache = os.path.join(out_dir, "footprinting_background.RData")
    job_file = os.path.join(out_dir, "footprinting_background.slurm.sh")

    # Build job
    cmd = tk.slurmHeader(
        "footprinting_background", os.path.join(out_dir, "slurm.log"),
        cpusPerTask=24, queue="longq", time="7-12:00:00", memPerCpu=8000
    )

    # stupid PIQ hard-coded links
    cmd += """
    cd /home/arendeiro/workspace/piq-single/
    """
    # prepare cache
    cmd += piq_prepare_bams([bam_file], out_cache)

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
    promoter_and_genesbody = pybedtools.BedTool("data/external/ensembl.promoter_and_genesbody.bed")
    # read in TSS info
    tss = pybedtools.BedTool("data/external/ensembl.tss.bed")
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


def collect_networks(foots_dir, motif_numbers, label):
    gene_assignments, interactions, stats = piq_to_network(foots_dir, motif_numbers)

    # Write gene assignments to disk
    gene_assignments.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.gene_assignments.tsv"), sep="\t", index=False)
    # Write stats to disk
    stats.index = [number2tf[tf] for tf in stats.index]
    stats.to_csv(os.path.join(foots_dir, label + ".piq.TF-gene_interactions.stats.tsv"), sep="\t", index=True)

    # Drop TFBS with no gene in the same chromosome (random chroms) - very rare cases
    interactions = interactions[interactions['gene'] != '.']

    # Get original TF name and gene symbol
    interactions['TF'] = [number2tf[tf] for tf in interactions['TF']]
    interactions['gene'] = [ensembl2gene[gene] for gene in interactions['gene']]

    # Add interaction type
    interactions['interaction_type'] = "pd"

    # reorder columns
    interactions = interactions[["TF", "gene", "interaction_type", "interaction_score"]]

    # Save all TF-gene interactions
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
scratch_dir = os.path.join("/scratch/users/arendeiro/piq")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# Start project
prj = Project("metadata/project_config.yaml")
prj.add_sample_sheet()

# annotated samples with a few more things:
prj.samples = annotate_samples(prj.samples, prj.sheet.df.columns.tolist())
samples = prj.samples

# MOTIFS FOR PIQ
motifs_file = "data/external/jaspar_human_motifs.txt"
n_motifs = 366

# prepare motifs for footprinting (done once)
cmds = piq_prepare_motifs(motifs_file, n_motifs)
for cmd in cmds:
    cmd2 = tk.slurmHeader("PIQ_preparemotifs", os.path.join("/home/arendeiro/", "piq_preparemotifs.slurm.log"), cpusPerTask=1, queue="shortq")

    # stupid PIQ hard-coded links
    cmd2 += """
    cd /home/arendeiro/workspace/piq-single/
    """

    # add the piq command
    cmd2 += cmd

    # write job to file
    with open("/home/arendeiro/tmp.sh", 'w') as handle:
        handle.writelines(textwrap.dedent(cmd2))

    tk.slurmSubmitJob("/home/arendeiro/tmp.sh")


# MERGE BAM FILES FROM SAMPLES, PREPARE R CACHE FOR PIQ
# 'automatable' features
features = {
    # "patient_gender": ("F", "M"),  # gender
    "mutated": (True, False),  # ighv mutation
}

jobs = list()

# PREPARE R CACHE FROM BACKGROUND READS
jobs.append(prepare_background("/scratch/users/arendeiro/PGA1Nextera/PGA_0001_Nextera-2.bam"))

# all CLL samples
jobs.append(run_merged("all", [sample.filtered for sample in samples], "all"))

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

# individual samples
# for each sample create R cache with bam file
for sample in samples:
    if sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(data_dir, "footprints")
    r_data = os.path.abspath(os.path.join(foots_dir, sample.name + ".filteredshifted.RData"))
    if os.path.isfile(r_data):
        continue
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)

    tmp_dir = os.path.join(scratch_dir, sample.name)
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    job_file = os.path.join(foots_dir, "slurm_job.sh")

    # prepare slurm job header
    cmd = tk.slurmHeader(sample.name + "_PIQ_prepareBam", os.path.join(foots_dir, "slurm.log"), cpusPerTask=2, queue="shortq")

    # stupid PIQ hard-coded links
    cmd += """
    cd /home/arendeiro/workspace/piq-single/
    """

    # prepare bams
    cmd += piq_prepare_bams([sample.filteredshifted], r_data)

    # slurm footer
    cmd += tk.slurmFooter()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    # append file to jobs
    jobs.append(job_file)

# submit jobs (to create bam file caches)
for job in jobs:
    tk.slurmSubmitJob(job)


# FOOTPRINT
# get motifs
n_motifs = 366
motif_numbers = range(1, n_motifs + 1)
motifs_mapping = "data/external/jaspar_human_motifs.id_mapping.txt"
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

# individual samples
# for each sample launch several jobs (366) to footprint
for sample in samples:
    print(sample.name)
    if sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    foots_dir = os.path.join(sample.dirs.sampleRoot, "footprints")
    files = os.listdir(foots_dir)
    # os.system("rm %s/*" % os.path.join(foots_dir))
    if not os.path.exists(foots_dir):
        os.mkdir(foots_dir)
    r_data = os.path.abspath(os.path.join("data", "footprints", sample.name + ".filteredshifted.RData"))

    for motif in motif_numbers:
        # Check R cache for this motif exists
        if not os.path.exists("/scratch/users/arendeiro/piq/motif.matches/%i.pwmout.RData" % motif):
            print("Missing cache for motif %i" % motif)
            continue

        # if footprint files exist, skip motif:
        if len([files[i] for i in range(len(files)) if re.search("^%i-.*\.pdf$" % motif, files[i]) is not None]) == 2:
            continue

        t_dir = os.path.join(scratch_dir, sample.name)
        if not os.path.exists(t_dir):
            os.mkdir(t_dir)
        tmp_dir = os.path.join(scratch_dir, sample.name, str(motif))
        # if not os.path.exists(tmp_dir):
        #     os.mkdir(tmp_dir)
        job_file = os.path.join(foots_dir, "slurm_job_motif%i.sh" % motif)
        slurm_log = os.path.join(foots_dir, "slurm_motif%i.log" % motif)

        # prepare slurm job header
        cmd = tk.slurmHeader(sample.name + "_PIQ_footprinting_motif%i" % motif, slurm_log, cpusPerTask=2, queue="shortq", memPerCpu=8000)

        # stupid PIQ hard-coded links
        cmd += """
    cd /home/arendeiro/workspace/piq-single/
        """

        # footprint
        cmd += piq_footprint_single(r_data, motif, tmp_dir, results_dir=foots_dir)

        # slurm footer
        cmd += tk.slurmFooter()

        # write job to file
        with open(job_file, 'w') as handle:
            handle.writelines(textwrap.dedent(cmd))

        # append file to jobs
        jobs.append(job_file)

# submit jobs (to footprint)
try:
    for job in jobs:
        import time
        from subprocess import Popen, PIPE

        process = Popen(["squeue"], stdout=PIPE)
        (output, err) = process.communicate()

        if len(output.split("\n")) - 2 > 2000:
            time.sleep(60)
            if len(output.split("\n")) - 2 > 2000:
                continue

        tk.slurmSubmitJob(job)
except KeyboardInterrupt("Interrupted"), e:
    raise e


# BUILD NETWORKS
ensembl2gene = pd.read_table("data/GRCh37_hg19_ensembl_genes.bed", sep="\t", header=None)
ensembl2gene = dict(zip(ensembl2gene[3], ensembl2gene[6]))

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

# individual samples
# connect each motif to a gene
for sample in samples:
    if sample.technique != "ATAC-seq" or sample.cellLine != "CLL":
        continue

    if os.path.exists(os.path.join(data_dir, "footprints", sample.name + ".piq.TF-gene_interactions.tsv")):
        continue

    print(sample.name)
    collect_networks(os.path.join(sample.dirs.sampleRoot, "footprints"), motif_numbers, sample.name)


# CIS-REGULATORY MODULE USAGE BETWEEN CLL TYPES
# Get CLL groups
# Select genome region of relevance
# Get reads for regions at e.g. 1kb resolution
# Correlate genome positions in each group

# Get differential between groups


# Examples:
# Open chromatin defined by DNaseI and FAIREidentifies regulatory elements that shape cell-type identity
# Lingyun Song et al.
# Figure 6A

# Patterns of regulatory activity across diverse human cell types predict tissue identity, transcription factor binding, and long-range interactions
# Nathan C. Sheffield et al.
