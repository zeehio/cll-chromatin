#!/usr/bin/env python

"""
This script calls supernenhncers for the cll-patients project.
"""

import os
from pipelines.models import Project
import pipelines.toolkit as tk


def macs2CallPeaks(treatmentBam, outputDir, sampleName, genome, controlBam=None, broad=False):
    """
    Use MACS2 to call peaks.
    """

    sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}

    if not broad:
        cmd = "macs2 callpeak -t {0}".format(treatmentBam)
        if controlBam is not None:
            cmd += " -c {0}".format(controlBam)
        cmd += " --bw 200 -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)
        # --fix-bimodal --extsize 180
    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        cmd = "macs2 callpeak -t {0}".format(treatmentBam)
        if controlBam is not None:
            cmd += " -c {0}".format(controlBam)
        cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3 -g {0} -n {1} --outdir {2}".format(
            sizes[genome], sampleName, outputDir
        )

    return cmd


def rose_superenhancers(
        input_bams, ranking_bam, input_bed,
        source_code_dir, output_dir, control_bam=None, genome="HG19"):
    """
    Use ROSE to call superenhancers.
    """
    import os

    def narrowPeak_to_gff(bed_file, gff_file):
        """
        Convert BED6 format to Gff suitable for ROSE input.
        """
        import pandas as pd
        try:
            bed = pd.read_csv(bed_file, sep="\t")
            bed.columns = ["chrom", "start", "end", "name", "score", "strand"] + ["_"] + ["."] * 3
            bed["name"] = ["p_%i" % i for i in range(bed.shape[0])]
            gff = bed[["chrom", "name", "_", "start", "end", "_", "strand", "_", "name"]]
            gff.to_csv(gff_file, sep="\t", index=False, header=False)
        except ValueError:
            print("Input file is empty")
            os.system("touch %s" % gff_file)

    # produce gff file from bed file
    if ".gff" in input_bed:
        print("Assuming input is in gff format.")
        gff_file = input_bed
    else:
        print("Assuming input is in narrowPeak format.")
        gff_file = os.path.join(os.path.dirname(input_bed), os.path.basename(input_bed.split(".")[0]) + ".gff")
        narrowPeak_to_gff(input_bed, gff_file)
        print(gff_file)

    # support multiple bam files as input
    if type(input_bams) is str:
        input_bams = [input_bams]

    # make output folder
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # assemble ROSE command
    cmd = """
    cd {}
    """.format(source_code_dir)
    cmd += """
    python ROSE_main.py \\
    -i {} \\
    -r {} \\
    -o {} \\
    -g {} \\
    -b {}""".format(
        gff_file,
        ranking_bam,
        output_dir,
        genome.upper(),
        ",".join(input_bams))

    # add control bam
    if control_bam is not None:
        cmd += """\\
    -c {}""".format(control_bam)

    # cd back to previous directory
    cmd += """

    cd -
    """

    return cmd


def call_peaks(sample, igg_sample):
    """
    Call H3K27ac peaks for pairs of sample and igg control.
    """
    import textwrap

    output_dir = sample.paths.peaks
    job_file = os.path.join(output_dir, "slurm.sh")
    # prepare slurm job header
    cmd = tk.slurmHeader(
        sample.name + "_enhancer_calling", os.path.join(output_dir, "slurm.log"),
        cpusPerTask=8,
        queue="shortq")

    # footprint
    cmd += macs2CallPeaks(
        sample.filtered,
        output_dir,
        sample.name,
        sample.genome,
        controlBam=igg_sample.filtered,
        broad=False)

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


def call_superenhancers(sample, igg_sample):
    """
    Call superenhancers for pairs of sample and igg control.
    """
    import textwrap

    output_dir = os.path.join(sample.paths.sample_root, "superenhancers")
    job_file = os.path.join(output_dir, "slurm.sh")
    # prepare slurm job header
    cmd = tk.slurmHeader(
        sample.name + "_superenhancer_calling", os.path.join(output_dir, "slurm.log"),
        cpusPerTask=8,
        queue="longq")

    # footprint
    cmd += rose_superenhancers(
        sample.filtered,
        sample.filtered,
        sample.peaks,
        "/home/arendeiro/workspace/rose/",
        output_dir,
        control_bam=igg_sample.filtered,
        genome=sample.genome.upper()
    )

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurmSubmitJob(job_file)


# Start project
prj = Project("metadata/project_config.yaml")
prj.add_sample_sheet()

# Start analysis object
# only with ATAC-seq samples that passed QC
samples_to_exclude = ["CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19", "CLL_ATAC-seq_AKH13_M3152_ATAC40s21_hg19"]
samples = [sample for sample in prj.samples if sample.cell_line == "CLL" and sample.name not in samples_to_exclude and sample.library == "ChIPmentation"]

# pair ChIP and IgG
for _id in set([s.sample_id for s in samples]):
    ip = [s for s in samples if s.sample_id == _id and s.ip == "H3K27ac"][0]
    igg = [s for s in samples if s.sample_id == _id and s.ip == "IgG"][0]

    print(_id, ip.name, igg.name)

    # Call peaks for all H3K27ac ChIPmentation samples
    call_peaks(ip, igg)

    # Call superenhancers for all ChIPmentation samples
    call_superenhancers(ip, igg)
