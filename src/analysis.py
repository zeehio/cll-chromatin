#!/usr/bin/env python

"""
This is the main script of the cll-patients project.
"""

# %logstart  # log ipython session

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['backend'] = "Agg"
import recipy
from argparse import ArgumentParser
import os
import sys
from pipelines.models import Project
import pybedtools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
import multiprocessing
import parmap
import pysam
import numpy as np
import pandas as pd
# import dask.dataframe as dd
try:  # stupid error, importing it twice works
    from sklearn import cross_validation
except AttributeError:
    from sklearn import cross_validation
from sklearn.preprocessing import normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from scipy.cluster.hierarchy import fcluster
import cPickle as pickle
from collections import Counter


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white")
sns.set_style({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica']})
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex=False)


def pickle_me(function):
    """
    Decorator for some methods of Analysis class.
    """
    def wrapper(obj, *args):
        function(obj, *args)
        pickle.dump(obj, open(obj.pickle_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    return wrapper


class Analysis(object):
    """
    Class to hold functions and data from analysis.
    """

    def __init__(self, data_dir, plots_dir, samples, pickle_file):
        self.data_dir = data_dir
        self.plots_dir = plots_dir
        self.samples = samples
        self.pickle_file = pickle_file

    @pickle_me
    def to_pickle(self):
        pass

    def from_pickle(self):
        return pickle.load(open(self.pickle_file, 'rb'))

    @pickle_me
    def get_consensus_sites(self, samples):
        """Get consensus (union) sites across samples"""

        for i, sample in enumerate(samples):
            print(sample.name)
            # Get summits of peaks and window around them
            peaks = pybedtools.BedTool(sample.peaks)  # .slop(g=pybedtools.chromsizes('hg19'), b=250)
            # Merge overlaping peaks within a sample
            peaks = peaks.merge()
            if i == 0:
                sites = peaks
            else:
                # Concatenate all peaks
                sites = sites.cat(peaks)

        # Merge overlaping peaks across samples
        sites = sites.merge()

        # Filter
        # remove blacklist regions
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "wgEncodeDacMapabilityConsensusExcludable.bed"))
        # remove chrM peaks and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(self.data_dir, "cll_peaks.bed"))

        # Read up again
        self.sites = pybedtools.BedTool(os.path.join(self.data_dir, "cll_peaks.bed"))

    def estimate_peak_saturation(self, samples, n=1000):
        """
        Randomizes sample order n times and measures the number of unique peaks after adding each sample after the other.
        Plots this.
        """
        import random
        from scipy import stats

        peak_count = np.zeros((len(samples), n))
        for i in range(n):
            samples_copy = samples[:]
            random.shuffle(samples_copy)
            for j, sample in enumerate(samples_copy):
                print(sample.name)
                # Get summits of peaks and window around them
                peaks = pybedtools.BedTool(sample.peaks)  # .slop(g=pybedtools.chromsizes('hg19'), b=250)
                # Merge overlaping peaks within a sample
                peaks = peaks.merge()
                if j == 0:
                    sites = peaks
                else:
                    # Concatenate all peaks
                    sites = sites.cat(peaks)
                # Let's keep track of the number of new peaks found with each new sample
                peak_count[j, i] = len(sites)

        # Make data frame
        df = pd.DataFrame(peak_count)

        # Calculate confidence intervals
        ci = df.apply(lambda x: pd.Series(stats.norm.interval(0.95, loc=x.mean(), scale=x.std() / np.sqrt(len(x)))), axis=1)
        ci.columns = ['low', 'high']

        # Add sample number
        ci['samples'] = np.array(range(len(samples))) + 1

        # Plot
        fig, axis = plt.subplots(1)
        axis.plot(ci['samples'], ci['low'])
        axis.plot(ci['samples'], ci['high'])
        axis.set_xlabel("# samples")
        axis.set_ylabel("# peaks")
        axis.set_ylim((0, 120000))
        plt.legend(loc='best')
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.cum_peak_count.95CI.svg"), bbox_inches="tight")

        # Store
        df.to_csv(os.path.join(self.data_dir, "cll_peaks.cum_peak_count.95CI.csv"), index=False)

    def estimate_reads_sequenced(self, n=100):
        """
        Randomizes sample order n times and measures the number of unique peaks after adding each sample after the other.
        Plots this.
        """
        import random
        from scipy import stats

        # GET CONSENSUS SITES ACROSS CLL ATAC-SEQ SAMPLES
        samples = [sample for sample in self.samples if sample.cell_line == "CLL" and sample.library == "ATAC-seq"]

        read_count = np.zeros((len(samples), n))
        for i in range(n):
            samples_copy = samples[:]
            random.shuffle(samples_copy)
            for j, sample in enumerate(samples_copy):
                print(sample.name)
                # Count number of filtered read pairs
                read_count[j, i] = read_count[j - 1, i] + int(pysam.Samfile(sample.filtered).mapped)

        # Make data frame
        df = pd.DataFrame(read_count)

        # Calculate confidence intervals
        ci = df.apply(lambda x: pd.Series(stats.norm.interval(0.95, loc=x.mean(), scale=x.std() / np.sqrt(len(x)))), axis=1)
        ci.columns = ['low', 'high']

        # Add sample number
        ci['samples'] = np.array(range(len(samples))) + 1

        # Plot
        fig, axis = plt.subplots(1)
        axis.plot(ci['samples'], ci['low'])
        axis.plot(ci['samples'], ci['high'])
        axis.set_xlabel("# samples")
        axis.set_ylabel("# reads")
        plt.legend(loc='best')
        fig.savefig(os.path.join(self.plots_dir, "cll_samples.cum_read_count.95CI.svg"), bbox_inches="tight")

        # Store
        df.to_csv(os.path.join(self.data_dir, "cll_samples.cum_read_count.95CI.csv"), index=False)

    def calculate_peak_support(self, samples):
        samples = [s for s in self.samples if s.library == "ATAC-seq" and s.cell_line == "CLL"]
        # calculate support (number of samples overlaping each merged peak)
        for i, sample in enumerate(samples):
            print(sample.name)
            if i == 0:
                support = self.sites.intersect(sample.peaks, wa=True, c=True)
            else:
                support = support.intersect(sample.peaks, wa=True, c=True)
        support = support.to_dataframe()
        support = support.reset_index()
        support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support[range(len(samples))].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(self.samples)), axis=1)
        # save
        support.to_csv(os.path.join(self.data_dir, "cll_peaks.support.csv"), index=False)

        self.support = support

    def get_peak_gene_annotation(self):
        """
        Annotates peaks with closest gene.
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with hg19 TSS positions
        hg19_ensembl_tss = pybedtools.BedTool(os.path.join(self.data_dir, "ensembl_tss.bed"))
        # get closest TSS of each cll peak
        closest = self.sites.closest(hg19_ensembl_tss, d=True).to_dataframe()[['chrom', 'start', 'end', 'thickStart', 'blockCount']]

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        gene_annotation = closest.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        gene_annotation.columns = ['chrom', 'start', 'end', 'ensembl_transcript_id', 'distance']

        # add gene name and ensemble_gene_id
        ensembl_gtn = pd.read_table(os.path.join(self.data_dir, "ensemblToGeneName.txt"), header=None)
        ensembl_gtn.columns = ['ensembl_transcript_id', 'gene_name']

        ensembl_gtp = pd.read_table(os.path.join(self.data_dir, "ensGtp.txt"), header=None)[[0, 1]]
        ensembl_gtp.columns = ['ensembl_gene_id', 'ensembl_transcript_id']

        gene_annotation = pd.merge(gene_annotation, ensembl_gtn, how="left")
        self.gene_annotation = pd.merge(gene_annotation, ensembl_gtp, how="left")

        # save to disk
        self.gene_annotation.to_csv(os.path.join(self.data_dir, "cll_peaks.gene_annotation.csv"), index=False)

        # save distances to all TSSs (for plotting)
        self.closest_tss_distances = closest['blockCount'].tolist()
        pickle.dump(self.closest_tss_distances, open(os.path.join(self.data_dir, "cll_peaks.closest_tss_distances.pickle"), 'wb'))

    def get_peak_genomic_location(self):
        """
        Annotates peaks with its type of genomic location.
        Needs files downloaded by prepare_external_files.py
        """
        regions = [
            "ensembl_genes.bed", "ensembl_tss2kb.bed",
            "ensembl_utr5.bed", "ensembl_exons.bed", "ensembl_introns.bed", "ensembl_utr3.bed"]

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        for i, region in enumerate(regions):
            region_name = region.replace(".bed", "").replace("ensembl_", "")
            r = pybedtools.BedTool(os.path.join(self.data_dir, region))
            if region_name == "genes":
                region_name = "intergenic"
                df = self.sites.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
                dfb = background.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
            else:
                df = self.sites.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
                dfb = background.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
            df['genomic_region'] = region_name
            dfb['genomic_region'] = region_name
            if i == 0:
                region_annotation = df
                region_annotation_b = dfb
            else:
                region_annotation = pd.concat([region_annotation, df])
                region_annotation_b = pd.concat([region_annotation_b, dfb])

        # sort
        region_annotation.sort(['chrom', 'start', 'end'], inplace=True)
        region_annotation_b.sort(['chrom', 'start', 'end'], inplace=True)
        # remove duplicates (there shouldn't be anyway)
        region_annotation = region_annotation.reset_index(drop=True).drop_duplicates()
        region_annotation_b = region_annotation_b.reset_index(drop=True).drop_duplicates()
        # join various regions per peak
        self.region_annotation = region_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        self.region_annotation_b = region_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.region_annotation.to_csv(os.path.join(self.data_dir, "cll_peaks.region_annotation.csv"), index=False)
        self.region_annotation_b.to_csv(os.path.join(self.data_dir, "cll_peaks.region_annotation_background.csv"), index=False)

    def get_peak_chromatin_state(self):
        """
        Annotates peaks with chromatin states.
        (For now states are from CD19+ cells).
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with CD19 chromatin states
        states_cd19 = pybedtools.BedTool(os.path.join(self.data_dir, "E032_15_coreMarks_mnemonics.bed"))

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        # intersect with cll peaks, to create annotation, get original peaks
        chrom_state_annotation = self.sites.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        chrom_state_annotation_b = background.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

        # aggregate annotation per peak, concatenate various annotations (comma-separated)
        self.chrom_state_annotation = chrom_state_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation.columns = ['chrom', 'start', 'end', 'chromatin_state']

        self.chrom_state_annotation_b = chrom_state_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation_b.columns = ['chrom', 'start', 'end', 'chromatin_state']

        # save to disk
        self.chrom_state_annotation.to_csv(os.path.join(self.data_dir, "cll_peaks.chromatin_state.csv"), index=False)
        self.chrom_state_annotation_b.to_csv(os.path.join(self.data_dir, "cll_peaks.chromatin_state_background.csv"), index=False)

    @pickle_me
    def measure_coverage(self, samples):
        # Count reads with pysam
        # make strings with intervals
        sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in self.sites]
        # count, create dataframe
        self.coverage = pd.DataFrame(
            map(
                lambda x:
                    pd.Series(x),
                    parmap.map(
                        count_reads_in_intervals,
                        [sample.filtered for sample in samples],
                        sites_str,
                        parallel=True
                    )
            ),
            index=[sample.name for sample in samples]
        ).T

        # Add interval description to df
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            self.coverage.index
        )
        self.coverage["chrom"] = [x[0] for x in ints]
        self.coverage["start"] = [int(x[1]) for x in ints]
        self.coverage["end"] = [int(x[2]) for x in ints]

        # save to disk
        self.coverage.to_csv(os.path.join(self.data_dir, "cll_peaks.raw_coverage.tsv"), sep="\t", index=True)

    @pickle_me
    def normalize_coverage_quantiles(self, samples):
        # Normalize by quantiles
        to_norm = self.coverage.iloc[:, :len(samples)]
        self.coverage_qnorm = pd.DataFrame(
            normalize_quantiles_r(np.array(to_norm)),
            index=to_norm.index,
            columns=to_norm.columns
        )
        # Log2 transform
        self.coverage_qnorm = np.log2(1 + self.coverage_qnorm)

        self.coverage_qnorm = self.coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        self.coverage_qnorm.to_csv(os.path.join(self.data_dir, "cll_peaks.coverage_qnorm.log2.tsv"), sep="\t", index=False)

    @pickle_me
    def annotate(self, samples):
        # add closest gene
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm,
            self.gene_annotation[['chrom', 'start', 'end', 'gene_name']], on=['chrom', 'start', 'end'])
        # add genomic location
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.region_annotation[['chrom', 'start', 'end', 'genomic_region']], on=['chrom', 'start', 'end'])
        # add chromatin state
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.chrom_state_annotation[['chrom', 'start', 'end', 'chromatin_state']], on=['chrom', 'start', 'end'])
        # add support to coverage - this is added here because support was calculated prior to coverage
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        # calculate mean coverage
        self.coverage_qnorm_annotated['mean'] = self.coverage_qnorm_annotated[[sample.name for sample in samples]].apply(lambda x: np.mean(x), axis=1)
        # calculate coverage variance
        self.coverage_qnorm_annotated['variance'] = self.coverage_qnorm_annotated[[sample.name for sample in samples]].apply(lambda x: np.var(x), axis=1)
        # calculate std deviation (sqrt(variance))
        self.coverage_qnorm_annotated['std_deviation'] = np.sqrt(self.coverage_qnorm_annotated['variance'])
        # calculate dispersion (variance / mean)
        self.coverage_qnorm_annotated['dispersion'] = self.coverage_qnorm_annotated['variance'] / self.coverage_qnorm_annotated['mean']
        # calculate qv2 (std / mean) ** 2
        self.coverage_qnorm_annotated['qv2'] = (self.coverage_qnorm_annotated['std_deviation'] / self.coverage_qnorm_annotated['mean']) ** 2

        # calculate "fold-change" (max - min)
        self.coverage_qnorm_annotated['fold_change'] = self.coverage_qnorm_annotated[[sample.name for sample in samples]].apply(lambda x: x.max() - x.min(), axis=1)

        self.coverage_qnorm_annotated.to_csv(os.path.join(self.data_dir, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index=False)

    def variability(self):
        # Variable region analysis
        df = self.coverage_qnorm_annotated
        # separate in three groups
        out = dict()
        out[1] = df[df['dispersion'] < 0.35]
        out[2] = df[
            (df['dispersion'] > 0.35) &
            (df['dispersion'] < 0.9)]
        out[3] = df[df['dispersion'] > 0.9]
        out[4] = df[df['dispersion'] > 1.5]

        # get regions of groups, write bed file
        out[1][['chrom', 'start', 'end']].to_csv("f1.bed", index=False, sep="\t", header=False)
        out[2][['chrom', 'start', 'end']].to_csv("f2.bed", index=False, sep="\t", header=False)
        out[3][['chrom', 'start', 'end']].to_csv("f3.bed", index=False, sep="\t", header=False)
        out[4][['chrom', 'start', 'end']].to_csv("f4.bed", index=False, sep="\t", header=False)

        # universe
        df[['chrom', 'start', 'end']].to_csv("universe.bed", index=False, sep="\t", header=False)

        # run lola
        lola("f4.bed", "universe.bed", "data/variability-4")
        lola("f3.bed", "universe.bed", "data/variability-3")
        lola("f2.bed", "universe.bed", "data/variability-2")
        lola("f1.bed", "universe.bed", "data/variability-1")

        # get genes of groups
        out2 = dict()
        out2[1] = out[1]['gene_name'].unique()
        out2[2] = out[2]['gene_name'].unique()
        out2[3] = out[3]['gene_name'].unique()
        out2[4] = out[4]['gene_name'].unique()
        # write gene names to file
        for i, d in out2.items():
            with open("f%i-genes.txt" % i, 'w') as handle:
                handle.writelines("\n".join(d.tolist()))

    def gene_oppeness_across_samples(self, samples):
        """
        Annotates peaks with closest gene.
        Needs files downloaded by prepare_external_files.py
        """
        from collections import OrderedDict
        from scipy.stats import ks_2samp

        sns.set(style="white", palette="pastel", color_codes=True)

        # genes of interest
        sel_genes = OrderedDict({
            # main B cell surface markers
            "CD19": "ENSG00000177455",
            "CD21": "ENSG00000117322",
            "CD22": "ENSG00000012124",
            "FCGR2B": "ENSG00000072694",  # CD32

            # other B cell surface markers
            "CD20": "ENSG00000156738",
            "CD24": "ENSG00000272398",
            "CD38": "ENSG00000004468",
            "CD72": "ENSG00000137101",
            "CD74": "ENSG00000019582",  # MHC2
            "CD79A": "ENSG00000105369",
            "CD79B": "ENSG00000007312",
            "CD83": "ENSG00000112149",
            "CD86": "ENSG00000114013",

            # signal transduction molecules
            "SYK": "ENSG00000165025",
            "LYN": "ENSG00000254087",
            "BTK": "ENSG00000010671",
            "BLNK": "ENSG00000095585",
            "BLK": "ENSG00000136573",

            # signaling important for B cells
            "NOTCH1": "ENSG00000148400",
            "NFKB1": "ENSG00000109320",

            # downstream of BCR signaling / proliferation
            "AKT1": "ENSG00000142208",
            "AIMP2": "ENSG00000106305",  # p38
            "mTOR": "ENSG00000198793",
            "ERK1": "ENSG00000102882",
            "PIK3CA": "ENSG00000121879",

            #
            "MYC": "ENSG00000136997",
            "MYCN": "ENSG00000134323",

            # transcriptional regulators
            "EBF1": "ENSG00000164330",
            "PAX5": "ENSG00000196092",
            "POU2AF1": "ENSG00000110777",
            "SPIB": "ENSG00000269404",
            "SPI1": "ENSG00000066336",
            "IRF4": "ENSG00000137265",
            "IRF8": "ENSG00000140968",
            "CEBPB": "ENSG00000172216",
            "BCL6": "ENSG00000113916",

            # others
            "ZAP70": "ENSG00000115085",
            "IL2": "ENSG00000109471",

            # CLL mutated genes
            # notch pathway
            "FBXW7": "ENSG00000109670",
            "SPEN": "ENSMUSG00000040761",
            "CREBBP": "ENSG00000005339",
            # b cell signaling
            "TLR2": "ENSG00000137462",
            "BCOR": "ENSG00000183337",
            "KLHL6": "ENSG00000172578",
            "IKZF3": "ENSG00000161405",
            # dna repair
            "ATM": "ENSG00000149311",
            "ATR": "ENSG00000175054",
            "POT1": "ENSG00000128513",
            "TP53": "ENSG00000141510",
            # chromatin
            "ARID1A": "ENSG00000117713",
            "SETD1A": "ENSG00000099381",
            "HIST1H1B": "ENSG00000168298",
            "ZMYM3": "ENSG00000147130",
            "SETD2": "ENSG00000181555",
            "KMT2D": "ENSG00000167548",  # MLL2
            "CHD2": "ENSG00000173575",
            "SYNE1": "ENSG00000234577",
            "ASXL1": "ENSG00000171456",
            # cell cycle
            "PTPN11": "ENSG00000179295",
            "KRAS": "ENSG00000133703",
            "NRAS": "ENSG00000213281",
            "BRAF": "ENSG00000157764",
            "CDKN1B": "ENSG00000111276",
            "CDKN2A": "ENSG00000147889",
            "CCND2": "ENSG00000118971",
            # apoptosis
            "BAX": "ENSG00000087088",
            "ANKHD1": "ENSG00000254996 ",
            # rna metabolism
            "MGA": "ENSG00000174197",
            "CNOT3": "ENSG00000088038 ",
            "MED12": "ENSG00000184634 ",
            "NXF1": "ENSG00000162231 ",
            "ZNF292": "ENSG00000188994",
            "SF3B1": "ENSG00000115524",
            "DDX3X": "ENSG00000215301",
            "XPO1": "ENSG00000082898",
            "TRAF3": "ENSG00000131323",
            "BIRC3": "ENSG00000023445",
            "NFKB2": "ENSG00000077150 ",
            "EGR2": "ENSG00000122877 ",
            "NFKBIE": "ENSG00000146232",
            "NKAP": "ENSG00000101882",
        })

        # get distance to gene and ensembl gene id annotation in whole matrix
        df = pd.merge(self.coverage_qnorm_annotated, self.gene_annotation, on=['chrom', 'start', 'end', 'gene_name'])

        # GET 1(!) element per gene
        # get peaks around promoter (+/- 1kb)
        df2 = df[df['distance'] <= 2500]
        # promoters
        promoter_index = df2.groupby(["ensembl_gene_id"]).apply(lambda x: np.argmin((x['distance'])))
        promoters = df2.ix[promoter_index]

        # get peaks away from promoters (> 1kb)
        df2 = df[df['distance'] > 2500]
        # enhancers
        enhancer_index = df2.groupby(["ensembl_gene_id"]).apply(lambda x: np.argmin((x['distance'])))
        enhancers = df2.ix[enhancer_index]

        # Figure 2a - variability
        # 81 genes on top of all genes
        genes_str = "|".join(sel_genes.values())
        p_values = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)]['variance']
        e_values = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)]['variance']

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.distplot(np.log2(1 + promoters['variance']), ax=axis[0], bins=100)
        sns.distplot(np.log2(1 + enhancers['variance']), ax=axis[1], bins=100)
        sns.distplot(np.log2(1 + p_values), ax=axis[0], rug=True, bins=20)
        sns.distplot(np.log2(1 + e_values), ax=axis[1], rug=True, bins=20)
        fig.savefig("prom-enh.variance.log2.svg", bbox_inches="tight")

        # test difference in distributions
        D, p = ks_2samp(promoters['variance'], p_values)
        D, p = ks_2samp(enhancers['variance'], e_values)

        # Plot distributions of amplitude (fold_change)
        fig, axis = plt.subplots(1, figsize=(15, 10))
        sns.distplot(promoters['amplitude'], color="b", ax=axis)
        sns.distplot(enhancers['amplitude'], color="y", ax=axis)
        fig.savefig(os.path.join("results", "plots", "all_genes.distplot.svg"), bbox_inches='tight')

        # plot aditional boxplots for selected genes
        gene_values = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)][[sample.name for sample in samples]].T
        gene_values.columns = promoters.ix[gene_values.columns]['gene_name']
        promoter_data = pd.melt(gene_values, var_name="gene", value_name="openness")
        promoter_data['region'] = 'promoter'

        gene_values = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)][[sample.name for sample in samples]].T
        gene_values.columns = enhancers.ix[gene_values.columns]['gene_name']
        enhancer_data = pd.melt(gene_values, var_name="gene", value_name="openness")
        enhancer_data['region'] = 'enhancer'

        boxplot_data = pd.concat([promoter_data, enhancer_data])

        fig, axis = plt.subplots(1, figsize=(45, 10))
        sns.violinplot(data=boxplot_data.sort('openness'), x="gene", y="openness", hue="region", split=True, inner="quart", palette={"promoter": "b", "enhancer": "y"}, jitter=True, ax=axis)
        fig.savefig(os.path.join("results", "plots", "relevant_genes.full.violinplot.svg"), bbox_inches='tight')

        # sort by predefined order (intensity/functional classes)
        sorterIndex = dict(zip(sel_genes.keys(), range(len(sel_genes.keys()))))
        boxplot_data['order'] = boxplot_data['gene'].map(sorterIndex)
        boxplot_data.sort('order', ascending=False, inplace=True)

        fig, axis = plt.subplots(1, figsize=(45, 10))
        sns.violinplot(data=boxplot_data, x="gene", y="openness", hue="region", split=True, inner="quart", palette={"promoter": "b", "enhancer": "y"}, jitter=True, ax=axis)
        fig.savefig(os.path.join("results", "plots", "relevant_genes.full.violinplot.funcorder.svg"), bbox_inches='tight')

    def inspect_variability(self, samples):
        """
        Investigate variability within sample groups.
        """
        from statsmodels.sandbox.stats.multicomp import multipletests

        def f_test_pvalue_rpy2(d1, d2):
            import rpy2.robjects as robjects

            rd1 = (robjects.FloatVector(d1))
            rd2 = (robjects.FloatVector(d2))
            rvtest = robjects.r['var.test']
            return rvtest(rd1, rd2)[2][0]

        # plot mean vs qv2
        plt.scatter(self.coverage_qnorm_annotated["mean"], self.coverage_qnorm_annotated["qv2"])
        plt.savefig(os.path.join(self.plots_dir, "mean_qv2.svg"), bbox_inches="tight")
        # plot mean vs dispersion
        plt.scatter(self.coverage_qnorm_annotated["mean"], self.coverage_qnorm_annotated["std"])
        plt.savefig(os.path.join(self.plots_dir, "mean_dispersion.svg"), bbox_inches="tight")

        # divide samples per IGHV status
        ighv_u = [s.name for s in samples if not s.mutated]
        ighv_m = [s.name for s in samples if s.mutated]

        df_u = self.coverage_qnorm_annotated[ighv_u]
        df_m = self.coverage_qnorm_annotated[ighv_m]

        # for each element calculate variability
        u = pd.DataFrame()
        m = pd.DataFrame()

        u["mean"] = df_m.apply(np.mean, axis=1)
        u["std"] = df_m.apply(np.std, axis=1)
        u["var"] = df_m.apply(np.var, axis=1)
        u["dispersion"] = u["var"] / u["mean"]
        u["qv2"] = (u["std"] / u["mean"]) ** 2

        m["mean"] = df_u.apply(np.mean, axis=1)
        m["std"] = df_u.apply(np.std, axis=1)
        m["var"] = df_u.apply(np.var, axis=1)
        m["dispersion"] = m["var"] / m["mean"]
        m["qv2"] = (m["std"] / m["mean"]) ** 2

        u_m = pd.melt(u)
        u_m["subtype"] = "u"
        m_m = pd.melt(m)
        m_m["subtype"] = "m"

        df = u_m.append(m_m, ignore_index=True)
        # plot mean, dispersion, variance, qv2 in uCLL vs mCLL
        # boxplot
        g = sns.FacetGrid(df, col="variable", sharey=False)
        g.map(sns.violinplot, "subtype", "value", palette="pastel")
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.violin.svg"), bbox_inches="tight")
        plt.close("all")

        # hexbin plots
        g = sns.FacetGrid(df, col="variable", sharex=False, sharey=False)
        for i, variable in enumerate(df["variable"].unique()):
            a = df[
                (df["variable"] == variable) &
                (df["subtype"] == "u")
            ]["value"]
            b = df[
                (df["variable"] == variable) &
                (df["subtype"] == "m")
            ]["value"]
            if i > 1:
                a = np.log2(1 + a)
                b = np.log2(1 + b)

            g.axes[0][i].hexbin(
                a,
                b,
                bins="log"
            )
            # x=y line
            lims = [
                np.min([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # min of both axes
                np.max([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # max of both axes
            ]
            g.axes[0][i].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
            g.axes[0][i].plot()
            g.axes[0][i].set_aspect('equal')
            g.axes[0][i].set_xlim(lims)
            g.axes[0][i].set_ylim(lims)
        g.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.hexbin.svg"), bbox_inches="tight")
        plt.close("all")

        # get significantly variable regions
        # F-test
        p_values = dict()
        for i in df_u.index:
            p_values[i] = f_test_pvalue_rpy2(df_u.ix[i], df_m.ix[i])
        # plot uncorrected p_values
        sns.distplot(-np.log10(np.array(p_values.values())), hist=False)
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.p_values.svg"), bbox_inches="tight")
        plt.close("all")
        # volcano plot
        plt.scatter(np.log2(u['dispersion'] / m['dispersion']), -np.log10(np.array(p_values.values())))
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.volcano.svg"), bbox_inches="tight")
        plt.close("all")

        # correct p-values
        p_values = dict(zip(p_values.keys(), multipletests(p_values.values())[1]))
        # plot p-value distribution
        sns.distplot(-np.log10(np.array(p_values.values())), hist=False)
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.p_values-corrected.svg"), bbox_inches="tight")
        plt.close("all")
        # volcano plot
        plt.scatter(np.log2(u['dispersion'] / m['dispersion']), -np.log10(np.array(p_values.values())))
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.volcano-corrected.svg"), bbox_inches="tight")
        plt.close("all")

        # significantly variable regions
        # by threshold
        indexes = [i for i, p in p_values.items() if -np.log10(p) > 1.3]

        # by n-most variable
        n = 1000
        indexes = [i for i, k in sorted(p_values.items(), key=lambda x: x[1])[:n]]

        # filter out regions with mean across all samples lower than 1
        indexes = [i for i in indexes if self.coverage_qnorm_annotated.ix[i]["mean"] > 1]

        # hexbin plots
        g = sns.FacetGrid(df, col="variable", sharex=False, sharey=False)
        for i, variable in enumerate(df["variable"].unique()):
            a = df[
                (df["variable"] == variable) &
                (df["subtype"] == "u")
            ]["value"]
            b = df[
                (df["variable"] == variable) &
                (df["subtype"] == "m")
            ]["value"]
            if i > 1:
                a = np.log2(1 + a)
                b = np.log2(1 + b)

            g.axes[0][i].hexbin(
                a,
                b,
                bins="log"
            )

            # significant as scatter
            if i > 1:
                g.axes[0][i].scatter(np.log2(1 + u.ix[indexes][variable]), np.log2(1 + m.ix[indexes][variable]), color="orange")
            else:
                g.axes[0][i].scatter(u.ix[indexes][variable], m.ix[indexes][variable], color="orange")

            # x=y line
            lims = [
                np.min([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # min of both axes
                np.max([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # max of both axes
            ]
            g.axes[0][i].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
            g.axes[0][i].set_aspect('equal')
            # g.axes[0][i].set_xlim(lims)
            # g.axes[0][i].set_ylim(lims)

            g.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.hexbin.significant.svg"), bbox_inches="tight")
            plt.close("all")

        # investigate:
        # get bed files for u/mCLL specifically
        v_u = list()
        v_m = list()
        for i in indexes:
            if (self.coverage_qnorm_annotated.ix[i]["mean"] > 1):
                if u.ix[i]['dispersion'] > m.ix[i]['dispersion']:
                    v_u.append(i)
                else:
                    v_m.append(i)
        self.coverage_qnorm_annotated.ix[v_u][['chrom', 'start', 'end']].to_csv("ighv_mutation_variable_comparison.significant-uCLL.bed", index=False, sep="\t", header=False)
        self.coverage_qnorm_annotated.ix[v_m][['chrom', 'start', 'end']].to_csv("ighv_mutation_variable_comparison.significant-mCLL.bed", index=False, sep="\t", header=False)

        universe = "data/cll_peaks.bed"

        # run lola
        lola("ighv_mutation_variable_comparison.significant-uCLL.bed", universe, "data/ighv_mutation_variable_comparison/uCLL")
        lola("ighv_mutation_variable_comparison.significant-mCLL.bed", universe, "data/ighv_mutation_variable_comparison/mCLL")

        # get gene names
        out = dict()
        out['u'] = self.coverage_qnorm_annotated.ix[v_u]['gene_name'].unique()
        out['m'] = self.coverage_qnorm_annotated.ix[v_m]['gene_name'].unique()
        # write gene names to file
        for c, genes in out.items():
            with open("data/ighv_mutation_variable_comparison/%sCLL/genes.txt" % c, 'w') as handle:
                handle.writelines("\n".join(genes.tolist()))

    def correlate_expression(self):
        # get expression
        expression_matrix = pd.read_csv(os.path.join("data", "CLL.geneReadcount.txt"), sep=" ")
        expression_matrix.index = expression_matrix['geneid']
        # get values from numbered samples (n=98)
        expression_values = expression_matrix[[n for n in expression_matrix.columns if "X" == n[0]]]
        # average across all samples
        expression_mean = expression_values.apply(np.mean, axis=1).reset_index()
        expression_mean.columns = ['ensembl_gene_id', 'rna']
        # log expression
        expression_mean['rna'] = np.log2(1 + expression_mean['rna'])
        # get only expressed genes
        # expression_mean = expression_mean[expression_mean['rna'] > 3]

        # get oppenness
        # get only autosomes
        openness = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated['chrom'].str.contains("chr[^X|Y]")]

        # get closest gene info with ensembl ids
        g = pd.read_csv("data/cll_peaks.gene_annotation.csv")
        # get only genes within 5kb
        g = g[g['distance'] < 5000]
        openness = pd.merge(openness, g)
        # weight atac signal with distance
        # openness['weighted_mean'] = openness.apply(lambda x: x['mean'] * np.exp(-x['distance'] / 1000.), axis=1)

        # average 'oppenness' for various sites for each gene
        openness = openness[['mean', 'ensembl_gene_id']].groupby('ensembl_gene_id').aggregate(np.mean).reset_index()
        openness.columns = ['ensembl_gene_id', 'atac']

        # merge expression and openness
        m = pd.merge(expression_mean, openness)

        # plot
        sns.jointplot(m['atac'], m['rna'], kind='scatter', s=3, linewidth=1, alpha=.3)
        plt.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.scatter.svg"), bbox_inches="tight")
        sns.jointplot(m['atac'], m['rna'], kind='kde')
        plt.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.kde.svg"), bbox_inches="tight")

        # os.path.join("data", "ferreira_2012_gr.ighv_mutation_status.csv")

    def plot_peak_characteristics(self):
        # Loop at summary statistics:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites], bins=300, kde=False, ax=axis)
        axis.set_xlim(0, 2000)  # cut it at 2kb
        axis.set_xlabel("peak width (bp)")
        axis.set_ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.lengths.svg"), bbox_inches="tight")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis)
        axis.set_ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.support.svg"), bbox_inches="tight")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots()
        sns.distplot(self.closest_tss_distances, bins=200, ax=axis)
        axis.set_xlim(0, 100000)  # cut it at 100kb
        axis.set_xlabel("distance to nearest TSS (bp)")
        axis.set_ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.tss_distance.svg"), bbox_inches="tight")

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = [item for sublist in self.region_annotation['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]
        all_region_annotation_b = [item for sublist in self.region_annotation_b['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_region_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_region_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        # plot individually
        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], (data[1] / background[1])]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("genomic region")
        axis[2].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.genomic_regions.svg"), bbox_inches="tight")

        # Plot chromatin states
        # get long list of chromatin states (for plotting)
        all_chrom_state_annotation = [item for sublist in self.chrom_state_annotation['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]
        all_chrom_state_annotation_b = [item for sublist in self.chrom_state_annotation_b['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_chrom_state_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_chrom_state_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], (data[1] / background[1])]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("chromatin state")
        axis[2].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.chromatin_states.svg"), bbox_inches="tight")

        # distribution of count attributes
        data = self.coverage_qnorm_annotated.copy()

        sns.distplot(data["mean"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.mean.distplot.svg"), bbox_inches="tight")
        plt.close()

        sns.distplot(data["qv2"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.qv2.distplot.svg"), bbox_inches="tight")
        plt.close()

        sns.distplot(data["dispersion"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.dispersion.distplot.svg"), bbox_inches="tight")
        plt.close()

        # this is loaded now
        df = pd.read_csv(os.path.join(self.data_dir, "cll_peaks.support.csv"))
        sns.distplot(df["support"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.support.distplot.svg"), bbox_inches="tight")
        plt.close()

    def plot_coverage(self):
        data = self.coverage_qnorm_annotated.copy()
        # (rewrite to avoid putting them there in the first place)
        for variable in ['gene_name', 'genomic_region', 'chromatin_state']:
            d = data[variable].str.split(',').apply(pd.Series).stack()  # separate comma-delimited fields
            d.index = d.index.droplevel(1)  # returned a multiindex Series, so get rid of second index level (first is from original row)
            data = data.drop([variable], axis=1)  # drop original column so there are no conflicts
            d.name = variable
            data = data.join(d)  # joins on index

        # Plot
        data_melted = pd.melt(
            data,
            id_vars=["chrom", "start", "end", 'mean', 'dispersion', 'qv2',
                     'gene_name', 'genomic_region', 'chromatin_state', 'support'], var_name="sample", value_name="norm_counts")

        # Together in same violin plot
        sns.violinplot("genomic_region", "norm_counts", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.per_genomic_region.violinplot.svg"), bbox_inches="tight")
        plt.close()
        sns.violinplot("genomic_region", "norm_counts", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.chromatin_state.violinplot.svg"), bbox_inches="tight")
        plt.close()
        # dispersion
        sns.violinplot("genomic_region", "dispersion", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.per_genomic_region.violinplot.svg"), bbox_inches="tight")
        plt.close()
        sns.violinplot("genomic_region", "dispersion", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.chromatin_state.violinplot.svg"), bbox_inches="tight")
        plt.close()
        # dispersion
        sns.violinplot("genomic_region", "dispersion", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.per_genomic_region.violinplot.svg"), bbox_inches="tight")
        plt.close()
        sns.violinplot("genomic_region", "dispersion", data=data_melted)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.chromatin_state.violinplot.svg"), bbox_inches="tight")
        plt.close()

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.mean.per_genomic_region.distplot.svg"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.mean.chromatin_state.distplot.svg"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.per_genomic_region.distplot.svg"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.dispersion.chromatin_state.distplot.svg"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.support.per_genomic_region.distplot.svg"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts.support.chromatin_state.distplot.svg"), bbox_inches="tight")
        plt.close()

    def plot_variance(self):
        sns.jointplot('mean', "dispersion", data=self.coverage_qnorm_annotated, kind="kde")
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.dispersion.svg"), bbox_inches="tight")
        plt.close("all")

        sns.jointplot('mean', "qv2", data=self.coverage_qnorm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.qv2_vs_mean.svg"), bbox_inches="tight")
        plt.close("all")

        sns.jointplot('support', "qv2", data=self.coverage_qnorm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.support_vs_qv2.svg"), bbox_inches="tight")
        plt.close("all")

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated[[sample.name for sample in self.samples]].apply(max, axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.dispersion.filtered.svg"), bbox_inches="tight")
        plt.close("all")
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.support_vs_qv2.filtered.svg"), bbox_inches="tight")

    def plot_qnorm_comparison(self):
        # Compare raw counts vs qnormalized data
        fig, axis = plt.subplots(2, sharex=True)
        [sns.distplot(np.log2(1 + self.coverage[[sample.name]]), ax=axis[0], hist=False) for sample in self.samples if sample.cell_line != "PBMC"]
        [sns.distplot(self.coverage_qnorm[[sample.name]], ax=axis[1], hist=False) for sample in self.samples if sample.cell_line != "PBMC"]
        axis[0].set_title("Raw counts")
        axis[1].set_title("Quantile normalized counts")
        axis[1].set_xlabel("log2(1 + x)")
        fig.savefig(os.path.join(self.plots_dir, "coverage_vs_coverage_qnorm.svg"), bbox_inches="tight")

    def plot_qv2_fit(self):
        from scipy.optimize import curve_fit
        from scipy import stats

        def neg_exponential(x, a, b, c):
            """
            Negative exponential function with 3 parameters.
            """
            return a * np.exp(-b * x) + c

        X = np.array(self.rpkm['mean'])
        Y = np.array(self.rpkm['qv2'])

        ci = 0.99
        # Convert to percentile point of the normal distribution.
        # See: https://en.wikipedia.org/wiki/Standard_score
        pp = (1. + ci) / 2.

        nstd = stats.norm.ppf(pp)

        # Find best fit.
        parameters, covariance_matrix = curve_fit(neg_exponential, X, Y)
        # Standard deviation errors on the parameters.
        perr = np.sqrt(np.diag(covariance_matrix))
        # Add nstd standard deviations to parameters to obtain the upper confidence
        # interval.
        popt_up = parameters + (nstd * perr)
        popt_dw = parameters - (nstd * perr)

        fig, axis = plt.subplots(2, sharex=True)
        # Plot data and best fit curve.
        axis[0].scatter(X, Y)
        x = np.linspace(0, 6.5, 100)
        axis[0].plot(x, neg_exponential(x, *parameters), c='g', lw=2.)
        axis[0].plot(x, neg_exponential(x, *popt_up), c='r', lw=2.)
        axis[0].plot(x, neg_exponential(x, *popt_dw), c='r', lw=2.)
        axis[0].set_title("fit")

        # get residuals
        residuals = Y - neg_exponential(X, *parameters)
        # get squared sum of residuals
        # fres = sum(residuals ** 2)

        axis[1].scatter(X, residuals)
        axis[1].set_title("residuals")

        axis[0].set_xlabel("mean")
        axis[0].set_ylabel("qv2")
        axis[1].set_xlabel("mean")
        axis[1].set_ylabel("residuals")

        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.qv2_vs_mean.fit_residuals.svg"), bbox_inches="tight")

    def plot_sample_correlations(self):
        # get colors depending on IGHV mut
        df = pd.DataFrame([sample.asSeries() for sample in self.samples])

        df = pd.merge(df, self.clinical, left_on="sample_id", right_on="sample_id")
        mut_s = {"1.0": "red", "2.0": "black", "nan": "grey"}
        sex_s = {"F": "red", "M": "blue", "nan": "grey"}
        muts = [self.clinical.loc[self.clinical['sample_id'] == sample.sample_id, "igvh_mutation_status"] for sample in self.samples]
        sex = [self.clinical.loc[self.clinical['sample_id'] == sample.sample_id, "patient_gender"] for sample in self.samples]
        mut_colors = [mut_s[str(x.get(x.index[0]))] if len(x.tolist()) > 0 else "grey" for x in muts]
        sex_colors = [sex_s[str(x.get(x.index[0]))] if len(x.tolist()) > 0 else "grey" for x in sex]

        # Correlation
        cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

        correlations = self.rpkm[[sample.name for sample in self.samples]].corr()
        correlations.index = map(name_to_repr, correlations.index)
        correlations.columns = map(name_to_repr, correlations.columns)

        sns.clustermap(
            correlations,
            method="complete",
            cmap=cmap,
            row_colors=mut_colors,
            col_colors=sex_colors,
            square=True, annot=False,
            figsize=(20, 16)
        )
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.correlation_clustering.svg"), bbox_inches="tight")
        plt.close("all")


def add_args(parser):
    """
    Options for project and pipelines.
    """
    # Behaviour
    parser.add_argument("-g", "--generate", dest="generate", action="store_true",
                        help="Should we generate data and plots? Default=False")

    return parser


def count_reads_in_intervals(bam, intervals):
    """
    Counts reads in a iterable holding strings
    representing genomic intervals of the type chrom:start-end.
    """
    counts = dict()

    bam = pysam.Samfile(bam, 'rb')

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX"]

    for interval in intervals:
        if interval.split(":")[0] not in chroms:
            continue
        counts[interval] = bam.count(region=interval)
    bam.close()

    return counts


def normalize_by_interval_length(series):
    """
    Divides a pandas.Series with numerical values by the length of the interval
    (defined in "chr", "start", "end") in a bam file (RPK).
    """
    length = float(series["end"] - series["start"])

    rpk = (series.drop(["chrom", "start", "end"]) / length) * 1e3

    return series[["chrom", "start", "end"]].append(rpk)


def normalize_by_library_size(series, samples, rm_mt=True):
    """
    Divides a pandas.Series with numerical values by the total number of mapped
    reads in a bam file.
    """
    # get sample which this series belongs to
    sample = [sample for sample in samples if sample.name == series.name][0]

    # get number of aligned reads for that sample
    size = float(pysam.AlignmentFile(sample.filtered).mapped)

    if rm_mt:
        mt = pysam.AlignmentFile(sample.filtered).count(reference="chrM")
        return (series / (size - mt)) * 1e6
    else:
        return (series / size) * 1e6


def normalize_quantiles(array):
    """
    array with samples in the columns and probes across the rows
    """
    n_array = np.zeros_like(array)

    sort_order = np.argsort(array, axis=0)

    n_array[sort_order, np.arange(array.shape[1])] = np.mean(array[sort_order, np.arange(array.shape[1])], axis=1)[:, np.newaxis]

    return n_array


def normalize_quantiles_r(array):
    # install package
    # R
    # source('http://bioconductor.org/biocLite.R')
    # biocLite('preprocessCore')

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("preprocessCore")')
    normq = robjects.r('normalize.quantiles')
    return np.array(normq(array))


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    """This returns joined patient and sample IDs"""
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_patient_id(name):
    return name.split("_")[2]


def name_to_sample_id(name):
    return name.split("_")[3]


def samples_to_color(samples, method="mutation"):
    # dependent on ighv status
    if method == "mutation":
        # This uses sns colorblind pallete
        colors = list()
        for sample in samples:
            if sample.mutated is True:
                colors.append(sns.color_palette("colorblind")[0])  # blue #0072b2
            elif sample.mutated is False:
                colors.append(sns.color_palette("colorblind")[2])  # vermillion #d55e00
            elif sample.mutated is None:
                if sample.cell_line == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
        return colors
    # dependent on ighv status
    if method == "homology":
        # This uses sns summer colormap
        cmap = plt.get_cmap('summer')
        # scale colormap to min and max ighv homology
        norm = matplotlib.colors.Normalize(vmin=86, vmax=100)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = list()
        for sample in samples:
            if 86 < sample.ighv_homology < 100:
                colors.append(m.to_rgba(sample.ighv_homology))
            else:
                if sample.cell_line == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
                return colors
    # unique color per patient
    elif method == "unique":
        # per patient
        patients = set(sample.patient_id for sample in samples)
        color_dict = cm.Paired(np.linspace(0, 1, len(patients)))
        color_dict = dict(zip(patients, color_dict))
        return [color_dict[sample.patient_id] for sample in samples]
    # rainbow (unique color per sample)
    elif method == "unique_sample":
        return cm.Paired(np.linspace(0, 1, len(samples)))
    elif method == "gender":
        colors = list()
        for sample in samples:
            if sample.patient_gender is "F":
                colors.append('red')
            elif sample.patient_gender is "M":
                colors.append('blue')
            elif sample.patient_gender is None:
                if sample.cell_line == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
        return colors
    elif method == "treatment":
        drugs = ['Alemtuz', "Ibrutinib"]
        colors = list()
        for sample in samples:
            if sample.under_treatment is False and sample.relapse is False:
                colors.append('peru')  # untreated samples
            elif sample.under_treatment is True and sample.treatment_regimen not in drugs:
                colors.append('black')  # chemo
            elif sample.under_treatment is True and sample.treatment_regimen in drugs:
                colors.append('green')  # antibodies
            elif sample.under_treatment is True and sample.treatment_regimen == "":
                colors.append('grey')  # no info
            else:
                colors.append('grey')  # no info
        return colors
    elif method == "disease":
        colors = list()
        for sample in samples:
            if sample.diagnosis_disease == "CLL":
                colors.append('#A6CEE3')
            elif sample.diagnosis_disease == "MBL":
                colors.append('#F17047')
            elif sample.diagnosis_disease == "SLL":
                colors.append('#482115')
            else:
                colors.append('grey')
        return colors
    else:
        raise ValueError("Method %s is not valid" % method)


def samples_to_symbol(samples, method="unique"):
    from itertools import cycle
    valid = ['D', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v']
    c = cycle([x for x in matplotlib.markers.MarkerStyle.markers.items() if x[0] in valid])

    # unique color per patient
    if method == "unique":
        # per patient
        patients = set(sample.patient_id for sample in samples)
        symbol_dict = [c.next()[0] for _ in range(len(patients))]
        symbol_dict = dict(zip(patients, symbol_dict))
        return [symbol_dict[sample.patient_id] for sample in samples]
    # rainbow (unique color per sample)
    elif method == "unique_sample":
        return [c.next()[0] for sample in samples]
    else:
        raise ValueError("Method %s is not valid" % method)

# def all_sample_colors(samples, order=""):
#     return [
#         samples_to_color(samples, "mutation"),
#         samples_to_color(samples, "homology"),
#         # samples_to_color(samples, "treatment"),
#         # samples_to_color(samples, "disease")
#         # samples_to_color(samples, "unique"),
#         samples_to_color(samples, "gender")
#     ]


def annotate_clinical_traits(samples):
    # Annotate traits
    # Gender
    t = [1 if s.patient_gender == "M" else 0 if s.patient_gender is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.gender = t[i]
    # IGHV mutation status
    t = [1 if s.igvh_mutation_status else 0 if s.igvh_mutation_status is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.IGHV = t[i]
    # CD38 expression
    t = [1 if s.CD38_positive else 0 if s.CD38_positive is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.CD38 = t[i]
    # ZAP70 expression
    t = [1 if s.ZAP70_positive else 0 if s.ZAP70_positive is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.ZAP70 = t[i]
    # ZAP70 mono-allelic expression
    t = [1 if s.ZAP70_monoallelic_methylation else 0 if s.ZAP70_monoallelic_methylation is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.ZAP70_monoallelic_methylation = t[i]
    # Disease at Diagnosis - comparison in untreated samples
    # Primary vs Secondary CLL (progressed from MBL and SLL)
    t = [1 if s.diagnosis_disease == "CLL" else 0 if s.diagnosis_disease is not pd.np.nan and not s.under_treatment and s.timepoint == 1 is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.primary_CLL = t[i]
    # Treatments: Under treatment
    t = [1 if s.under_treatment else 0 if s.under_treatment is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.treated = t[i]
    # Relapse
    # "relapse", ("True", "False", rerun=True), # relapse or before relapse
    t = [1 if s.relapse else 0 if s.relapse in [1, 0] is not pd.np.nan else pd.np.nan for s in samples]
    for i, s in enumerate(samples): s.relapse = t[i]
    # # Early vs Late
    # # "progression", ("True", "False", rerun=True), # early vs late timepoint
    # t = [1 if s.diagnosis_collection else 0 if type(s.diagnosis_collection) is bool is not pd.np.nan else pd.np.nan for s in samples]
    # for i, s in enumerate(samples): s.early = t[i]

    # Annotate samples which are under treament but with different types
    chemo_drugs = ["Chlor", "Chlor R", "B Of", "BR", "CHOPR"]  # Chemotherapy
    target_drugs = ["Alemtuz", "Ibrutinib"]  # targeted treatments
    muts = ["del13", "del11", "tri12", "del17"]  # chrom abnorms
    muts += ["SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "TP53", "MYD88", "CHD2", "NFKIE"]  # mutations
    for sample in samples:
        sample.chemo_treated = True if sample.under_treatment and sample.treatment_regimen in chemo_drugs else False
        sample.target_treated = True if sample.under_treatment and sample.treatment_regimen in target_drugs else False
        for mut in muts:
            setattr(sample, mut, True if sample.mutations is not pd.np.nan and mut in sample.mutations else False)

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
            sample.primary_CLL = True if sample.diagnosis_disease == "CLL" else False  # binary label useful for later

            # Get time since diagnosis
            sample.time_since_diagnosis = sample.collection_date - sample.diagnosis_date

            # Annotate treatment type, time since treatment
            if sample.under_treatment:
                sample.time_since_treatment = sample.collection_date - string_to_date(sample.treatment_date)

        # Append sample
        new_samples.append(sample)
    return new_samples


def annotate_samples(samples):
    new_samples = list()
    for sample in samples:
        # If any attribute is not set, set to NaN
        attrs = [
            "sample_name", "cell_line", "number_cells", "library", "ip", "patient_id", "sample_id", "condition",
            "experiment_name", "technical_replicate", "organism", "flowcell", "lane", "BSF_name", "clinical_centre", "sample_type",
            "cell_type", "sample_origin", "original_patient_id", "timepoint", "timepoint_name", "patient_gender", "patient_birth_date",
            "patient_death_date", "patient_last_checkup_date", "diagnosis_date", "diagnosis_disease", "hospital", "diagnosis_change_date",
            "diagnosis_stage_rai", "diagnosis_stage_binet", "under_treatment", "treatment_date", "treatment_regimen", "treatment_response",
            "treatment_end_date", "relapse", "treatment_1_date", "treatment_1_regimen", "treatment_1_duration", "treatment_1_response",
            "treatment_2_date", "treatment_2_regimen", "treatment_2_duration", "treatment_2_response", "treatment_3_date", "treatment_3_regimen",
            "treatment_3_duration", "treatment_3_response", "treatment_4_date", "treatment_4_regimen", "treatment_4_response", "treatment_end_date.1",
            "igvh_gene", "igvh_homology", "igvh_mutation_status", "CD38_cells_percentage", "CD38_positive", "CD38_measurement_date", "CD38_changed",
            "ZAP70_cells_percentage", "ZAP70_positive", "ZAP70_monoallelic_methylation", "sample_collection_date", "storage_condition", "lymp_count",
            "mutations", "sample_shippment_batch", "sample_cell_number", "sample_experiment_name", "sample_processing_date", "sample_viability"
            "diagnosis_collection", "diagnosis_date", "diagnosis_disease", "time_since_treatment", "treatment_regimen",
            "treatment_response", "treated", "previous_treatment_date", "previous_response", "relapse"]
        for attr in attrs:
            if not hasattr(sample, attr):
                setattr(sample, attr, pd.np.nan)
        new_samples.append(sample)

    # read in file with IGHV group of samples selected for ChIPmentation
    selected = pd.read_csv(os.path.join("metadata", "selected_samples.tsv"), sep="\t").astype(str)
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


def state_enrichment_overlap(n=100):
    cll_peaks = "~/cll_peaks.bed"
    all_states = "all_states_all_lines.bed"

    # states of interest:
    # get names of all states
    states = pd.read_csv(all_states, sep="\t", header=None)[3].unique().tolist()

    # loop through states, merge intervals, count number intersepting CLL peaks, and not intersepting
    cll_ints = pybedtools.BedTool(cll_peaks)

    df = pd.DataFrame()
    for state in states[-3:]:
        state_bed = "{0}.bed".format(state)
        os.system("grep {0} {1} > {2}".format(state, all_states, state_bed))

        # merge all intervals (of the same type across cell types)
        state_ints = pybedtools.BedTool(state_bed).sort().merge()

        total = len(state_ints)
        pos = len(state_ints.intersect(cll_ints))

        # get mean of `n` shuffled cll sites
        background = list()
        for i in range(n):
            background.append(len(state_ints.intersect(cll_ints.shuffle(genome='hg19', chrom=True))))

        # append to df
        df = df.append(pd.Series([total, pos, np.round(np.mean(background))]), ignore_index=True)
    df.index = states
    df.columns = ['total', 'pos', 'background']

    df['state'] = df.index

    df.to_csv("chrom_state_overlap_all.csv", index=False)

    df2 = pd.melt(df, id_vars='state')

    df2.sort(['variable', 'value'], inplace=True)

    fig, axis = plt.subplots(1)
    sns.barplot(data=df2, x='state', y='value', hue='variable', ax=axis)
    fig.savefig("chrom_state_overlap_all.svg", bbox_inches='tight')

    # fraction of total
    df['posF'] = df['pos'] / df['total']
    df['backgroundF'] = df['background'] / df['total']
    df3 = pd.melt(df[["state", "posF", "backgroundF"]], id_vars='state')
    df3.sort(['variable', 'value'], inplace=True)

    fig, axis = plt.subplots(1)
    sns.barplot(data=df3, x='state', y='value', hue='variable', ax=axis)
    fig.savefig("chrom_state_overlap_all.fraction_total.svg", bbox_inches='tight')

    # fraction of total enriched over background
    df['foldF'] = df['posF'] / df['backgroundF']
    df4 = pd.melt(df[["state", "foldF"]], id_vars='state')
    df4.sort(['variable', 'value'], inplace=True)

    fig, axis = plt.subplots(1)
    sns.barplot(data=df4, x='state', y='value', hue='variable', ax=axis)
    fig.savefig("chrom_state_overlap_all.fraction_total.enriched.svg", bbox_inches='tight')

    # same with log2
    df['foldFlog'] = np.log2(df['foldF'])
    df5 = pd.melt(df[["state", "foldFlog"]], id_vars='state')
    df5.sort(['variable', 'value'], inplace=True)

    fig, axis = plt.subplots(1)
    sns.barplot(data=df5, x='state', y='value', hue='variable', ax=axis)
    fig.savefig("chrom_state_overlap_all.fraction_total.enriched.log.svg", bbox_inches='tight')


def bed_to_fasta(bed_file, fasta_file):
    cmd = "bedtools getfasta -fi ~/resources/genomes/hg19/hg19.fa -bed {0} -fo {1}".format(bed_file, fasta_file)
    os.system(cmd)


def pca_r(x, colors, output_pdf):
    import rpy2.robjects as robj
    import pandas.rpy.common as com

    # save csvs for pca
    pd.DataFrame(x).T.to_csv('pca_file.csv', index=False)
    pd.Series(colors).to_csv('colors.csv', index=False)

    result = com.convert_robj(robj.r("""
    df = read.csv('pca_file.csv')

    colors = read.csv('colors.csv', header=FALSE)

    df.pca <- prcomp(df,
                     center = TRUE,
                     scale. = TRUE)
    return(df.pca)
    """))
    x = result['x']
    variance = result['sdev']

    # plot PC1 vs PC2
    fig, axis = plt.subplots(nrows=1, ncols=2)
    fig.set_figheight(10)
    fig.set_figwidth(25)

    # 1vs2 components
    for i in range(1, x.shape[0] + 1):
        axis[0].scatter(
            x.loc[i, 'PC1'], x.loc[i, 'PC2'],
            color=colors[i - 1],
            s=50
        )
    axis[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
    axis[0].set_ylabel("PC2 - {0}% variance".format(variance[1]))

    # plot PC1 vs PC3
    for i in range(1, x.shape[0] + 1):
        axis[1].scatter(
            x.loc[i, 'PC1'], x.loc[i, 'PC3'],
            color=colors[i - 1],
            s=50
        )
    axis[1].set_xlabel("PC1 - {0}% variance".format(variance[0]))
    axis[1].set_ylabel("PC3 - {0}% variance".format(variance[2]))

    fig.savefig(output_pdf, bbox_inches='tight')


def lola(bed_files, universe_file, output_folder):
    """
    Performs location overlap analysis (LOLA) on bedfiles with regions sets.
    """
    import rpy2.robjects as robj

    run = robj.r("""
        function(bedFiles, universeFile, outputFolder) {
            library("LOLA")

            userUniverse  <- LOLA::readBed(universeFile)

            dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/"
            dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/hg19/"
            regionDB = loadRegionDB(c(dbPath1, dbPath2))

            if (typeof(bedFiles) == "character") {
                userSet <- LOLA::readBed(bedFiles)
                lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                lolaResults[order(support, decreasing=TRUE), ]
                writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
            } else if (typeof(bedFiles) == "double") {
                for (bedFile in bedFiles) {
                    userSet <- LOLA::readBed(bedFile)
                    lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                    lolaResults[order(support, decreasing=TRUE), ]
                    writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
                }
            }
        }
    """)

    # convert the pandas dataframe to an R dataframe
    run(bed_files, universe_file, output_folder)


def seq2pathway(bed_file, go_term_mapping):
    """
    Performs seq2pathway analysis on a bedfile with a region set.
    """
    import rpy2.robjects as robj
    import pandas.rpy.common as com

    run = robj.r("""
        function(bed_file) {
            library("seq2pathway")

            input_df = read.table(bed_file)

            results = runseq2pathway(
                input_df, search_radius=150000, promoter_radius=2000,
                genome=c("hg19"),
                adjacent=FALSE, SNP=FALSE, PromoterStop=TRUE, NearestTwoDirection=TRUE,
                DataBase=c("GOterm", "MsigDB_C5", 'MsigDB_C6'), FAIMETest=FALSE, FisherTest=TRUE,
                collapsemethod=c("MaxMean", "function", "ME", "maxRowVariance", "MinMean", "absMinMean", "absMaxMean", "Average"),
                alpha=5, B=100, na.rm=F, min_Intersect_Count=5
            )
            # extract content
            BP = do.call(rbind.data.frame, results$gene2pathway_result.FET$GO_BP)
            MF = do.call(rbind.data.frame, results$gene2pathway_result.FET$GO_MF)
            CC = do.call(rbind.data.frame, results$gene2pathway_result.FET$GO_CC)

            return(results)
        }
    """)
    # convert to Python objects
    results = com.convert_robj(run(bed_file))
    # get only pathway results
    results = results['gene2pathway_result.FET']
    # concate results from different  ontologies
    df2 = pd.DataFrame()
    for ontology in ['MF', 'CC', 'BP']:
        df = results["GO_" + ontology]
        df['ontology'] = ontology
        df2 = pd.concat([df2, df], ignore_index=True)

    # intersect with GO term ID and name
    names = pd.read_csv(go_term_mapping)
    names.columns = ["Name", "GOID"]

    df2 = df2.merge(names)

    return df2


def goverlap(genes_file, universe_file, output_file):
    """
    from https://github.com/endrebak/goverlap
    """
    cmd = """goverlap -a {0} -s hsap -n 12 -l 0.10 -x {1} > {2}
    """.format(genes_file, universe_file, output_file)

    try:
        os.system(cmd)
    except:
        pass


def meme_ame(input_fasta, output_dir, background_fasta=None):
    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 10 -f {0} > {1}
        """.format(input_fasta, shuffled)
        os.system(cmd)

    cmd = """
    ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
    --control {0} -o {1} {2} ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv9.meme
    """.format(background_fasta if background_fasta is not None else shuffled, output_dir, input_fasta)
    os.system(cmd)

    os.system("rm %s" % shuffled)


def parse_ame(ame_dir):

    with open(os.path.join(ame_dir, "ame.txt"), 'r') as handle:
        lines = handle.readlines()

    output = list()
    for line in lines:
        # skip header lines
        if line[0] not in [str(i) for i in range(10)]:
            continue

        # get motif string and the first half of it (simple name)
        motif = line.strip().split(" ")[5].split("_")[0]
        # get corrected p-value
        q_value = line.strip().split(" ")[-2]
        # append
        output.append((motif, q_value))

    return output


def characterize_regions_structure(df, prefix, universe_df=None, plots_dir="results/plots"):
    # use all cll sites as universe
    if universe_df is None:
        universe_df = pd.read_csv(os.path.join("data", "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")

    # compare amplitude and support
    fig, axis = plt.subplots(2, 1, figsize=(5, 10))
    for i, var in enumerate(['fold_change', 'support']):
        sns.distplot(df[[var]], ax=axis[i])
        sns.distplot(universe_df[[var]], ax=axis[i])
        axis[i].set_title(var)
    fig.savefig(os.path.join(plots_dir, "%s_regions.amplitude_support.svg" % prefix), bbox_inches="tight")

    # compare 4 measurements of variability
    fig, axis = plt.subplots(2, 2, figsize=(10, 10))
    axis = axis.flatten()
    for i, var in enumerate(['variance', 'std_deviation', 'dispersion', 'qv2']):
        sns.distplot(df[[var]], ax=axis[i])
        sns.distplot(universe_df[[var]], ax=axis[i])
        axis[i].set_title(var)
    fig.savefig(os.path.join(plots_dir, "%s_regions.variability.svg" % prefix), bbox_inches="tight")

    # compare genomic regions and chromatin_states
    for i, var in enumerate(['genomic_region', 'chromatin_state']):
        # prepare:
        # separate comma-delimited fields:
        df_count = Counter(df[var].str.split(',').apply(pd.Series).stack().tolist())
        df_universe_count = Counter(universe_df[var].str.split(',').apply(pd.Series).stack().tolist())

        # divide by total:
        df_count = {k: v / float(len(df)) for k, v in df_count.items()}
        df_universe_count = {k: v / float(len(universe_df)) for k, v in df_universe_count.items()}

        # join data, sort by subset data
        both = pd.DataFrame([df_count, df_universe_count], index=['subset', 'all']).T
        both = both.sort("subset")
        both['region'] = both.index
        data = pd.melt(both, var_name="set", id_vars=['region']).replace(np.nan, 0)

        # sort for same order
        data.sort('region', inplace=True)

        g = sns.FacetGrid(col="region", data=data, col_wrap=3, sharey=True)
        g.map(sns.barplot, "set", "value")
        plt.savefig(os.path.join(plots_dir, "%s_regions.%s.svg" % (prefix, var)), bbox_inches="tight")

        fc = pd.DataFrame(np.log2(both['subset'] / both['all']), columns=['value'])
        fc['variable'] = var

        if i == 0:
            df2 = fc
        else:
            df2 = df2.append(fc)

    return df2


def characterize_regions_function(df, output_dir, prefix, data_dir="data", universe_file=None):
    # use all cll sites as universe
    if universe_file is None:
        universe_file = os.path.join(data_dir, "cll_peaks.bed")
    # get go term mapping
    go_term_mapping = os.path.join(data_dir, "goID_goName.csv")

    # make output dirs
    lola_output = os.path.join(output_dir, "lola")
    meme_output = os.path.join(output_dir, "meme")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(lola_output):
        os.makedirs(lola_output)
    if not os.path.exists(meme_output):
        os.makedirs(meme_output)

    # save to bed
    bed_file = os.path.join(output_dir, "%s_regions.bed" % prefix)
    df[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)

    tsv_file = os.path.join(output_dir, "%s_regions.tsv" % prefix)
    df['index'] = df.index
    df[['index', 'chrom', 'start', 'end']].to_csv(tsv_file, sep="\t", header=False, index=False)

    # export ensembl gene names
    gene_annotation = pd.read_csv(os.path.join('data', "cll_peaks.gene_annotation.csv"))
    ensembl_genes = pd.merge(df, gene_annotation, on=['chrom', 'start', 'end'])['ensembl_gene_id']
    ensembl_genes.to_csv(os.path.join(output_dir, "%s_genes.txt" % prefix), index=False)

    # Lola
    try:
        lola(bed_file, universe_file, lola_output)
    except:
        print("LOLA analysis for %s failed!" % prefix)

    # seq2pathway
    try:
        results = seq2pathway(tsv_file, go_term_mapping)
        results_file = os.path.join(output_dir, "%s_regions.seq2pathway.csv" % prefix)
        results.to_csv(results_file, index=False)
    except:
        print("seq2pathway analysis for %s failed!" % prefix)

    # GO Terms
    # write all gene names to file
    universe_genes_file = os.path.join(output_dir, "%s_regions.universe.closest_genes.txt" % prefix)
    with open(universe_genes_file, 'w') as handle:
        for gene in df['gene_name']:
            handle.write(gene + "\n")
    # write gene names to file
    genes_file = os.path.join(output_dir, "%s_regions.closest_genes.txt" % prefix)
    with open(genes_file, 'w') as handle:
        for gene in df['gene_name']:
            handle.write(gene + "\n")
    output_file = os.path.join(output_dir, "%s_regions.goverlap.tsv" % prefix)
    # test enrichements of closest gene function: GO, KEGG, OMIM
    try:
        goverlap(genes_file, universe_genes_file, output_file)
    except:
        print("Goverlap analysis for %s failed!" % prefix)

    # Motifs
    # de novo motif finding - enrichment
    fasta_file = os.path.join(output_dir, "%s_regions.fa" % prefix)
    bed_to_fasta(bed_file, fasta_file)

    # meme_ame(fasta_file, meme_output)


def classify_samples(analysis, sel_samples, labels, trait, rerun=False):
    """
    Use a machine learning approach for sample classification based on known sample attributes.
    Extract features most important to separate samples and investigate those.
    """
    print("Trait:%s" % trait)
    print("%i samples with trait annotated" % len(sel_samples))
    print(Counter(labels))

    dataframe_file = os.path.join(
        analysis.data_dir,
        "trait_specific",
        "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
    if os.path.exists(dataframe_file) and not rerun:  # Load up
        dataframe = pd.read_csv(dataframe_file, sep="\t")
    else:  # Run analysis
        # Get all CLL ATAC-seq samples for validation
        all_samples = [s for s in analysis.samples if s.cell_line == "CLL" and s.library == "ATAC-seq"]

        # Get colors depending on this feature label (trait) (True == green; False == redish)
        palette = sns.color_palette("colorblind")
        trait_colors = [palette[1] if l else palette[2] for l in labels]
        all_samples_colors = [trait_colors[sel_samples.index(s)] if s in sel_samples else "gray" for s in all_samples]
        cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)  # for later usage in heatmaps

        # BINARY CLASSIFICATION ON ALL CLL OPEN CHROMATIN REGIONS
        # get features and labels
        X = normalize(analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]].T)
        y = np.array(labels)

        loo = cross_validation.LeaveOneOut(len(X))

        for i, (train_index, test_index) in enumerate(loo):
            # Remove samples from the same patient that we are predicting during training
            # get current patient_id
            _pid = [s.patient_id for s in sel_samples][test_index]
            # get indexes of samples from current patient
            _p_indexes = [index for index, s in enumerate(sel_samples) if s.patient_id == _pid]
            # remove indexes from training
            train_index = np.delete(train_index, _p_indexes)
            train_index = np.delete(train_index, _p_indexes)

            # Slice accordingly
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            # Train, predict
            classifier = RandomForestClassifier(n_estimators=100, n_jobs=-1)
            y_score = classifier.fit(X_train, y_train).predict_proba(X_test)

            if i == 0:
                y_all_test = y_test
                y_all_scores = y_score
                importance = classifier.feature_importances_
            else:
                y_all_test = np.vstack([y_all_test, y_test])
                y_all_scores = np.vstack([y_all_scores, y_score])
                importance = np.vstack([importance, classifier.feature_importances_])

        # Metrics
        binary_labels = [0 if x == classifier.classes_[0] else 1 for x in y_all_test]
        binary_scores = [0 if x > 0.5 else 1 for x in y_all_scores[:, 0]]
        # Specificity (TN / N)
        tn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 1) and (binary_scores[i] == 1)])
        n = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 1])
        TNR = tn / float(n)
        # Sensitivity (TP / P)
        tp = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 0) and (binary_scores[i] == 0)])
        p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 0])
        TPR = tp / float(p)
        # FPR (FP / P)
        fn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 0) and (binary_scores[i] == 1)])
        p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 0])
        FPR = fn / float(p)
        # FNR (FN / P)
        fn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 1) and (binary_scores[i] == 0)])
        p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 1])
        FNR = fn / float(p)

        # Compute ROC curve and ROC area for each class
        fpr, tpr, _ = roc_curve(y_all_test, y_all_scores[:, 1], pos_label=1)
        roc_auc = auc(fpr, tpr, reorder=True)
        # Compute Precision-Recall and average precision
        precision, recall, _ = precision_recall_curve(y_all_test, y_all_scores[:, 1], pos_label=1)
        binary_labels = [0 if x == classifier.classes_[0] else 1 for x in y_all_test]
        aps = average_precision_score(binary_labels, y_all_scores[:, 1])

        # Plot ROC and PRC curves
        fig, axis = plt.subplots(1, 2, figsize=(12, 5))
        axis[0].plot(fpr, tpr, label='ROC (AUC = {0:0.2f}; TNR = {1:0.2f}; TPR = {2:0.2f})'.format(roc_auc, TNR, TPR))
        axis[1].plot(recall, precision, label='PRC (AUC = {0:0.2f})'.format(aps))
        axis[0].plot((0, 1), (0, 1), '--', color='gray')
        axis[0].set_xlim([-0.05, 1.0])
        axis[0].set_ylim([0.0, 1.05])
        axis[0].set_xlabel('False Positive Rate')
        axis[0].set_ylabel('True Positive Rate')
        axis[0].legend(loc="lower right")
        axis[1].set_xlim([-0.05, 1.0])
        axis[1].set_ylim([0.0, 1.05])
        axis[1].set_xlabel('Recall')
        axis[1].set_ylabel('Precision')
        axis[1].legend(loc="lower right")
        # plot specificity (tpr) and sensitivity (1-tnr)
        axis[0].plot(1 - TNR, TPR, 'o', color='gray')  # , s=50)
        fig.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.ROC_PRC.svg" % trait), bbox_inches="tight")

        # Display training and prediction of pre-labeled samples of most informative features:
        # average feature importance across iterations
        mean_importance = importance.mean(axis=0)

        # visualize feature importance
        fig, axis = plt.subplots(1)
        sns.distplot(mean_importance, ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.mean_importance.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # get important features
        # n = 500; x = matrix.loc[np.argsort(mean_importance)[-n:], :] # get n top features
        # or
        # Get most informative features
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]]
        x = matrix.loc[[i for i, j in enumerate(mean_importance > 1e-4) if j == True], :]  # get features on the tail of the importance distribution
        sites_cluster = sns.clustermap(
            x,
            cmap=cmap,
            standard_scale=0,
            col_colors=trait_colors,
            yticklabels=False)
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # SEE ALL SAMPLES
        # Display most informative features for ALL samples:
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in all_samples]]
        # get important features
        x = matrix.loc[[i for i, j in enumerate(mean_importance > 1e-4) if j == True], :]  # get features on the tail of the importance distribution

        # get colors for each cluster
        group_number = 4 if trait == "IGHV" else 2  # 4 IGHV groups, everything else 2
        # cluster all samples first
        samples_cluster = sns.clustermap(x.corr())
        # get cluster labels for samples
        Z = samples_cluster.dendrogram_col.linkage
        clusters = fcluster(Z, group_number, criterion="maxclust")
        # get cluster colors
        cluster_colors = dict(zip(np.unique(clusters), sns.color_palette("colorblind")))
        colors = [cluster_colors[c] for c in clusters]

        # sample correlation dendrogram
        sns.clustermap(
            x.corr(),
            col_colors=all_samples_colors,  # all_sample_colors(all_samples),
            row_colors=all_samples_colors)
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sample_correlation.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # pca on these regions
        pca_r(x, all_samples_colors, os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.pca.sample_labels.svg" % trait))

        # REGIONS
        # Split in major groups
        if trait == "gender":
            gr = 2
        elif trait == "IGHV":
            gr = 4
        else:
            gr = 5
        Z = sites_cluster.dendrogram_row.linkage
        clusters = fcluster(Z, gr, criterion="maxclust")
        # # visualize  cluster site attribution
        # # get cluster colors
        cluster_colors = dict(zip(np.unique(clusters), sns.color_palette("colorblind")))
        colors = [cluster_colors[c] for c in clusters]
        sns.clustermap(
            x,
            cmap=cmap,
            standard_scale=0,
            col_colors=all_samples_colors,
            row_colors=colors,
            yticklabels=False)
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sites_labels.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # Export info
        dataframe = analysis.coverage_qnorm_annotated.loc[x.index, :]
        # add difference of standardized openness between positive and negative groups used in classification
        df2 = dataframe[[s.name for s in sel_samples]].apply(lambda j: (j - j.min()) / (j.max() - j.min()), axis=0)
        dataframe["change"] = df2.icol([i for i, l in enumerate(labels) if l == 1]).mean(axis=1) - df2.icol([i for i, l in enumerate(labels) if l == 0]).mean(axis=1)
        # add feature importance
        dataframe["importance"] = mean_importance[x.index]
        # get cluster labels for sites
        dataframe["cluster"] = clusters
        # Save whole dataframe as csv
        dataframe_file = os.path.join(
            analysis.data_dir,
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
        dataframe.to_csv(dataframe_file, sep="\t", index=False)

        # Save as bed
        bed_file = os.path.join(
            analysis.data_dir,
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.sites.bed" % trait)
        dataframe[["chrom", "start", "end"]].to_csv(bed_file, sep="\t", header=False, index=False)


def characterize_regions(analysis, traits):
    """
    Characterize structural-, functionally and in the chromatin regions trait-specific regions.
    """
    for trait in traits:
        # Load dataframe with trait-specific regions
        dataframe = pd.read_csv(
            os.path.join(
                analysis.data_dir,
                "trait_specific",
                "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait
            ), sep="\t")

        # ALL REGIONS
        # region's structure
        characterize_regions_structure(df=dataframe, prefix=trait)
        # plot chromosome distribution of regions
        chrom_count = Counter(dataframe[['chrom', 'start', 'end']].chrom)
        chrom_count = np.array(sorted(chrom_count.items(), key=lambda x: x[1], reverse=True))
        fig, axis = plt.subplots(1, figsize=(10, 5))
        sns.barplot(chrom_count[:, 0], chrom_count[:, 1].astype(int), ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.sites_location.svg" % trait), bbox_inches="tight")

        # region's function
        output_dir = os.path.join(analysis.data_dir, "trait_specific", "%s_peaks" % trait)
        characterize_regions_function(df=dataframe, output_dir=output_dir, prefix=trait)

        # CLUSTERS OF REGIONS
        # Region characterization of each cluster of regions
        for i, cluster in enumerate(np.unique(dataframe['cluster'])):
            # GET REGIONS FROM CLUSTER
            df = dataframe[dataframe['cluster'] == cluster]

            # ignore clusters with less than 20 regions
            if len(df) < 20:
                continue

            # output folder
            outdir = os.path.join("%s_cluster%i" % (trait, cluster))
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            # region's structure
            regions = characterize_regions_structure(df=df, prefix="%s_cluster%i" % (trait, cluster), universe_df=analysis.coverage_qnorm_annotated)
            regions['cluster'] = cluster

            if i == 0:
                df3 = regions
            else:
                df3 = df3.append(regions)

            # region's function
            output_dir = os.path.join(analysis.data_dir, "trait_specific", "%s_peaks_cluster%i" % (trait, cluster))
            characterize_regions_function(df=df, output_dir=output_dir, prefix="%s_cluster%i" % (trait, cluster))

            # # parse meme-ame output
            # motifs = pd.DataFrame(parse_ame(os.path.join(output_dir, "meme")), columns=['motifs', 'q_values'])
            # motifs['cluster'] = cluster
            # if i == 0:
            #     df2 = motifs
            # else:
            #     df2 = pd.concat([df2, motifs])

        # Plot region enrichment
        df3['region'] = df3.index
        df3 = df3.sort(['region'])
        df3 = df3.replace(np.nan, 0)

        g = sns.FacetGrid(col="region", data=df3[df3['variable'] == 'genomic_region'], col_wrap=3, sharey=True)
        g.map(sns.barplot, "cluster", "value")
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "%s_regions.region_enrichment.svg" % trait), bbox_inches="tight")

        g = sns.FacetGrid(col="region", data=df3[df3['variable'] == 'chromatin_state'], col_wrap=3, sharey=True)
        g.map(sns.barplot, "cluster", "value")
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "trait_specific", "%s_regions.chromatin_state_enrichment.svg" % trait), bbox_inches="tight")

        # save data
        df3.to_csv(os.path.join(
            analysis.plots_dir,
            "trait_specific", "%s_regions.region_enrichment.csv" % trait), index=False)

        # # Plot motif enrichments
        # df2['q_values'] = -np.log10(np.array(df2['q_values'], dtype=float))
        # # get highest (worst) p-value from motifs of the same TF
        # df3 = df2.groupby(['motifs', 'cluster']).aggregate(max).reset_index()
        # # spread again for each cluster
        # df3 = df3.pivot('motifs', 'cluster', 'q_values')
        # # replace nan with 0
        # df3 = df3.replace(np.nan, 0)
        # # fix types
        # for i in [1, 2]:
        #     df3[i] = df3[i].astype(float)

        # # sort
        # df3 = df3.sort([1])
        # # plot heatmap of p-values
        # sns.clustermap(df3, z_score=1, figsize=(8, 24), linewidths=0, cmap=plt.cm.YlGn)
        # plt.savefig(os.path.join(
        #     analysis.plots_dir,
        #     "trait_specific", "%s_regions.motif_enrichment.svg" % trait), bbox_inches="tight")
        # plt.close("all")

        # sns.clustermap(df3[(df3.icol(0) > 30) & (df3.icol(1) > 30)], cmap=plt.cm.YlGn)
        # plt.savefig(os.path.join(
        #     analysis.plots_dir,
        #     "trait_specific", "%s_regions.motif_enrichment.highest.svg" % trait), bbox_inches="tight")
        # plt.close("all")

        # # save data
        # df3.to_csv(os.path.join(
        #     analysis.plots_dir,
        #     "trait_specific", "%s_regions.motif_enrichment.csv" % trait), index=False)

        # # get n_max most different motifs
        # n_max = 25
        # p_value = 3  # -log10

        # # filter p-values
        # df4 = df3[(df3[1] > p_value) | (df3[2] > p_value)]

        # # floor n_max
        # if len(df4) < n_max:
        #     n_max = len(df4)

        # # get n most dissimilar features
        # s = abs(df4[1] - df4[2])
        # # sort by similarity
        # s.sort(ascending=False)
        # # get top
        # index = s.irow(range(n_max)).index

        # # plot matrix of clusters/terms with p_values
        # sns.clustermap(df4.ix[index], col_cluster=False, cmap=plt.cm.YlGn)
        # plt.savefig(os.path.join(
        #     analysis.plots_dir,
        #     "trait_specific", "%s_regions.motif_enrichment.difference.svg" % trait), bbox_inches="tight")
        # plt.close("all")


def classification_validation(analysis, train_samples, train_labels, val_samples, val_labels, comparison):
    """
    Use a machine learning approach for sample classification based on known sample attributes.
    Extract features most important to separate samples and investigate those.
    """
    print("Trait:%s" % comparison)
    print("Train:")
    print("%i samples with trait annotated" % len(train_samples))
    print(Counter(train_labels))
    print("Validation:")
    print("%i samples with trait annotated" % len(val_samples))
    print(Counter(val_labels))

    # ALL CLL OPEN CHROMATIN REGIONS
    matrix = pd.DataFrame(
        normalize(
            analysis.coverage_qnorm_annotated[[s.name for s in analysis.samples if s.library == "ATAC-seq"]]
        ),
        columns=[[s.name for s in analysis.samples if s.library == "ATAC-seq"]])
    matrix_train = matrix[[s.name for s in train_samples]]
    matrix_val = matrix[[s.name for s in val_samples]]

    # BINARY CLASSIFICATION
    # get features and train_labels
    X_train = matrix_train.T
    y_train = np.array(train_labels)
    X_val = matrix_val.T
    y_val = np.array(val_labels)

    # Train, predict
    classifier = RandomForestClassifier(n_estimators=100, n_jobs=-1)
    y_score = classifier.fit(X_train, y_train).predict_proba(X_val)

    # Metrics
    binary_labels = [0 if x == classifier.classes_[0] else 1 for x in y_val]
    binary_scores = [0 if x > 0.5 else 1 for x in y_score[:, 0]]
    # Specificity (TN / N)
    tn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 1) and (binary_scores[i] == 1)])
    n = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 1])
    TNR = tn / float(n)
    # Sensitivity (TP / P)
    tp = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 0) and (binary_scores[i] == 0)])
    p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 0])
    TPR = tp / float(p)
    # FPR (FP / P)
    fn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 0) and (binary_scores[i] == 1)])
    p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 0])
    FPR = fn / float(p)
    # FNR (FN / P)
    fn = len([1 for i in range(len(binary_scores)) if (binary_labels[i] == 1) and (binary_scores[i] == 0)])
    p = len([1 for i in range(len(binary_scores)) if binary_labels[i] == 1])
    FNR = fn / float(p)

    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(y_val, y_score[:, 1], pos_label=1)
    roc_auc = auc(fpr, tpr, reorder=True)
    # Compute Precision-Recall and average precision
    precision, recall, _ = precision_recall_curve(y_val, y_score[:, 1], pos_label=1)
    binary_labels = [0 if x == classifier.classes_[0] else 1 for x in y_val]
    aps = average_precision_score(binary_labels, y_score[:, 1])

    # Plot ROC and PRC curves
    fig, axis = plt.subplots(1, 2, figsize=(12, 5))
    axis[0].plot(fpr, tpr, label='ROC (AUC = {0:0.2f}; TNR = {1:0.2f}; TPR = {2:0.2f})'.format(roc_auc, TNR, TPR))
    axis[1].plot(recall, precision, label='PRC (AUC = {0:0.2f})'.format(aps))
    axis[0].plot((0, 1), (0, 1), '--', color='gray')
    axis[0].set_xlim([-0.05, 1.0])
    axis[0].set_ylim([0.0, 1.05])
    axis[0].set_xlabel('False Positive Rate')
    axis[0].set_ylabel('True Positive Rate')
    axis[0].legend(loc="lower right")
    axis[1].set_xlim([-0.05, 1.0])
    axis[1].set_ylim([0.0, 1.05])
    axis[1].set_xlabel('Recall')
    axis[1].set_ylabel('Precision')
    axis[1].legend(loc="lower right")
    # plot specificity (tpr) and sensitivity (1-tnr)
    axis[0].plot(1 - TNR, TPR, 'o', color='gray')  # , s=50)
    fig.savefig(os.path.join(
        analysis.plots_dir,
        "trait_specific", "cll_peaks.%s_significant.classification_validation.random_forest.loocv.ROC_PRC.svg" % comparison), bbox_inches="tight")


def trait_analysis(analysis, samples, traits):
    """
    Run trait classification (with independent validation if possible for that trait)
    on all samples with known annotation of the trait.
    """
    for trait in traits:
        if trait != "ibrutinib":
            # cross-validated classification
            sel_samples = [s for s in samples if getattr(s, trait) is not pd.np.nan]
            labels = np.array([getattr(s, trait) for s in sel_samples])
            classify_samples(analysis, sel_samples, labels, trait, rerun=True)
            # classification with independent validation
            if trait in ["gender", "IGHV", "CD38", "treated"]:
                train_samples = [s for s in samples if getattr(s, trait) is not pd.np.nan and s.hospital != "AKH"]
                train_labels = np.array([getattr(s, trait) for s in train_samples])
                val_samples = [s for s in samples if getattr(s, trait) is not pd.np.nan and s.hospital == "AKH"]
                val_labels = np.array([getattr(s, trait) for s in val_samples])
                classification_validation(analysis, train_samples, train_labels, val_samples, val_labels, trait)
        # for ibrutinib, train only on AKH samples
        else:
            sel_samples = [s for s in samples if s.hospital == "AKH"]
            labels = np.array([1 if s.under_treatment else 0 for s in sel_samples])
            classify_samples(analysis, sel_samples, labels, trait, rerun=True)


def join_trait_specific_regions(analysis, traits):
    # Put trait-specific chromatin regions in one matrix
    features = pd.DataFrame()
    for trait in traits:
        file_name = os.path.join(
            "data",
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
        try:
            df = pd.read_csv(file_name, sep="\t")
            df['trait'] = trait
            features = features.append(df, ignore_index=True)
        except IOError:
            print("Trait %s did not generate any associated regions" % trait)

    # add direction of chromatin feature association with clinical feature
    features['direction'] = features['change'].apply(lambda x: 1 if x > 0 else -1)

    # # Save whole dataframe as csv
    features.to_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.csv"), sep="\t", index=False)


def characterize_regions_chromatin(analysis, traits):
    """
    For each trait-associated region, get ratios of active/repressed and poised/repressed chromatin,
    across all patients or groups of patients with same IGHV mutation status.

    Visualize histone mark ratios dependent on the positive- or negative-association of regions with
    accessibility signal.
    """
    def get_intensities(df, subset, mark):
        samples = [s for s in analysis.samples if s.library == "ChIPmentation"]
        # subset samples according to requested group (all, uCLL, iCLL or mCLL)
        if subset == "all":
            samples = [s for s in samples if s.ighv_group is not pd.np.nan]
        else:
            samples = [s for s in samples if s.ighv_group == subset]

        # mean of that mark across the subset of samples
        return df[[s.name for s in samples if s.ip == mark]].mean(axis=1)

    def calculate_ratio(df, subset, kind):
        samples = [s for s in analysis.samples if s.library == "ChIPmentation"]
        # subset samples according to requested group (all, uCLL, iCLL or mCLL)
        if subset == "all":
            samples = [s for s in samples if s.ighv_group is not pd.np.nan]
        else:
            samples = [s for s in samples if s.ighv_group == subset]

        # get ratio of the requested histone marks
        if kind == "A:R":
            a = df[[s.name for s in samples if s.ip == "H3K27ac"]]
            r = df[[s.name for s in samples if s.ip == "H3K27me3"]]
        elif kind == "P:R":
            a = df[[s.name for s in samples if s.ip == "H3K4me1"]]
            r = df[[s.name for s in samples if s.ip == "H3K27me3"]]
        return np.log2(((1 + a.values) / (1 + r.values))).mean(axis=1)

    # read in dataframe
    features = pd.read_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.csv"), sep="\t")

    # get intensity of chromatin marks
    # across all patients
    for group in ["all", "uCLL", "iCLL", "mCLL"]:
        for mark in ["H3K27ac", "H3K4me1", "H3K27me3"]:
            features["%s_intensity_%s" % (group, mark)] = get_intensities(features, group, mark)

    # compute ratios of chromatin marks
    # across all patients and for each each IGHV class
    for group in ["all", "uCLL", "iCLL", "mCLL"]:
        for ratio in ["A:R", "P:R"]:
            features["%s_ratio_%s" % (group, ratio)] = calculate_ratio(features, group, ratio)

    # save dataframe with intensities/ratios of chromatin marks per peak
    features.to_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.histone_intensities_ratios.csv"), index=False)

    # Heatmap accessibility and histone marks
    # for each trait make heatmap with chroamtin marks in each
    for trait in traits:
        # get regions specific to this trait
        df = features[features["trait"] == trait]

        # get all ATAC-seq samples
        samples = [s for s in analysis.samples if s.library == "ATAC-seq"]

        # color samples according to mutation status
        palette = sns.color_palette("colorblind")
        trait_colors = {True: palette[1], False: palette[2]}
        sample_colors = [trait_colors[s.igvh_mutation_status] if s.igvh_mutation_status is not pd.np.nan else "gray" for s in samples]

        # colormap for heatmap values
        cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)

        # heatmap
        dend = sns.clustermap(
            df[[s.name for s in samples]],
            cmap=cmap,
            standard_scale=0,
            col_colors=sample_colors,
            yticklabels=False)
        plt.savefig("/home/arendeiro/ighv_atac.svg", bbox_inches="tight")

        # order rows by ATAC-seq dendrogram
        chrom = df[df.columns[df.columns.str.contains("intensity")]]
        order = dend.dendrogram_row.dendrogram['leaves']
        chrom = chrom.ix[order]

        # order columns by histone mark
        chrom = chrom[sorted(chrom.columns.tolist(), key=lambda x: x.split("_")[2])]
        sns.heatmap(chrom, yticklabels=False)
        plt.savefig("/home/arendeiro/ighv_histones.ordered_mark.svg", bbox_inches="tight")
        plt.close("all")

        # order columns by group
        chrom = chrom[sorted(chrom.columns.tolist(), key=lambda x: (x.split("_")[0], x.split("_")[2]))]
        sns.heatmap(chrom, yticklabels=False)
        plt.savefig("/home/arendeiro/ighv_histones.ordered_group.svg", bbox_inches="tight")
        plt.close("all")

        # Z-scored, order columns by histone mark
        chrom = chrom[sorted(chrom.columns.tolist(), key=lambda x: x.split("_")[2])]
        sns.heatmap(chrom.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig("/home/arendeiro/ighv_histones.ordered_mark.zscore_rows.svg", bbox_inches="tight")
        plt.close("all")

        # Z-scored, order columns by group
        chrom = chrom[sorted(chrom.columns.tolist(), key=lambda x: (x.split("_")[0], x.split("_")[2]))]
        sns.heatmap(chrom.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig("/home/arendeiro/ighv_histones.ordered_group.zscore_rows.svg", bbox_inches="tight")
        plt.close("all")

    # Plot values
    # subset data and melt dataframe for ploting
    features_slim = pd.melt(
        features[
            ["trait", "direction"] +
            ["%s_intensity_H3K27ac" % group for group in ["uCLL", "iCLL", "mCLL"]] +
            ["%s_intensity_H3K4me1" % group for group in ["uCLL", "iCLL", "mCLL"]] +
            ["%s_intensity_H3K27me3" % group for group in ["uCLL", "iCLL", "mCLL"]]
        ],
        id_vars=["trait", "direction"],
        var_name="mark_type",
        value_name="intensity"
    )
    features_slim["group"] = features_slim["mark_type"].apply(lambda x: x.split("_")[0])
    features_slim["mark"] = features_slim["mark_type"].apply(lambda x: x.split("_")[2])

    # traits as columns, group by positively- and negatively-associated regions
    for trait in ["IGHV", "CD38", "ZAP70"]:
        p = features_slim[features_slim['trait'] == trait].dropna()

        g = sns.FacetGrid(p, col="group", row="direction", legend_out=True, row_order=[-1, 1], margin_titles=True)
        g.map(sns.violinplot, "mark", "intensity", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_intensity.%s.mark_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(p, col="group", row="mark", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "intensity", split=True, order=[-1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_intensity.%s.direction_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(p, row="mark", col="direction", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "group", "intensity", split=True, order=["uCLL", "iCLL", "mCLL"])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_intensity.%s.group_centric.svg" % trait),
            bbox_inches='tight')

    # Plot ratios
    # subset data and melt dataframe for ploting
    features_slim = pd.melt(
        features[
            ["trait", "direction"] +
            ["%s_ratio_A:R" % group for group in ["all", "uCLL", "iCLL", "mCLL"]] +
            ["%s_ratio_P:R" % group for group in ["all", "uCLL", "iCLL", "mCLL"]]
        ],
        id_vars=["trait", "direction"],
        var_name="ratio_type",
        value_name="ratio"
    )
    features_slim["group"] = features_slim["ratio_type"].apply(lambda x: x.split("_")[0])
    features_slim["type"] = features_slim["ratio_type"].apply(lambda x: x.split("_")[2])

    # traits as columns, group by positively- and negatively-associated regions
    for trait in ["IGHV", "CD38", "ZAP70"]:
        p = features_slim[features_slim['trait'] == trait].dropna()

        g = sns.FacetGrid(p, col="group", row="direction", legend_out=True)
        g.map(sns.violinplot, "type", "ratio", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_ratios.%s.mark_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(p, col="group", row="type", legend_out=True)
        g.map(sns.violinplot, "direction", "ratio", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_ratios.%s.direction_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(p, row="type", col="direction", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "group", "ratio", split=True, order=["uCLL", "iCLL", "mCLL"])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "trait_specific", "cll.trait-specific_regions.chromatin_ratios.%s.group_centric.svg" % trait),
            bbox_inches='tight')


def generate_signature_matrix(array, n=101, bounds=(0, 0)):
    """
    :param np.array: 2D np.array
    """
    def get_score(i, j, p, n):
        """Get signature score between p% of the values of the two groups."""
        return ((float(i) * p) + (float(j) * (n - p))) / n

    matrix = np.zeros([array.shape[0], n])
    for x in range(array.shape[0]):
        for y, p in enumerate(np.linspace(0 + bounds[0], n + bounds[1], n)):
            matrix[x, y] = get_score(array[x, 0], array[x, 1], p, n)

    return matrix


def best_signature_matrix(array, matrix):
    """
    :param np.array: 2D np.array
    """
    from scipy.stats import pearsonr

    cors = dict()
    for i in range(matrix.shape[1]):
        cors[i] = pearsonr(array, matrix[:, i])

    return cors.values().index(max(cors.values()))  # index
    # (
    # cors.values().index(max(cors.values())),  # index
    # max(cors.values())  # highest correlation value


def pairwise(iterable):
            from itertools import tee, izip
            "s -> (s0,s1), (s1,s2), (s2, s3), ..."
            a, b = tee(iterable)
            next(b, None)
            return izip(a, b)


def get_signatures(analysis, traits):
    """
    Assign samples to a trait-related signature
    """
    def biplot(x, y, samples):
        # Z-score variables to get them into same space
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)
        yy = y.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        # Plot samples and variables in same space (biplot)
        fig, axis = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)

        sample_colors = samples_to_color(samples, "unique")
        sample_symbols = samples_to_symbol(samples, "unique")

        var_colors = sns.color_palette("Paired", yy.shape[0])

        # Plot observations (samples)
        samples = pd.DataFrame([pd.Series(sample.__dict__) for sample in samples])
        for i in range(len(xx)):
            axis.scatter(xx.ix[i][0], xx.ix[i][1], s=50, color=sample_colors[i], marker=sample_symbols[i], label=x.index[i])

        # Plot variables (clinical features)
        for i, trait in enumerate(yy.index):
            axis.plot((0, yy.ix[trait][0]), (0, yy.ix[trait][1]), '-o', color=var_colors[i], label=x.columns[i])

        # # add dashed line between patient's timepoints
        # for patient, indexes in samples.groupby("patient_id").groups.items():
        #     for t1, t2 in pairwise(sorted(samples.ix[indexes]["timepoint"])):
        #         tt1 = samples[(samples["timepoint"] == sorted([t1, t2])[0]) & (samples["patient_id"] == patient)].index
        #         tt2 = samples[(samples["timepoint"] == sorted([t1, t2])[1]) & (samples["patient_id"] == patient)].index
        #         if len(tt1) == 1 and len(tt2) == 1:  # this is a problem due to one patient has two timepoints with same number :grrrr:
        #             axis.annotate(
        #                 "",  # samples.ix[t1]["patient_id"],
        #                 xy=(xx.ix[tt1][0].squeeze(), xx.ix[tt1][1].squeeze()), xycoords="data",
        #                 xytext=(xx.ix[tt2][0].squeeze(), xx.ix[tt2][1].squeeze()), textcoords="data",
        #                 arrowprops=dict(arrowstyle="fancy", color="0.5", shrinkB=5, connectionstyle="arc3,rad=0.3",),
        #             )
        axis.legend()

    # Read in openness values in regions associated with clinical traits
    features = pd.read_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.csv"), sep="\t")

    # Position each patient within the trait-specific chromatin signature
    samples = [s for s in analysis.samples if s.cell_line == "CLL" and s.library == "ATAC-seq"]
    muts = ['del11', 'tri12', "del17", 'TP53']
    sigs = pd.DataFrame(index=[s.name for s in samples])
    for i, trait in enumerate(traits):
        print(trait)

        # Get trait-specific signature
        # 1. get median accessibility of each group
        x = features[features["trait"] == trait].drop_duplicates(['chrom', 'start', 'end'])
        if trait not in muts:
            x1 = x[[s.name for s in samples if getattr(s, trait) is 1]].apply(np.median, axis=1)
            x2 = x[[s.name for s in samples if getattr(s, trait) is 0]].apply(np.median, axis=1)
        else:
            x1 = x[[s.name for s in samples if getattr(s, trait)]].apply(np.median, axis=1)
            x2 = x[[s.name for s in samples if not getattr(s, trait)]].apply(np.median, axis=1)

        # 2. get signature matrix
        # here, bounds are set to (-20, 20) so that the extremes of the signature represent -20% or 120% of the signature
        # this is done because the extreme values (0 and 100%) values correspond to the median value within each group,
        # meaning that some samples are expected to surpass those values.
        sign = generate_signature_matrix(np.vstack([x1, x2]).T, n=101, bounds=(-20, 20))

        # 3. get signature value of each patient
        # x[[s.name for s in samples]].apply(best_signature_matrix, matrix=sign, axis=1)
        trait_sigs = list()
        for name in [s.name for s in samples]:
            if name in x.columns:  # this is temporary
                trait_sigs.append(best_signature_matrix(array=x[name], matrix=sign))
        sigs[trait] = trait_sigs

    # Plot distribution of signature values
    # sns.distplot(trait_sigs, bins=100, kde=False, ax=axis[i])
    g = sns.FacetGrid(data=pd.melt(sigs), col="variable", col_wrap=4, sharey=False)
    g.map(sns.distplot, "value")
    # Save fig
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.trait_signatures.101.svg")
    plt.savefig(output_pdf, bbox_inches='tight')

    # correlate signatures
    sns.clustermap(sigs.corr())
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.trait_signatures.sigs_correlation.svg")
    plt.savefig(output_pdf, bbox_inches='tight')

    # correlate samples
    sns.clustermap(sigs.T.corr())
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.trait_signatures.samples_correlation.svg")
    plt.savefig(output_pdf, bbox_inches='tight')

    # Plots of signature enrichment vs signature rank
    fig, axis = plt.subplots(nrows=2, ncols=5, sharex=True, sharey=True, figsize=(20, 8))
    axis = axis.flatten()
    for i, trait in enumerate(traits):
        d = sigs[trait].copy()
        d.sort()

        sns.distplot(
            d, bins=50, vertical=True,
            rug=True, hist=False, kde=False,
            ax=axis[i])
        axis[i].set_xlim((0, .05))
        axis[i].set_ylim((-5, 105))

        for tl in axis[i].get_xticklabels():
            tl.set_color('blue')

        ax2 = axis[i].twiny()
        ax2.plot(
            d.rank(method='first'),
            d,
            "-o",
            color=sns.color_palette("colorblind")[1],
            label="rank"
        )
        ax2.set_ylim((-5, 105))
        ax2.set_xlim((-10, len(d) + 5))
        # axis[i].set_title(trait)
        for tl in ax2.get_xticklabels():
            tl.set_color('green')

    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.trait_signatures.rank_distribution.svg")
    fig.savefig(output_pdf, bbox_inches='tight')

    # # PCA on signatures (experimental)
    # from sklearn.decomposition import PCA
    # pca = PCA()
    # x_new = pca.fit_transform(sigs.T)

    # # Variable space
    # loadings = pd.DataFrame(pca.components_, columns=sigs.index, index=sigs.columns).T

    # # sum vectors for each clinical variable, dependent on the direction of the features (positive or negatively associated)
    # biplot(x=loadings, y=pd.DataFrame(x_new), samples=analysis.samples)


class gaussian_kde(object):
    """Representation of a kernel-density estimate using Gaussian kernels.

    Kernel density estimation is a way to estimate the probability density
    function (PDF) of a random variable in a non-parametric way.
    `gaussian_kde` works for both uni-variate and multi-variate data.   It
    includes automatic bandwidth determination.  The estimation works best for
    a unimodal distribution; bimodal or multi-modal distributions tend to be
    oversmoothed.

    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a scalar,
        this will be used directly as `kde.factor`.  If a callable, it should
        take a `gaussian_kde` instance as only parameter and return a scalar.
        If None (default), 'scott' is used.  See Notes for more details.
    weights : array_like, shape (n, ), optional, default: None
        An array of weights, of the same shape as `x`.  Each value in `x`
        only contributes its associated weight towards the bin count
        (instead of 1).

    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.
    d : int
        Number of dimensions.
    n : int
        Number of datapoints.
    neff : float
        Effective sample size using Kish's approximation.
    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.
    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).
    inv_cov : ndarray
        The inverse of `covariance`.

    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.
    kde(points) : ndarray
        Same as kde.evaluate(points)
    kde.pdf(points) : ndarray
        Alias for ``kde.evaluate(points)``.
    kde.set_bandwidth(bw_method='scott') : None
        Computes the bandwidth, i.e. the coefficient that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        .. versionadded:: 0.11.0
    kde.covariance_factor : float
        Computes the coefficient (`kde.factor`) that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        The default is `scotts_factor`.  A subclass can overwrite this method
        to provide a different method, or set it through a call to
        `kde.set_bandwidth`.

    Notes
    -----
    Bandwidth selection strongly influences the estimate obtained from the KDE
    (much more so than the actual shape of the kernel).  Bandwidth selection
    can be done by a "rule of thumb", by cross-validation, by "plug-in
    methods" or by other means; see [3]_, [4]_ for reviews.  `gaussian_kde`
    uses a rule of thumb, the default is Scott's Rule.

    Scott's Rule [1]_, implemented as `scotts_factor`, is::

        n**(-1./(d+4)),

    with ``n`` the number of data points and ``d`` the number of dimensions.
    Silverman's Rule [2]_, implemented as `silverman_factor`, is::

        (n * (d + 2) / 4.)**(-1. / (d + 4)).

    Good general descriptions of kernel density estimation can be found in [1]_
    and [2]_, the mathematics for this multi-dimensional implementation can be
    found in [1]_.

    References
    ----------
    .. [1] D.W. Scott, "Multivariate Density Estimation: Theory, Practice, and
           Visualization", John Wiley & Sons, New York, Chicester, 1992.
    .. [2] B.W. Silverman, "Density Estimation for Statistics and Data
           Analysis", Vol. 26, Monographs on Statistics and Applied Probability,
           Chapman and Hall, London, 1986.
    .. [3] B.A. Turlach, "Bandwidth Selection in Kernel Density Estimation: A
           Review", CORE and Institut de Statistique, Vol. 19, pp. 1-33, 1993.
    .. [4] D.M. Bashtannyk and R.J. Hyndman, "Bandwidth selection for kernel
           conditional density estimation", Computational Statistics & Data
           Analysis, Vol. 36, pp. 279-298, 2001.

    Examples
    --------
    Generate some random two-dimensional data:

    >>> from scipy import stats
    >>> def measure(n):
    >>>     "Measurement model, return two coupled measurements."
    >>>     m1 = np.random.normal(size=n)
    >>>     m2 = np.random.normal(scale=0.5, size=n)
    >>>     return m1+m2, m1-m2

    >>> m1, m2 = measure(2000)
    >>> xmin = m1.min()
    >>> xmax = m1.max()
    >>> ymin = m2.min()
    >>> ymax = m2.max()

    Perform a kernel density estimate on the data:

    >>> X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    >>> positions = np.vstack([X.ravel(), Y.ravel()])
    >>> values = np.vstack([m1, m2])
    >>> kernel = stats.gaussian_kde(values)
    >>> Z = np.reshape(kernel(positions).T, X.shape)

    Plot the results:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    ...           extent=[xmin, xmax, ymin, ymax])
    >>> ax.plot(m1, m2, 'k.', markersize=2)
    >>> ax.set_xlim([xmin, xmax])
    >>> ax.set_ylim([ymin, ymax])
    >>> plt.show()

    """
    def __init__(self, dataset, bw_method=None, weights=None):
        self.dataset = np.atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")
        self.d, self.n = self.dataset.shape

        if weights is not None:
            self.weights = weights / np.sum(weights)
        else:
            self.weights = np.ones(self.n) / self.n

        # Compute the effective sample size
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
        self.neff = 1.0 / np.sum(self.weights ** 2)

        self.set_bandwidth(bw_method=bw_method)

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.

        """
        from scipy.spatial.distance import cdist

        points = np.atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = np.reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d, self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis', VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights, axis=1) / self._norm_factor

        return result

    __call__ = evaluate

    def scotts_factor(self):
        return np.power(self.neff, -1. / (self.d + 4))

    def silverman_factor(self):
        return np.power(self.neff * (self.d + 2.0) / 4.0, -1. / (self.d + 4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.

        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.

        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth.  This can be
            'scott', 'silverman', a scalar constant or a callable.  If a
            scalar, this will be used directly as `kde.factor`.  If a callable,
            it should take a `gaussian_kde` instance as only parameter and
            return a scalar.  If None (default), nothing happens; the current
            `kde.covariance_factor` method is kept.

        Notes
        -----
        .. versionadded:: 0.11

        Examples
        --------
        >>> x1 = np.array([-7, -5, 1, 4, 5.])
        >>> kde = stats.gaussian_kde(x1)
        >>> xs = np.linspace(-10, 10, num=50)
        >>> y1 = kde(xs)
        >>> kde.set_bandwidth(bw_method='silverman')
        >>> y2 = kde(xs)
        >>> kde.set_bandwidth(bw_method=kde.factor / 3.)
        >>> y3 = kde(xs)

        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.plot(x1, np.ones(x1.shape) / (4. * x1.size), 'bo',
        ...         label='Data points (rescaled)')
        >>> ax.plot(xs, y1, label='Scott (default)')
        >>> ax.plot(xs, y2, label='Silverman')
        >>> ax.plot(xs, y3, label='Const (1/3 * Silverman)')
        >>> ax.legend()
        >>> plt.show()

        """
        from six import string_types
        if bw_method is None:
            pass
        elif bw_method == 'scott':
            self.covariance_factor = self.scotts_factor
        elif bw_method == 'silverman':
            self.covariance_factor = self.silverman_factor
        elif np.isscalar(bw_method) and not isinstance(bw_method, string_types):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            msg = "`bw_method` should be 'scott', 'silverman', a scalar " \
                  "or a callable."
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            # Compute the mean and residuals
            _mean = np.sum(self.weights * self.dataset, axis=1)
            _residual = (self.dataset - _mean[:, None])
            # Compute the biased covariance
            self._data_covariance = np.atleast_2d(np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor ** 2
        self.inv_cov = self._data_inv_cov / self.factor ** 2
        self._norm_factor = np.sqrt(np.linalg.det(2 * np.pi * self.covariance))  # * self.n


def create_clinical_epigenomic_space(analysis, traits):
    """
    """
    from sklearn.decomposition import PCA
    import cPickle as pickle
    import itertools

    def plot_pca(x_new):
        # plot PC1 vs PC2
        fig, axis = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)

        colors = samples_to_color(analysis.samples, "unique")
        symbols = samples_to_symbol(analysis.samples, "unique")
        names = [s.name for s in analysis.samples]

        for i in range(1, x_new.shape[0]):
            axis[0].scatter(  # plot PC1 vs PC3
                x_new[i, 0], x_new[i, 1],
                color=colors[i - 1],
                s=50,
                label=names[i - 1],
                marker=symbols[i - 1]
            )
            axis[1].scatter(  # plot PC1 vs PC3
                x_new[i, 0], x_new[i, 1],
                color=colors[i - 1],
                s=50,
                label=names[i - 1],
                marker=symbols[i - 1]
            )
        axis[0].set_xlabel("PC1 - {0}% variance".format(pca.explained_variance_[0]))
        axis[0].set_ylabel("PC2 - {0}% variance".format(pca.explained_variance_[1]))
        axis[1].set_xlabel("PC1 - {0}% variance".format(pca.explained_variance_[0]))
        axis[1].set_ylabel("PC3 - {0}% variance".format(pca.explained_variance_[2]))
        axis[1].legend()
        output_pdf = os.path.join(
            analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.svg")
        fig.savefig(output_pdf, bbox_inches='tight')

    def sum_cartesian_vectors(vectors):
        """
        Given A = (xi, yi, zi, ...) and  B = (xj, yj, zj, ...);
        sum A + B = (xi + xj, yi + yj, zi + zj, ...)
        """
        return vectors.sum(axis=0)

    def plot_space(x, y, scale=1.0, axis=None):
        standalone = True if axis is None else False

        # Z-score variables to get them into same space
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)
        yy = y.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        # Plot samples and variables in same space (biplot)
        if standalone:
            fig, axis = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)

        sample_colors = samples_to_color(analysis.samples, "unique")
        sample_symbols = samples_to_symbol(analysis.samples, "unique")

        var_colors = sns.color_palette("Paired", yy.shape[0])

        # Plot observations (samples)
        samples = pd.DataFrame([pd.Series(sample.__dict__) for sample in analysis.samples])
        for i in range(len(xx)):
            axis.scatter(xx.ix[i][0], xx.ix[i][1], s=50, color=sample_colors[i], marker=sample_symbols[i])

        # Plot variables (clinical features)
        for i, trait in enumerate(yy.index):
            axis.plot((0, yy.ix[trait][0]), (0, yy.ix[trait][1]), '-o', color=var_colors[i], label=trait)

        # add dashed line between patient's timepoints
        for patient, indexes in samples.groupby("patient_id").groups.items():
            for t1, t2 in pairwise(sorted(samples.ix[indexes]["timepoint"])):
                tt1 = samples[(samples["timepoint"] == sorted([t1, t2])[0]) & (samples["patient_id"] == patient)].index
                tt2 = samples[(samples["timepoint"] == sorted([t1, t2])[1]) & (samples["patient_id"] == patient)].index
                if len(tt1) == 1 and len(tt2) == 1:  # this is a problem due to one patient has two timepoints with same number :grrrr:
                    axis.annotate(
                        "",  # samples.ix[t1]["patient_id"],
                        xy=(xx.ix[tt1][0].squeeze(), xx.ix[tt1][1].squeeze()), xycoords="data",
                        xytext=(xx.ix[tt2][0].squeeze(), xx.ix[tt2][1].squeeze()), textcoords="data",
                        arrowprops=dict(arrowstyle="fancy", color="0.5", shrinkB=5, connectionstyle="arc3,rad=0.3",),
                    )
        axis.legend()

        if standalone:
            output_pdf = os.path.join(
                analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.biplot.all_clinical_feratures.svg")
            fig.savefig(output_pdf, bbox_inches='tight')
        else:
            return axis

    def add_survival(chrom_x, chrom_y, kind="survival"):
        # Get survival and hazard predictions
        # This comes from the "/src/survival_analysis.py" script
        survival = pd.read_csv(os.path.join("data", "survival_hazard_predictions.csv"))

        # get survival of each sample as weights
        s = pd.Series([survival[survival['sample_id'] == s.sampleID]['predicted_%s' % kind].squeeze() for s in analysis.samples])
        weights = s[[type(i) is np.float64 for i in s]].dropna()

        # remove points for patients without survival predictions
        xy = obs[[0, 1]].values
        xy = xy[weights.index, :]
        xy = np.transpose(xy)

        # kde space
        xmin, xmax = (xy[0].min() - 1, xy[0].max() + 1)
        ymin, ymax = (xy[1].min() - 1, xy[1].max() + 1)
        x = np.linspace(xmin, xmax, 100)  # kde resolution
        y = np.linspace(ymin, ymax, 100)  # kde resolution
        xx, yy = np.meshgrid(x, y)

        # Unweighted KDE
        pdf = gaussian_kde(xy)
        pdf.set_bandwidth(bw_method=pdf.factor / 1.5)  # kde bandwidth
        zz = pdf((np.ravel(xx), np.ravel(yy)))
        zz1 = np.reshape(zz, xx.shape)

        # weighted KDE
        pdf = gaussian_kde(xy, weights=10 * np.array(weights, np.float))
        pdf.set_bandwidth(bw_method=pdf.factor / 1.5)  # kde bandwidth
        zz2 = pdf((np.ravel(xx), np.ravel(yy)))
        zz2 = np.reshape(zz2, xx.shape)

        # PLot
        fig, axis = plt.subplots(2, figsize=(30, 20))
        bounds = [xmin, xmax, ymin, ymax]

        # plot unweighted
        # axis[0].scatter(obs[:, 0], obs[:, 1])
        axis[0] = plot_space(x=chrom_x, y=chrom_y, axis=axis[0])  # chrom_x[weights.index, :]
        cax = axis[0].imshow(np.rot90(zz1.T), cmap=plt.cm.Reds, extent=bounds, alpha=0.5)

        # plot weighted
        # axis[1].scatter(obs[:, 0], obs[:, 1])
        axis[1] = plot_space(x=chrom_x, y=chrom_y, axis=axis[1])  # chrom_x[weights.index, :]
        cax = axis[1].imshow(np.rot90(zz2.T), cmap=plt.cm.Reds, extent=bounds, alpha=0.5)

        # graphical stuff
        axis[0].legend()
        axis[1].legend()
        axis[0].set_title("unweighted %s" % kind)
        axis[1].set_title("weighted %s" % kind)
        cbar = fig.colorbar(cax, ticks=[zz2.min(), zz2.max()])
        cbar.ax.set_xticklabels(['Low', 'High'])  # horizontal colorbar
        output_pdf = os.path.join(
            analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.%s_weight-comparison.svg" % kind)
        fig.savefig(output_pdf, bbox_inches='tight')

        # compare with scatter values
        fig, axis = plt.subplots(1)
        plt.scatter(x_new[:, 0], x_new[:, 1], color=map(cm.Reds, 1 - weights))
        output_pdf = os.path.join(
            analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.raw_1-survival.svg")
        fig.savefig(output_pdf, bbox_inches='tight')

    # read in trait-specific regions
    features = pd.read_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.csv"), sep="\t")

    # Plot matrix of overlap between sets
    m = pd.DataFrame()  # build matrix of common features
    for t1, t2 in itertools.combinations(traits, 2):
        a = set(features[features['trait'] == t1][['start', 'end']].apply(sum, axis=1))
        b = set(features[features['trait'] == t2][['start', 'end']].apply(sum, axis=1))
        m.loc[t1, t2] = len(a.intersection(b)) / float(len(features[features['trait'] == t1].drop_duplicates()))
        m.loc[t2, t1] = len(a.intersection(b)) / float(len(features[features['trait'] == t2].drop_duplicates()))
    m.replace(np.nan, 1, inplace=True)

    # clustermap on values
    sns.clustermap(np.log2(m.sort(axis=0).sort(axis=1)))
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.common_trait_regions.values.svg")
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close("all")

    # clustermap on correlation
    sns.clustermap(np.log2(m).corr())
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.common_trait_regions.correlation.svg")
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close("all")

    # See how much the regions are positively or negatively correlated with the feature
    g = sns.FacetGrid(features, col="trait", col_wrap=4)
    g.map(sns.distplot, "change")
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.trait-region_association.svg")
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close("all")

    # Get numbers for fractions of regions associated with traits
    trait_associated = features[['chrom', 'start', 'end']].drop_duplicates().shape[0]  # unique trait-associated features
    total = len(analysis.sites)  # unique total features
    variable = analysis.coverage_qnorm_annotated[analysis.coverage_qnorm_annotated["dispersion"] > 0.5].shape[0]  # unique "variable" regions
    # fraction of total discovered regions associated with a clinical trait
    fraction_total = trait_associated / float(total)
    # fraction of variable regions associated with a clinical trait
    fraction_variable = trait_associated / float(variable)
    # write these numbers out
    fig, axis = plt.subplots(2)
    sns.stripplot(["total", "variable"], [total, variable], ax=axis[0])
    axis[0].set_xlim((0, total))
    sns.stripplot(["fraction_total", "fraction_variable"], [fraction_total, fraction_variable], ax=axis[1])
    axis[1].set_xlim((0, 1))
    fig.savefig(os.path.join(analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.fraction_regions_with_traits.svg"), bbox_inches='tight')

    # Build space
    fig, axis = plt.subplots(nrows=2, ncols=5, sharex=True, sharey=True, figsize=(60, 20))
    for i in range(1, len(traits) + 1):
        print(traits[:i])
        # PCA
        # get a matrix of unique features for all samples
        x = features[features['trait'].apply(lambda x: x in traits[:i])].drop_duplicates(['chrom', 'start', 'end'])[
            [s.name for s in analysis.samples if s.library == "ATAC-seq"] + ["trait", "direction"]
        ]
        # save csvs for pca
        pca = PCA()
        x_new = pca.fit_transform(x.T.drop(["trait", "direction"]).apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1))

        # Variable space
        # get one vector loading per feature (sum vectors of each chromatin feature)
        loadings = pd.DataFrame(pca.components_, columns=x.index).T
        loadings['trait'] = x['trait']
        loadings['direction'] = x['direction']

        # sum vectors for each clinical variable, dependent on the direction of the features (positive or negatively associated)
        loadings = pd.melt(loadings.groupby(['trait', 'direction']).sum().reset_index(), ['trait'] + range(x_new.shape[0]), value_name="direction")
        loadings.index = loadings['trait'] + "_" + loadings['direction'].astype(str)
        plot_space(x=pd.DataFrame(x_new), y=loadings.drop(["trait", "variable", "direction"], axis=1), axis=axis.flatten()[i - 1])

    # Save fig
    output_pdf = os.path.join(
        analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.biplot.traits.svg")
    fig.savefig(output_pdf, bbox_inches='tight')

    #

    # Integrate survival into epigenomics space
    # Build final space with all unique features
    x = features.drop_duplicates(['chrom', 'start', 'end'])[[s.name for s in analysis.samples] + ["trait", "direction"]]
    # save csvs for pca
    pca = PCA()
    x_new = pca.fit_transform(x.T.drop(["trait", "direction"]).apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1))

    # plot % explained variance per PC
    fig, axis = plt.subplots(1)
    axis.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_, 'o-')
    axis.set_xlabel("PC")
    axis.set_ylabel("% variance")
    fig.savefig(os.path.join(analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.pca_variance.svg"), bbox_inches='tight')

    # get loadings
    loadings = pd.DataFrame(pca.components_, columns=x.index).T
    loadings['trait'] = x['trait']
    loadings['direction'] = x['direction']

    # sum vectors for each clinical variable, dependent on the direction of the features (positive or negatively associated)
    loadings = pd.melt(loadings.groupby(['trait', 'direction']).sum().reset_index(), ['trait'] + range(x_new.shape[0]), value_name="direction")
    loadings.index = loadings['trait'] + "_" + loadings['direction'].astype(str)

    # plot once more
    plot_space(x=pd.DataFrame(x_new), y=loadings.drop(["trait", "variable", "direction"], axis=1), axis=None)

    # save coordinates on same space
    # z-score both observations and variables
    obs = pd.DataFrame(x_new).apply(lambda j: (j - j.mean()) / j.std(), axis=0)
    var = loadings.drop(["trait", "variable", "direction"], axis=1).apply(lambda j: (j - j.mean()) / j.std(), axis=0)
    # pickle that
    pickle.dump((obs, var), open(os.path.join(analysis.data_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.pickle"), "w"))
    (obs, var) = pickle.load(open(os.path.join(analysis.data_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.pickle"), "r"))

    add_survival(chrom_x=pd.DataFrame(x_new), chrom_y=loadings.drop(["trait", "variable", "direction"], axis=1))


def describe_patients_clinically(analysis):
    """
    """

    traits = ["IGHV", "CD38", "ZAP70", "primary_CLL", "treated", "chemo_treated", "target_treated"]
    muts = ['del11', 'tri12']  # abnormalities
    muts += ['TP53']  # mutations
    traits += muts

    # Here I am removing traits which had poor classification performance or too few patients associated
    # because they exist in either only one patient

    # Read trait-specific chromatin regions in one matrix
    features = pd.read_csv(os.path.join(analysis.data_dir, "trait_specific", "cll.trait-specific_regions.csv"), sep="\t")

    # add direction of chromatin feature association with clinical feature
    features['direction'] = features['change'].apply(lambda x: 1 if x > 0 else -1)

    # sum vectors for each clinical variable, dependent on the direction of the features (positive or negatively associated)
    df = features.groupby(['trait', 'direction'])[[s.name for s in analysis.samples]].aggregate(np.nanmedian)

    # Plot of patient changing across each pos or neg feature
    # df.groupby(level=0).plot()
    df2 = pd.melt(df.reset_index(), ['trait', 'direction'], var_name=["sample"])
    g = sns.FacetGrid(df2, col="trait", hue="sample", margin_titles=True, sharey=False, col_wrap=3, palette="Set2")
    g.map(plt.plot, "direction", "value")
    plt.savefig(os.path.join(analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.sample_change_in_feature.svg"), bbox_inches='tight')
    plt.close("all")

    # Plot rank vs value
    fig, axis = plt.subplots(1)
    num_colors = len(traits)
    cm = plt.get_cmap('Paired')
    axis.set_color_cycle([cm(1. * i / num_colors) for i in range(num_colors)])
    for trait in traits:
        a = df.ix[trait].ix[1] / df.ix[trait].ix[-1]
        axis.plot(range(len(a)), sorted(a, reverse=False), "-o", label=trait)
    axis.legend()
    fig.savefig(os.path.join(analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.sample_rank_features.svg"), bbox_inches='tight')

    # same but independently for each variable and with boxplot
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(
        (len(traits) / 3) * 2, 4,
        width_ratios=[3, 1, 3, 1],
        # height_ratios=[4, 4, 4, 4, 4]
    )

    fig = plt.figure(figsize=(10, 12))
    for i, trait in enumerate([val for val in traits for _ in (0, 1)]):
        a = df.ix[trait].ix[1] / df.ix[trait].ix[-1]
        axis = plt.subplot(gs[i])
        color = plt.get_cmap('Paired')(i)
        if i % 2 == 0:
            axis.plot(range(len(a)), sorted(a, reverse=False), "-o", color=color, label=trait)
            axis.set_title(trait)
        if i % 2 == 1:
            sns.violinplot(a, color=color, orient='v', ax=axis)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis.plots_dir, "trait_specific", "cll_peaks.medical_epigenomics_space.sample_rank_dist_features.svg"), bbox_inches='tight')


def compare_go_terms(cluster_labels, file_names, plot, p_value=0.05, n_max=35):
    """
    """
    for i, _ in enumerate(file_names):
        # read in file
        df = pd.read_csv(file_names[i])
        # label as cluster
        df['cluster'] = cluster_labels[i]
        # append / concatenate
        if i == 0:
            df2 = df
        else:
            df2 = df2.append(df)

    # make readable name
    df2['name'] = (df2['Name'] + " (" + df2['GOID'] + ")").tolist()
    # melt
    df3 = pd.melt(df2, value_vars=['FDR'], id_vars=['name', 'cluster'])
    # expand to two axis terms/cluster, fill with p=1
    df3 = df3.pivot('name', 'cluster', 'value').replace(np.nan, 1)

    # filter p-values
    df3 = df3[(df3[1] < p_value) | (df3[2] < p_value)]

    # floor n_max
    if len(df3) < n_max:
        n_max = len(df3)

    # get n most dissimilar features
    s = abs(df3[1] - df3[2])
    # sort by similarity
    s.sort(ascending=False)
    # get top
    index = s.irow(range(n_max)).index

    # plot matrix of clusters/terms with p_values
    sns.clustermap(df3.ix[index], col_cluster=False, cmap=plt.cm.YlGn_r)
    plt.savefig(plot, bbox_inches='tight')


def compare_pathway_enrichment(cluster_labels, file_names, plot, p_value=0.05, n_max=35):
    """
    """
    for i, _ in enumerate(file_names):
        # read in file
        df = pd.read_csv(file_names[i])
        # label as cluster
        df['cluster'] = cluster_labels[i]
        # append / concatenate
        if i == 0:
            df2 = df
        else:
            df2 = df2.append(df)

    # make readable name
    df2['name'] = (df2['#Term'] + " (" + df2['Database'] + " " + df2['ID'].astype(str) + ")").tolist()
    # melt
    df3 = pd.melt(df2, value_vars=['P-Value'], id_vars=['name', 'cluster'])
    df3.drop_duplicates(inplace=True)
    # expand to two axis terms/cluster, fill with p=1
    df3 = df3.pivot('name', 'cluster', 'value').replace(np.nan, 1)

    # filter p-values
    df4 = df3[(df3[1] < p_value) | (df3[2] < p_value)]

    # get n most dissimilar features
    s = abs(df4[1] - df4[2])
    # sort by similarity
    s.sort(ascending=False)
    # get top
    index = s.irow(range(n_max)).index

    # plot matrix of clusters/terms with p_values
    sns.clustermap(df4.ix[index], col_cluster=False, cmap=plt.cm.YlGn_r)
    plt.savefig(plot, bbox_inches='tight')


def compare_lola_enrichment(cluster_labels, file_names, plot, p_value=20, n_max=35):
    """
    """
    for i, _ in enumerate(file_names):
        # read in file
        df = pd.read_csv(file_names[i], sep="\t")
        # label as cluster
        df['cluster'] = cluster_labels[i]
        # append / concatenate
        if i == 0:
            df2 = df
        else:
            df2 = df2.append(df)

    # make readable name
    df2['name'] = (df2['description'] + " (" + df2['tissue'] + " " + df2['antibody'] + ")").tolist()
    # melt
    df3 = pd.melt(df2, value_vars=['pValueLog'], id_vars=['name', 'cluster'])
    df3.drop_duplicates(['name', 'cluster'], inplace=True)
    # expand to two axis terms/cluster, fill with p=1
    df3 = df3.pivot('name', 'cluster', 'value').replace(np.nan, 1)

    # filter p-values
    df4 = df3[(df3[1] > p_value) | (df3[2] > p_value)]

    # get n most dissimilar features
    s = abs(df3[1] - df3[2])
    # sort by similarity
    s.sort(ascending=True)
    # get top
    index = s.irow(range(n_max)).index

    # plot matrix of clusters/terms with p_values
    sns.clustermap(df4.ix[index], col_cluster=False, cmap=plt.cm.YlGn_r)
    plt.savefig(plot, bbox_inches='tight')


def main():
    # Parse arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = add_args(parser)
    args = parser.parse_args()

    # Start project
    prj = Project("metadata/project_config.yaml")
    prj.add_sample_sheet()

    # annotated samples with a few more things:
    prj.samples = annotate_samples(prj.samples)

    # Start analysis object
    # only with ATAC-seq samples that passed QC
    samples_to_exclude = ["CLL_ATAC-seq_4851_1-5-45960_ATAC29-6_hg19", "CLL_ATAC-seq_AKH13_M3152_ATAC40s21_hg19"]
    samples = [sample for sample in prj.samples if sample.cell_line == "CLL" and sample.name not in samples_to_exclude]

    # Start analysis object
    analysis = Analysis(
        data_dir=os.path.join(".", "data"),
        plots_dir=os.path.join(".", "results", "plots"),
        samples=samples,
        pickle_file=os.path.join(".", "data", "analysis.pickle")
    )
    # pair analysis and Project
    analysis.prj = prj

    atac_seq_samples = [sample for sample in analysis.samples if sample.cell_line == "CLL" and sample.library == "ATAC-seq"]

    # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
    if args.generate:
        # Get consensus peak set from all samples
        analysis.get_consensus_sites(atac_seq_samples)
        # estimate peak saturation among all samples
        analysis.estimate_peak_saturation(n=1000)
        # Calculate peak support
        analysis.calculate_peak_support(atac_seq_samples)
        # Annotate peaks with closest gene
        analysis.get_peak_gene_annotation()
        # Annotate peaks with genomic regions
        analysis.get_peak_genomic_location()
        # Annotate peaks with ChromHMM state from CD19+ cells
        analysis.get_peak_chromatin_state()
    else:
        analysis.sites = pybedtools.BedTool(os.path.join(analysis.data_dir, "cll_peaks.bed"))
        analysis.peak_count = pickle.load(open(os.path.join(analysis.data_dir, "cll_peaks.cum_peak_count.pickle"), "rb"))
        analysis.support = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.support.csv"))
        analysis.gene_annotation = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.gene_annotation.csv"))
        analysis.closest_tss_distances = pickle.load(open(os.path.join(analysis.data_dir, "cll_peaks.closest_tss_distances.pickle"), "rb"))
        analysis.region_annotation = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.region_annotation.csv"))
        analysis.region_annotation_b = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.region_annotation_background.csv"))
        analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.chromatin_state.csv"))
        analysis.chrom_state_annotation_b = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.chromatin_state_background.csv"))

    # Plot general peak set features
    analysis.plot_peak_characteristics()

    # GET CHROMATIN OPENNESS MEASUREMENTS, PLOT
    if args.generate:
        # Get coverage values for each peak in each sample of ATAC-seq and ChIPmentation
        analysis.measure_coverage([sample for sample in analysis.samples if sample.cell_line == "CLL"])
        # normalize coverage values
        analysis.normalize_coverage_quantiles([sample for sample in analysis.samples if sample.cell_line == "CLL"])
        # Annotate peaks with closest gene, chromatin state,
        # genomic location, mean and variance measurements across samples
        analysis.annotate(atac_seq_samples)
    else:
        analysis.coverage = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.raw_coverage.tsv"), sep="\t", index_col=0)
        analysis.coverage_qnorm = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.coverage_qnorm.log2.tsv"), sep="\t")
        analysis.coverage_qnorm_annotated = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")

    # Characterize all CLL regions as a whole
    # region's structure
    characterize_regions_structure(df=analysis.coverage_qnorm_annotated, prefix="all_regions")
    # region's function
    output_dir = os.path.join(analysis.data_dir, "%s_peaks" % "all_regions")
    characterize_regions_function(df=analysis.coverage_qnorm_annotated, output_dir=output_dir, prefix="all_regions")

    # Plot oppenness across peaks/samples
    analysis.plot_coverage()
    analysis.plot_variance()
    analysis.plot_sample_correlations()
    # Observe exponential fit to the coeficient of variation
    analysis.plot_qv2_fit()

    # INTER-PATIENT VARIATION (Figure 2)
    # cross-cohort variation at gene level
    analysis.gene_oppeness_across_samples(atac_seq_samples)
    # try seeing what most variable regions are with LOLA (no good results - not included)
    analysis.variability()
    # inter-group variation
    analysis.inspect_variability(atac_seq_samples)

    # TRAIT-SPECIFIC ANALYSIS (Figure 3)
    # all traits (~21)
    traits = ["IGHV", "CD38", "ZAP70", "treated", "primary_CLL", "ibrutinib", "chemo_treated", "target_treated"]
    muts = ['del11', 'tri12', "del17"]
    muts += ['TP53']  # mutations
    traits += muts
    trait_analysis(analysis, atac_seq_samples, traits)
    # Inspect which traits have good performance
    traits = ["IGHV", "CD38", "ZAP70", "treated", "primary_CLL", "chemo_treated", "target_treated", "del11", "tri12", "del17"]
    # here I am removing traits which had poor classification performance or too few patients associated
    # because they exist in either only one patient
    # Join all regions in one dataframe
    join_trait_specific_regions(analysis, traits)
    # characterize trait-specific regions
    # structurally, functionaly and in light of histone marks
    characterize_regions(analysis, traits)
    characterize_regions_chromatin(analysis, traits)
    # Clinical trait signatures and enrichment of samples in them
    get_signatures(analysis, traits)

    # MEDICAL EPIGENOMIC SPACE
    # Build it (biplot) + add survival
    create_clinical_epigenomic_space(analysis, traits)
    # Describe patient across clinical epigenomics axis
    describe_patients_clinically(analysis)

    #

    # Post-processing

    # Compare go terms/LOLA results from regions
    # compare_go_terms(
    #     [1, 2],
    #     ["data/mutation_peaks_cluster1/mutation_cluster1_regions.seq2pathway.csv', 'data/mutation_peaks_cluster2/mutation_cluster2_regions.seq2pathway.csv'],
    #     "results/plots/mutation_regions.goterm_enrichment.svg",
    #     p_value=0.05)

    # compare_pathway_enrichment(
    #     [1, 2],
    #     ['data/mutation_peaks_cluster1/mutation_cluster1_regions.pathway_enrichment.csv', 'data/mutation_peaks_cluster2/mutation_cluster2_regions.pathway_enrichment.csv'],
    #     "results/plots/mutation_regions.pathway_enrichment.svg",
    #     p_value=0.05)

    # compare_pathway_enrichment(
    #     [1, 2],
    #     ['data/mutation_peaks_cluster1/mutation_cluster1_regions.disease_enrichment.csv', 'data/mutation_peaks_cluster2/mutation_cluster2_regions.disease_enrichment.csv'],
    #     "results/plots/mutation_regions.disease_enrichment.svg",
    #     p_value=0.05)

    # compare_lola_enrichment(
    #     [1, 2],
    #     ['data/mutation_peaks_cluster1/lola/allEnrichments.txt', 'data/mutation_peaks_cluster2/lola/allEnrichments.txt'],
    #     "results/plots/mutation_regions.lola_enrichment.svg",
    #     p_value=0.05)

    # # save "digested" clinical sheet to disk
    # if args.generate:
    #     fields = [
    #         'sampleName', 'diagnosis_disease', 'diagnosis_date', 'collection_date', 'time_since_diagnosis',
    #         'diagnosis_collection', "under_treatment", 'previous_treatment_date', 'time_since_treatment',
    #         'treatment_regimen', 'treatment_response', 'relapse']
    #     prj.sheet.asDataFrame()[fields].drop_duplicates().to_csv("clinical_annotation_digested.csv", index=False)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
