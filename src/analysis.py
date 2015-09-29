#!/usr/bin/env python

"""
This is the main script of the cll-patients project.
"""

import recipy
from argparse import ArgumentParser
import os
import sys
from pipelines.models import Project, ATACseqSample
import pybedtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing
import parmap
import pysam
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.decomposition import RandomizedPCA
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests
import cPickle as pickle
from collections import Counter

# Set settings
pd.set_option("date_dayfirst", True)
sns.set_style("whitegrid")
sns.set_context("paper")


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
    def get_consensus_sites(self):
        # GET CONSENSUS SITES ACROSS CLL ATAC-SEQ SAMPLES
        samples = [sample for sample in self.samples if sample.cellLine == "CLL" and sample.technique == "ATAC-seq"]

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

    def estimate_peak_saturation(self, n=100):
        """
        Randomizes sample order n times and measures the number of unique peaks after adding each sample after the other.
        Plots this.
        """
        import random
        # GET CONSENSUS SITES ACROSS CLL ATAC-SEQ SAMPLES
        samples = [sample for sample in self.samples if sample.cellLine == "CLL" and sample.technique == "ATAC-seq"]

        peak_count = dict()
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
                if j not in peak_count.keys():
                    peak_count[j] = [len(sites)]
                else:
                    peak_count[j].append(len(sites))

        # Plot
        fig, axis = plt.subplots(1)
        df = pd.DataFrame(peak_count)
        df.apply(np.mean).plot(label="mean", ax=axis)
        df.apply(np.min).plot(label="min", ax=axis)
        df.apply(np.max).plot(label="max", ax=axis)
        axis.set_xlabel("# samples")
        axis.set_ylabel("# peaks")
        axis.set_ylim((0, 120000))
        plt.legend(loc='best')
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.cum_peak_count.svg"), bbox_inches="tight")

        # Store
        df.to_csv(os.path.join(self.data_dir, "cll_peaks.cum_peak_count.csv"), index=False)

    def calculate_peak_support(self):
        samples = [s for s in self.samples if s.technique == "ATAC-seq" and s.cellLine == "CLL"]
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

        gene_annotation = pd.merge(gene_annotation, ensembl_gtn)
        self.gene_annotation = pd.merge(gene_annotation, ensembl_gtp)

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
    def measure_coverage(self):
        # Select ATAC-seq samples
        samples = [s for s in self.prj.samples if type(s) == ATACseqSample]

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
    def normalize_coverage(self):
        # Normalize by feature length (Reads per kilobase)
        rpk = self.coverage.apply(normalize_by_interval_length, axis=1)

        # Normalize by library size - mapped reads (Reads per kilobase per million)
        rpkm = rpk[[sample.name for sample in self.samples]].apply(normalize_by_library_size, args=(self.samples, ), axis=0)

        # Save
        self.rpkm = pd.concat([rpk[["chrom", "start", "end"]], rpkm], axis=1)

        # calculate log
        self.log_rpkm()

        self.rpkm.to_csv(os.path.join(self.data_dir, "cll_peaks.rpkm.tsv"), sep="\t", index=False)

    def log_rpkm(self):
        # Log2 transform
        self.rpkm[[sample.name for sample in self.samples]] = np.log2(1 + self.rpkm[[sample.name for sample in self.samples]])

    @pickle_me
    def normalize_coverage_quantiles(self):
        # Normalize by quantiles
        to_norm = self.coverage.iloc[:, :len(self.samples)]
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
    def annotate(self):
        atacseq_samples = [sample for sample in self.samples if sample.cellLine == "CLL"]

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
        # add support to rpkm - this is added here because support was calculated prior to rpkms
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        # calculate mean rpkm
        self.coverage_qnorm_annotated['mean'] = self.coverage_qnorm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: np.mean(x), axis=1)
        # calculate variance
        self.coverage_qnorm_annotated['variance'] = self.coverage_qnorm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: np.var(x), axis=1)
        # calculate std deviation (sqrt(variance))
        self.coverage_qnorm_annotated['std_deviation'] = np.sqrt(self.coverage_qnorm_annotated['variance'])
        # calculate dispersion (variance / mean)
        self.coverage_qnorm_annotated['dispersion'] = self.coverage_qnorm_annotated['variance'] / self.coverage_qnorm_annotated['mean']
        # calculate qv2 (std / mean) ** 2
        self.coverage_qnorm_annotated['qv2'] = (self.coverage_qnorm_annotated['std_deviation'] / self.coverage_qnorm_annotated['mean']) ** 2

        # calculate "fold-change" (max - min)
        self.coverage_qnorm_annotated['fold_change'] = self.coverage_qnorm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: x.max() - x.min(), axis=1)

        self.coverage_qnorm_annotated.to_csv(os.path.join(self.data_dir, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index=False)

    def filter_rpkm(self, x, method):
        if method == "rpkm":
            self.rpkm_filtered = self.rpkm_annotated[self.rpkm_annotated['mean'] > x]
        elif method == "std":
            # this assumes x represents x times the value
            # of the standard deviation for a given peak
            self.rpkm_filtered = self.rpkm_annotated[self.rpkm_annotated['std_deviation'] > x]
        elif method == "support":
            # this assumes x represents a minimum of samples
            # therefore we need to calculate
            n = float(len([s for s in self.samples if (s.cellLine == "CLL" and s.technique == "ATAC-seq")]))
            self.rpkm_filtered = self.rpkm_annotated[self.rpkm_annotated['support'] > x / n]

    def pca_analysis(self, data):
        # PCA
        # normalize
        x = normalize(data)

        # random PCA
        self.pca = RandomizedPCA()
        self.pca_fit = self.pca.fit(x).transform(x)

    def mds_analysis(self, data):
        # normalize, get *n* most variable sites
        x = normalize(data)

        # convert two components as we're plotting points in a two-dimensional plane
        # "precomputed" because we provide a distance matrix
        # we will also specify `random_state` so the plot is reproducible.
        self.mds = MDS()  # n_components=2, dissimilarity="precomputed", random_state=1)

        self.mds_fit = self.mds.fit_transform(x)

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

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[1].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")

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

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[1].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")

        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plots_dir, "cll_peaks.chromatin_states.svg"), bbox_inches="tight")

    def plot_rpkm(self):
        # data = self.rpkm_annotated.copy()
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
        # rpkm
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

        #

        # Beware below!

        # rpkm density
        # all in one plot
        for sample in self.samples:
            sns.distplot(self.rpkm_annotated[[sample.name]], hist=False, label=sample.name)
        # plt.legend()
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.all.svg"), bbox_inches="tight")
        plt.close()

        # separately in one grid
        g = sns.FacetGrid(data_melted, col="sample", aspect=2, col_wrap=4)
        g.map(sns.distplot, "rpkm", hist=False)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.svg"), bbox_inches="tight")
        plt.close()

        # boxplot rpkm per sample
        # Plot the orbital period with horizontal boxes
        sns.boxplot(x="rpkm", y="sample", data=data_melted)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.boxplot.svg"), bbox_inches="tight")
        plt.close()

    def plot_variance(self):
        sns.jointplot('mean', "dispersion", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.dispersion.svg"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('mean', "qv2", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.qv2_vs_mean.svg"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('support', "qv2", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.support_vs_qv2.svg"), bbox_inches="tight")
        plt.close('all')

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.rpkm_annotated[self.rpkm_annotated[[sample.name for sample in self.samples]].apply(max, axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.dispersion.filtered.svg"), bbox_inches="tight")
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.support_vs_qv2.filtered.svg"), bbox_inches="tight")

    def plot_qnorm_comparison(self):
        # Compare raw counts vs qnormalized data
        fig, axis = plt.subplots(2, sharex=True)
        [sns.distplot(np.log2(1 + self.coverage[[sample.name]]), ax=axis[0], hist=False) for sample in self.samples if sample.cellLine != "PBMC"]
        [sns.distplot(self.coverage_qnorm[[sample.name]], ax=axis[1], hist=False) for sample in self.samples if sample.cellLine != "PBMC"]
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
        plt.close('all')


def add_args(parser):
    """
    Options for project and pipelines.
    """
    # Behaviour
    parser.add_argument("-g", "--generate", dest="stats", action="store_true",
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
        colors = list()
        for sample in samples:
            if sample.mutated is True:
                colors.append('yellow')
            elif sample.mutated is False:
                colors.append('green')
            elif sample.mutated is None:
                if sample.cellLine == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
        return colors
    # unique color per patient
    elif method == "unique":
        # per patient
        patients = set([sample.patientID for sample in samples])
        color_dict = cm.Paired(np.linspace(0, 1, len(patients)))
        color_dict = dict(zip(patients, colors))
        return [color_dict[sample.patientID] for sample in samples]
    # rainbow (unique color per sample)
    elif method == "unique":
        return cm.Paired(np.linspace(0, 1, len(samples)))
    elif method == "gender":
        colors = list()
        for sample in samples:
            if sample.patient_gender is "F":
                colors.append('red')
            elif sample.patient_gender is "M":
                colors.append('blue')
            elif sample.patient_gender is None:
                if sample.cellLine == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
        return colors
    elif method == "treatment":
        drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR', 'Alemtuz']
        colors = list()
        for sample in samples:
            if sample.treatment_active is False and sample.relapse is False:
                colors.append('dodgerblue')
            elif sample.treatment_active is True and sample.treatment_type in drugs:
                colors.append('peru')
            else:
                colors.append('grey')
        return colors
    elif method == "disease":
        colors = list()
        for sample in samples:
            if sample.diagnosis_disease == "CLL":
                colors.append('#A6CEE3')
            elif sample.diagnosis_disease == "MBL":
                colors.append('#F17047')
            else:
                colors.append('grey')
        return colors
    else:
        raise ValueError("Method %s is not valid" % method)


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


def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)


def get_cluster_classes(dendrogram, label='ivl'):
    from collections import defaultdict

    cluster_idxs = defaultdict(list)
    for c, pi in zip(dendrogram['color_list'], dendrogram['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [dendrogram[label][j] for j in l]
        cluster_classes[c] = i_l

    return cluster_classes


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


def meme(input_fasta, output_dir):
    """
    De novo motif finding with MEME-chip.
    """
    cmd = """meme-chip \\
    -meme-p 12 \\
    -ccut 147 \\
    -meme-minw 6 -meme-maxw 30 -meme-nmotifs 20 \\
    -dreme-e 0.05 \\
    -centrimo-score 5.0 \\
    -centrimo-ethresh 10.0 \\
    -ccut 147 \\
    -db ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv9.meme -meme-mod zoops \\
    -oc {0} \\
    {1}
    """.format(output_dir, input_fasta)
    return cmd


def group_analysis(analysis, sel_samples, feature, g1, g2, group1, group2, generate=True):
    """
    Analysis between two groups of samples.
    :param analysis: an Analysis object.
    :param sel_samples: a list of selected samples to perform analysis on.
    :param g1: a pandas.DataFrame with openness values.
    :param g2: a pandas.DataFrame with openness values.
    :param group1: a string representing the g1 group.
    :param group2: a string representing the g2 group.
    """
    go_term_mapping = os.path.join(analysis.prj.dirs.data, "goID_goName.csv")

    # get colors depending on feature (gender, mutation, drugs, etc...)
    if feature == "mutated":
        method = "mutation"
    elif feature == "patient_gender":
        method = "gender"
    elif feature == "untreated_vs_1stline":
        method = "treatment"
    elif feature == "CLL_vs_MBL":
        method = "disease"
    sample_colors = samples_to_color(sel_samples, method=method)

    # ANNOTATE
    # compute p-value, add to annotation
    fold_change = list()
    mean_a = list()
    mean_b = list()
    p_values = list()
    for i in range(len(analysis.coverage_qnorm_annotated)):
        # compute log2 fold-change
        a = g1.ix[i].mean()
        b = g2.ix[i].mean()
        mean_a.append(a)
        mean_b.append(b)
        fold_change.append(np.log2(a / b))
        # compute p-value
        p_values.append(mannwhitneyu(g1.ix[i], g2.ix[i])[1])

    # add to annotation
    analysis.coverage_qnorm_annotated["_".join(["fold_change", feature])] = fold_change
    analysis.coverage_qnorm_annotated["_".join(["mean", feature, str(group1)])] = mean_a
    analysis.coverage_qnorm_annotated["_".join(["mean", feature, str(group2)])] = mean_b
    # compute q values
    q_values = pd.Series(multipletests(p_values)[1])
    analysis.coverage_qnorm_annotated["_".join(["p_value", feature])] = p_values
    analysis.coverage_qnorm_annotated["_".join(["q_value", feature])] = q_values
    # save as csv
    csv = os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_mean_foldchange_pvalues.csv" % method)
    analysis.coverage_qnorm_annotated[[
        "_".join(["fold_change", feature]),
        "_".join(["mean", feature, str(group1)]),
        "_".join(["mean", feature, str(group2)]),
        "_".join(["q_value", feature]),
        "_".join(["p_value", feature])]].to_csv(csv, index=False)
    # pickle
    # analysis.to_pickle()

    # VISUALIZE
    # visualize distribution of fold-change, p-values
    # A vs B
    sns.jointplot(np.log2(1 + np.array(mean_a)), np.log2(1 + np.array(mean_b)))
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_fold_change.svg" % method), bbox_inches="tight")
    plt.close('all')
    # volcano plot (logfoldchange vs logpvalue)
    sns.jointplot(
        np.array(fold_change), -np.log10(np.array(p_values)),
        stat_func=None, space=0, xlim=(-2.5, 2.5)
    )
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_volcano.svg" % method), bbox_inches="tight")
    plt.close('all')

    # get significant sites
    if feature == "mutated":
        significant = analysis.coverage_qnorm_annotated[
            (analysis.coverage_qnorm_annotated["_".join(["p_value", feature])] < 0.0001) &
            (abs(analysis.coverage_qnorm_annotated["_".join(["fold_change", feature])]) > 1)
        ]
    elif feature == "patient_gender":
        significant = analysis.coverage_qnorm_annotated[
            (analysis.coverage_qnorm_annotated["_".join(["p_value", feature])] < 0.000001)
        ]
    else:
        significant = analysis.coverage_qnorm_annotated[
            (analysis.coverage_qnorm_annotated["_".join(["p_value", feature])] < 0.0001)
        ]

    # SAVE AS BED
    bed_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.bed" % method)
    significant[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)

    # EXPLORE
    # get normalized counts in significant sites only for CLL samples
    significant_values = significant[[sample.name for sample in sel_samples]]

    # get nice sample IDs
    significant_values.columns = map(name_to_repr, significant_values.columns)

    # cluster samples and sites
    # plot heatmap of differentialy open sites
    sites_cluster = sns.clustermap(
        significant_values,
        cmap=plt.get_cmap('YlGn'),
        annot=False,
        standard_scale=0,
        yticklabels=False,
        col_colors=sample_colors
    )
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.svg" % method), bbox_inches="tight")
    plt.close('all')

    # correlate samples on significantly different sites
    corr_cluster = sns.clustermap(
        significant_values.corr(),
        method="complete",
        annot=False,
        col_colors=sample_colors,
        row_colors=sample_colors
    )
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_correlation.svg" % method), bbox_inches="tight")
    plt.close('all')

    # pca
    from sklearn.preprocessing import normalize
    x = significant_values.values  # returns a numpy array
    x_scaled = normalize(x)

    output_pdf = os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.pca.svg" % method)
    pca_r(x_scaled, sample_colors, output_pdf)

    # CHARACTERIZE
    # plot features of all sites
    significant['length'] = significant.apply(lambda x: x['end'] - x['start'], axis=1)

    # plot interval lengths and support
    for variable in ['length', 'support']:
        fig, axis = plt.subplots()
        sns.distplot(significant[variable], bins=300, kde=False, ax=axis)
        fig.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.peak_%s.svg" % (method, variable)), bbox_inches="tight")

    # # Plot distance to nearest TSS
    # fig, axis = plt.subplots()
    # sns.distplot(significant['closest_tss_distances'], bins=200, ax=axis)
    # fig.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.closest_tss_distances.svg" % method), bbox_inches="tight")

    # plot genomic location and chromatin state
    for variable in ['genomic_region', 'chromatin_state']:
        fig, axis = plt.subplots()
        count = Counter([y for x in significant.apply(lambda x: x[variable].split(","), axis=1) for y in x])
        data = pd.DataFrame([count.keys(), count.values()]).T.reset_index(drop=True)
        data.sort([1], ascending=False, inplace=True)
        sns.barplot(x=0, y=1, data=data, ax=axis)
        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.%s.svg" % (method, variable)), bbox_inches="tight")

    # TODO:
    # Get CLL expression data
    # Divide significant in <0 and >0 fold_change
    # Correlate with expression

    # Lola
    # use all cll sites as universe
    universe_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.bed")
    output_folder = os.path.join(analysis.prj.dirs.data, "lola", "cll_peaks.%s_significant.clustering_sites" % method)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # run
    lola(bed_file, universe_file, output_folder)

    # seq2pathway
    # export file with ID, chrom, start, end
    tsv_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.tsv" % method)
    export = significant[['chrom', 'start', 'end']]
    export['index'] = export.index
    export[['index', 'chrom', 'start', 'end']].to_csv(tsv_file, sep="\t", header=False, index=False)

    results = seq2pathway(tsv_file, go_term_mapping)

    results_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.csv" % method)
    results.to_csv(results_file, index=False)

    # GO Terms
    # write gene names to file
    all_cll_genes_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.closest_genes.txt")
    with open(all_cll_genes_file, 'w') as handle:
        for gene in analysis.coverage_qnorm_annotated['gene_name']:
            handle.write(gene + "\n")
    feature_genes_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.closest_genes.txt" % method)
    with open(feature_genes_file, 'w') as handle:
        for gene in significant['gene_name']:
            handle.write(gene + "\n")
    output_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.closest_genes.go_enrichment.tsv" % method)
    # test enrichements of closest gene function: GO, KEGG, OMIM
    goverlap(feature_genes_file, all_cll_genes_file, output_file)

    # Motifs
    # de novo motif finding - enrichment
    bed_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.bed" % method)
    fasta_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.fa" % method)
    bed_to_fasta(bed_file, fasta_file)
    output_folder = os.path.join(analysis.prj.dirs.data, "meme", "cll_peaks.%s_significant" % method)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    meme(fasta_file, output_folder)

    # SIGNIFICANT SITE SUBSET CARACTERIZATION (CLUSTERS)
    # get dendrogram data and plot it,
    # determine height to separate clusters
    # this must be done empirically for each feature type
    thresholds = {"mutated": 12, "patient_gender": 8, "untreated_vs_1stline": 8, "CLL_vs_MBL": 8}

    dendr = dendrogram(sites_cluster.dendrogram_row.linkage, color_threshold=thresholds[feature], labels=significant_values.index)  # labels have to be reset for some reason... grrrr!
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.dendrogram.svg" % method), bbox_inches="tight")

    # assign a cluster to each peak
    site_cluster_dict = dict(zip(dendr['ivl'], dendr['color_list']))
    significant['cluster'] = [site_cluster_dict[x] if x in site_cluster_dict.keys() else None for x in significant.index]

    # for each cluster, test equality of features with all significant peaks
    # or count genomic regions/chromatin states
    p_values = pd.DataFrame()
    counts = pd.DataFrame()
    for cluster_name, cluster in enumerate(significant['cluster'].unique()):
        if cluster is None:
            continue

        # GET DATA FROM CLUSTER
        # grab data only from this cluster
        cluster_data = significant[significant['cluster'] == cluster]

        if len(cluster_data) < 30:
            continue

        # SAVE AS BED
        bed_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.cluster_%i.bed" % (method, cluster_name))
        cluster_data[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)

        # TEST DIFFERENCES IN PEAK CARACTERISTICS FROM WHOLE SET
        # continuous variables, test difference of means
        for variable in ['length', 'support']:
            # test and append to df
            series = pd.Series(
                [cluster_name, variable, mannwhitneyu(cluster_data[variable], significant[variable])[1]],
                index=['cluster', 'variable', 'p-value'])
            p_values = p_values.append(series, ignore_index=True)

        # categorical variables, plot side by side
        for variable in ['genomic_region', 'chromatin_state']:
            for data in ['significant', 'cluster_data']:
                significant_count = Counter([y for x in eval(data).apply(lambda x: x[variable].split(","), axis=1) for y in x])
                df = pd.DataFrame(significant_count.values())
                df['data'] = data
                df['values'] = significant_count.keys()
                df['cluster'] = cluster_name
                df['variable'] = variable
                counts = counts.append(df)

        # Lola
        # use all cll sites as universe
        universe_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.bed")
        output_folder = os.path.join(analysis.prj.dirs.data, "lola", "cll_peaks.%s_significant.clustering_sites.cluster_%i" % (method, cluster_name))
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # run
        lola(bed_file, universe_file, output_folder)

        # seq2pathway
        # export file with ID, chrom, start, end
        tsv_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.cluster_%i.tsv" % (method, cluster_name))
        export = significant[['chrom', 'start', 'end']]
        export['index'] = export.index
        export[['index', 'chrom', 'start', 'end']].to_csv(tsv_file, sep="\t", header=False, index=False)

        results = seq2pathway(tsv_file, go_term_mapping)

        results_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.cluster_%i.csv" % (method, cluster_name))
        results.to_csv(results_file, index=False)

        # GO Terms
        # write gene names to file
        all_cll_genes_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.closest_genes.txt")
        with open(all_cll_genes_file, 'w') as handle:
            for gene in analysis.coverage_qnorm_annotated['gene_name']:
                handle.write(gene + "\n")
        feature_genes_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.closest_genes.cluster_%i.txt" % (method, cluster_name))
        with open(feature_genes_file, 'w') as handle:
            for gene in significant['gene_name']:
                handle.write(gene + "\n")
        output_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.closest_genes.go_enrichment.cluster_%i.tsv" % (method, cluster_name))
        # test enrichements of closest gene function: GO, KEGG, OMIM
        goverlap(feature_genes_file, all_cll_genes_file, output_file)

        # Motifs
        # de novo motif finding - enrichment
        bed_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.cluster_%i.bed" % (method, cluster_name))
        fasta_file = os.path.join(analysis.prj.dirs.data, "cll_peaks.%s_significant.clustering_sites.cluster_%i.fa" % (method, cluster_name))
        bed_to_fasta(bed_file, fasta_file)
        output_folder = os.path.join(analysis.prj.dirs.data, "meme", "cll_peaks.%s_significant.cluster_%i" % (method, cluster_name))
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        meme(fasta_file, output_folder)

    # plot enrichments
    df = df.rename(columns={0: "counts"})
    g = sns.FacetGrid(counts, col="cluster", row="variable", hue="data", sharex=False, margin_titles=True)
    g.map(sns.barplot, "values", 0)
    plt.legend(loc="best")
    g.set_axis_labels(x_var="cluster #", y_var="count")
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.clusters.enrichments.svg" % method), bbox_inches="tight")

    # plot p-values
    p_values['p-value'] = -np.log10(p_values['p-value'])
    p_values.sort('p-value', ascending=False, inplace=True)
    g = sns.FacetGrid(p_values, col="variable", margin_titles=True)
    g.map(sns.barplot, "cluster", "p-value")
    # add sig line
    for axis in g.axes[0]:
        axis.axhline(-np.log10(0.05), linestyle='- -', color='black')
    g.set_axis_labels(x_var="cluster #", y_var="-log10(pvalue)")
    plt.savefig(os.path.join(analysis.prj.dirs.plots, "cll_peaks.%s_significant.clustering_sites.clusters.length_support_p-values.svg" % method), bbox_inches="tight")


def main():
    # Parse arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = add_args(parser)

    # Parse
    args = parser.parse_args()

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

    # save "digested" clinical sheet to disk
    if args.generate:
        fields = [
            'sampleName', 'diagnosis_disease', 'diagnosis_date', 'collection_date', 'time_since_diagnosis',
            'diagnosis_collection', "treatment_active", 'previous_treatment_date', 'time_since_treatment',
            'treatment_type', 'treatment_response', 'relapse']
        prj.sheet.asDataFrame()[fields].drop_duplicates().to_csv("clinical_annotation_digested.csv", index=False)

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
    samples = [sample for sample in prj.samples if sample.technique == "ATAC-seq" and sample.name not in samples_to_exclude]

    analysis = Analysis(
        data_dir, plots_dir, samples,
        pickle_file=os.path.join(data_dir, "analysis.pickle")
    )
    analysis.prj = prj
    analysis.clinical = clinical

    # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
    # GET CHROMATIN OPENNESS MEASUREMENTS
    # PLOT STUFF
    if args.generate:
        # Get consensus peak set from all samples
        analysis.get_consensus_sites()
        # estimate peak saturation among all samples
        analysis.estimate_peak_saturation(n=25)
        # Calculate peak support
        analysis.calculate_peak_support()
        # Annotate peaks with closest gene
        analysis.get_peak_gene_annotation()
        # Annotate peaks with genomic regions
        analysis.get_peak_genomic_location()
        # Annotate peaks with ChromHMM state from CD19+ cells
        analysis.get_peak_chromatin_state()

        # WORK WITH "OPENNESS"
        # Get coverage values for each peak in each sample
        analysis.measure_coverage()
        # normalize coverage values
        analysis.normalize_coverage()
        analysis.normalize_coverage_quantiles()
        # Annotate peaks with closest gene, chromatin state,
        # genomic location, mean and variance measurements across samples
        analysis.annotate()

        # Plots
        # plot general peak set features
        analysis.plot_peak_characteristics()
        # Plot rpkm features across peaks/samples
        analysis.plot_rpkm()
        analysis.plot_variance()
        analysis.plot_sample_correlations()
        # Observe exponential fit to the coeficient of variation
        analysis.plot_qv2_fit()
        # PCA & MDS analysis:
        # currently not working
    else:
        analysis.sites = pybedtools.BedTool(os.path.join(analysis.prj.dirs.data, "cll_peaks.bed"))
        analysis.peak_count = pickle.load(open(os.path.join(analysis.prj.dirs.data, "cll_peaks.cum_peak_count.pickle"), 'rb'))

        analysis.support = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.support.csv"))

        analysis.gene_annotation = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.gene_annotation.csv"))
        analysis.closest_tss_distances = pickle.load(open(os.path.join(analysis.prj.dirs.data, "cll_peaks.closest_tss_distances.pickle"), 'rb'))

        analysis.region_annotation = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.region_annotation.csv"))
        analysis.region_annotation_b = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.region_annotation_background.csv"))

        analysis.chrom_state_annotation = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.chromatin_state.csv"))
        analysis.chrom_state_annotation_b = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.chromatin_state_background.csv"))

        analysis.coverage = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.raw_coverage.tsv"), sep="\t", index_col=0)
        # analysis.rpkm = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.rpkm.tsv"), sep="\t")
        analysis.coverage_qnorm = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.coverage_qnorm.log2.tsv"), sep="\t")

        analysis.coverage_qnorm_annotated = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")
        # analysis.rpkm_annotated = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.rpkm.annotated.tsv"), sep="\t")

    # TRAIT-SPECIFIC ANALYSIS
    # Use these samples only
    sel_samples = [sample for sample in analysis.samples if sample.cellLine == "CLL" and sample.name not in samples_to_exclude and sample.technique == "ATAC-seq"]

    # "gender" and "mutated"
    features = {
        "patient_gender": ("F", "M"),  # gender
        "mutated": (True, False),  # ighv mutation
    }
    for i, (feature, (group1, group2)) in enumerate(features.items()):
        # example : i, (feature, (group1, group2)) = (1, (features.items()[1]))
        # get dataframe subset with groups
        g1 = analysis.coverage_qnorm_annotated[[sample.name for sample in sel_samples if getattr(sample, feature) == group1]]
        g2 = analysis.coverage_qnorm_annotated[[sample.name for sample in sel_samples if getattr(sample, feature) == group2]]
        group_analysis(analysis, sel_samples, feature, g1, g2, group1, group2)

    # untreated vs 1st line chemotherapy +~ B cell antibodies
    g1 = analysis.coverage_qnorm_annotated[
        [sample.name for sample in sel_samples if not sample.treatment_active and not sample.relapse]
    ]
    drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR', 'Alemtuz']
    g2 = analysis.coverage_qnorm_annotated[
        [sample.name for sample in sel_samples if sample.treatment_active and sample.treatment_type in drugs]
    ]
    group_analysis(analysis, sel_samples, "untreated_vs_1stline", g1, g2, "untreated", "1stlinetreatment")

    # Disease at Diagnosis - comparison in untreated samples
    # CLL vs MBL
    g1 = analysis.coverage_qnorm_annotated[
        [sample.name for sample in sel_samples if sample.diagnosis_disease == "CLL" and not sample.treatment_active and not sample.relapse]
    ]
    g2 = analysis.coverage_qnorm_annotated[
        [sample.name for sample in sel_samples if sample.diagnosis_disease == "MBL" and not sample.treatment_active and not sample.relapse]
    ]
    group_analysis(analysis, sel_samples, "CLL_vs_MBL", g1, g2, "CLL", "MBL")

    # "relapse", ("True", "False"), # relapse or before relapse
    # "treatment_1st", ("untreated", "Chlor"),  # untreated vs 1st line chemotherapy
    # "treatment_2nd", ("untreated", "Ibrut"),  # untreated vs ibrutinib
    # "diagnosis_start", ("CLL", "MBL"),
    # possibly other groups:
    # ['SF3B1', 'ATM', 'del13', 'del11q', 'tri12', 'NOTCH1', 'BIRC3', 'BCL2', 'TP53', 'MYD88', 'CHD2', 'NFKIE']

    # iCLL ANALYSIS
    # get cluster assignments from linkage matrix
    # select 'intermediate' cluster
    # Repeat again independence test and
    # get all differential sites (from 3 comparisons)

    # MINIMUM ELEMENT ANALYSIS
    # Subsample peaks or reads and see the minimum required to form the clusters previously
    # See which regions explain most of variability for each cluster

    # SAMPLE CLASSIFICATION
    # Stratify patients on above groups
    # Train classifiers on groups
    # Predict for all samples
    # Assess: ROC, AUC

    # COMPARE CLL WITH NORMAL CELLS
    # we'll have matched B-cells from the patients as well

    # DE NOVO/CLL-SPECIFIC ENHANCERS
    # Find unique enhancers across CLL samples compared with normal B-cells
    # Search other cell types for overlaping enhancers:
    # - if possitive -> enhancer activation
    # - if negative -> de novo enhancer -> explore mechanism
    # validate with H3K27ac ChIP-seq
    # validate with RNA expression


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
