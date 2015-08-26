#!/usr/bin/env python

"""
This script makes plots to ilustrate the stratification between CLL patients.
"""

import os
from pipelines.models import Project, ATACseqSample
import pybedtools
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
# from statsmodels.sandbox.stats.multicomp import multipletests
import itertools
import cPickle as pickle
from collections import Counter


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
        # GET CONSENSUS SITES ACROSS SAMPLES
        peak_count = dict()
        for i, sample in enumerate(self.samples):
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
            # Let's keep track of the number of new peaks found with each new sample
            peak_count[i] = len(sites)

        # Merge overlaping peaks across samples
        sites = sites.merge()

        # Filter
        # remove blacklist regions
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "wgEncodeDacMapabilityConsensusExcludable.bed"))
        # remove chrM peaks and save
        sites = sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(self.data_dir, "all_sample_peaks.concatenated.bed"))
        sites = pybedtools.BedTool(os.path.join(self.data_dir, "all_sample_peaks.concatenated.bed"))

        # Store
        self.sites = sites
        self.peak_count = peak_count

    def calculate_peak_support(self):
        # calculate support (number of samples overlaping each merged peak)
        for i, sample in enumerate(self.samples):
            print(sample.name)
            if i == 0:
                support = self.sites.intersect(sample.peaks, wa=True, c=True)
            else:
                support = support.intersect(sample.peaks, wa=True, c=True)
        support = support.to_dataframe()
        support = support.reset_index()
        support.columns = ["chrom", "start", "end"] + [sample.name for sample in self.samples]
        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support[range(len(self.samples))].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(self.samples)), axis=1)
        # save
        support.to_csv(os.path.join(self.data_dir, "all_sample_peaks.support.csv"), index=False)
        support = pd.read_csv(os.path.join(self.data_dir, "all_sample_peaks.support.csv"))

        self.support = support

    def get_peak_gene_annotation(self):
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
        self.gene_annotation.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.gene_annotation.csv"), index=False)

        # save distances to all TSSs (for plotting)
        self.closest_tss_distances = closest['blockCount'].tolist()

    def get_peak_genomic_location(self):
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

        # save to disk
        self.region_annotation.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.region_annotation.csv"), index=False)

        # keep long list of regions for plotting later
        self.all_region_annotation = region_annotation['genomic_region']
        self.all_region_annotation_backround = region_annotation_b['genomic_region']

    def get_peak_chromatin_state(self):
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

        # save to disk
        self.chrom_state_annotation.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.chromatin_state.csv"), index=False)

        # # keep long list of chromatin states (for plotting)
        self.all_chrom_state_annotation = chrom_state_annotation['thickStart'].tolist()
        self.all_chrom_state_annotation_background = chrom_state_annotation_b['thickStart'].tolist()

    @pickle_me
    def measure_chromatin_openness(self):
        # Select ATAC-seq samples
        samples = [s for s in self.prj.samples if type(s) == ATACseqSample]

        # Count reads with pysam
        # make strings with intervals
        sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in self.sites]
        # count, create dataframe
        coverage = pd.DataFrame(
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
        coverage.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.raw_coverage.tsv"), sep="\t", index=True)
        coverage = pd.read_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.raw_coverage.tsv"), sep="\t", index_col=0)

        # Add interval description to df
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            coverage.index
        )
        coverage["chrom"] = [x[0] for x in ints]
        coverage["start"] = [int(x[1]) for x in ints]
        coverage["end"] = [int(x[2]) for x in ints]

        # Normalize by feature length (Reads per kilobase)
        rpk = coverage.apply(normalize_by_interval_length, axis=1)
        rpk.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpk.tsv"), sep="\t", index=True)
        rpk = pd.read_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpk.tsv"), sep="\t", index_col=0)

        # Normalize by library size - mapped reads (Reads per kilobase per million)
        rpkm = rpk[[sample.name for sample in samples]].apply(normalize_by_library_size, args=(samples, ), axis=0)

        # Save
        self.rpkm = pd.concat([rpk[["chrom", "start", "end"]], rpkm], axis=1)

        # calculate log
        self.log_rpkm()

        self.rpkm.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t", index=False)

    def log_rpkm(self):
        # Log2 transform
        self.rpkm[[sample.name for sample in self.samples]] = np.log2(1 + self.rpkm[[sample.name for sample in self.samples]])

    @pickle_me
    def annotate_rpkm(self):
        atacseq_samples = [sample for sample in self.samples if sample.technique == "ATAC-seq"]

        # add closest gene
        self.rpkm_annotated = pd.merge(
            self.rpkm,
            self.gene_annotation[['chrom', 'start', 'end', 'gene_name']], on=['chrom', 'start', 'end'])
        # add genomic location
        self.rpkm_annotated = pd.merge(
            self.rpkm_annotated,
            self.region_annotation[['chrom', 'start', 'end', 'genomic_region']], on=['chrom', 'start', 'end'])
        # add chromatin state
        self.rpkm_annotated = pd.merge(
            self.rpkm_annotated,
            self.chrom_state_annotation[['chrom', 'start', 'end', 'chromatin_state']], on=['chrom', 'start', 'end'])
        # add support to rpkm - this is added here because support was calculated prior to rpkms
        self.rpkm_annotated = pd.merge(
            self.rpkm_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        # calculate mean rpkm
        self.rpkm_annotated['mean'] = self.rpkm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: np.mean(x), axis=1)
        # calculate dispersion (variance / mean)
        self.rpkm_annotated['dispersion'] = self.rpkm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: np.var(x) / np.mean(x), axis=1)
        # calculate qv2
        self.rpkm_annotated['qv2'] = self.rpkm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: (np.std(x) / np.mean(x)) ** 2, axis=1)

        self.rpkm_annotated.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpkm.annotated.tsv"), sep="\t", index=False)

    def filter_rpkm(self, x, method="rpkm"):
        if method == "rpkm":
            self.rpkm_filtered = self.rpkm[self.rpkm_annotated['mean'] > x]
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
        # Plot cumulative number of peaks
        fig, axis = plt.subplots()
        axis.plot(self.peak_count.keys(), self.peak_count.values(), 'o')
        axis.ylim(0, max(self.peak_count.values()) + 5000)
        axis.title("Cumulative peaks per sample")
        axis.xlabel("Number of samples")
        axis.ylabel("Total number of peaks")
        fig.savefig(os.path.join(self.plots_dir, "total_peak_count.per_patient.pdf"), bbox_inches="tight")

        # Loop at summary statistics:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites], bins=300, kde=False, ax=axis)
        axis.xlabel("peak width (bp)")
        axis.ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "all_sample_peaks.lengths.pdf"), bbox_inches="tight")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis)
        axis.ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "all_sample_peaks.support.pdf"), bbox_inches="tight")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots()
        sns.distplot(self.closest_tss_distances, bins=200, ax=axis)
        axis.set_xlabel("distance to nearest TSS (bp)")
        axis.set_ylabel("frequency")
        fig.savefig(os.path.join(self.plots_dir, "all_sample_peaks.tss_distance.pdf"), bbox_inches="tight")

        # Plot genomic regions
        # count region frequency
        count = Counter(self.all_region_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(self.all_region_annotation_backround)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(2, sharex=True)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[1].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")

        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plots_dir, "all_sample_peaks.genomic_regions.pdf"), bbox_inches="tight")

        # Plot chromatin states
        # count region frequency
        count = Counter(self.all_chrom_state_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(self.all_chrom_state_annotation_background)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(2, sharex=True)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[1].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")

        fig.autofmt_xdate()
        fig.tight_layout()
        fig.savefig(os.path.join(self.plots_dir, "all_sample_peaks.chromatin_states.pdf"), bbox_inches="tight")

    def plot_rpkm(self):
        data = self.rpkm_annotated.copy()
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
                     'gene_name', 'genomic_region', 'chromatin_state', 'support'], var_name="sample", value_name="rpkm")

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.mean.per_genomic_region.distplot.pdf"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.mean.chromatin_state.distplot.pdf"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.dispersion.per_genomic_region.distplot.pdf"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.dispersion.chromatin_state.distplot.pdf"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.support.per_genomic_region.distplot.pdf"), bbox_inches="tight")
        plt.close()

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        plt.savefig(os.path.join(self.plots_dir, "rpkm.support.chromatin_state.distplot.pdf"), bbox_inches="tight")
        plt.close()

        #

        # Beware below!

        # rpkm density
        # all in one plot
        for sample in self.samples:
            sns.distplot(self.rpkm_annotated[[sample.name]], hist=False, label=sample.name)
        # plt.legend()
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.all.pdf"), bbox_inches="tight")
        plt.close()

        # separately in one grid
        g = sns.FacetGrid(data_melted, col="sample", aspect=2, col_wrap=4)
        g.map(sns.distplot, "rpkm", hist=False)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.pdf"), bbox_inches="tight")
        plt.close()

        # boxplot rpkm per sample
        # Plot the orbital period with horizontal boxes
        sns.boxplot(x="rpkm", y="sample", data=data_melted)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.boxplot.pdf"), bbox_inches="tight")
        plt.close()

    def plot_variance(self):
        sns.jointplot('mean', "dispersion", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.dispersion.pdf"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('mean', "qv2", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.qv2_vs_mean.pdf"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('support', "qv2", data=self.rpkm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.support_vs_qv2.pdf"), bbox_inches="tight")
        plt.close('all')

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.rpkm_annotated[self.rpkm_annotated[[sample.name for sample in self.samples]].apply(max, axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.dispersion.filtered.pdf"), bbox_inches="tight")
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.support_vs_qv2.filtered.pdf"), bbox_inches="tight")

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

        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.qv2_vs_mean.fit_residuals.pdf"), bbox_inches="tight")

    def plot_sample_correlations(self):
        # get colors depending on IGVH mut
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
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.correlation_clustering.pdf"), bbox_inches="tight")
        plt.close('all')

    def plot_pca(self, suffix=""):
        # get variance explained by each component
        variance = [np.round(i * 100, 0) for i in self.pca.explained_variance_ratio_]

        # dependent on igvh status
        colors = samples_to_color(samples)

        # plot
        fig, axis = plt.subplots(nrows=2, ncols=2)
        fig.set_figheight(10)
        fig.set_figwidth(10)
        axis = axis.flatten()

        # 1vs2 components
        for i, sample in enumerate(self.samples):
            axis[0].scatter(
                self.pca.components_[i, 0], self.pca.components_[i, 1],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        axis[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
        axis[0].set_ylabel("PC2 - {0}% variance".format(variance[1]))

        # 2vs3 components
        for i, sample in enumerate(self.samples):
            axis[1].scatter(
                self.pca.components_[i, 1], self.pca.components_[i, 2],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        axis[1].set_xlabel("PC2 - {0}% variance".format(variance[1]))
        axis[1].set_ylabel("PC3 - {0}% variance".format(variance[2]))

        # 3vs4 components
        for i, sample in enumerate(self.samples):
            axis[2].scatter(
                self.pca.components_[i, 2], self.pca.components_[i, 3],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        axis[2].set_xlabel("PC3 - {0}% variance".format(variance[2]))
        axis[2].set_ylabel("PC4 - {0}% variance".format(variance[3]))

        # 4vs5 components
        for i, sample in enumerate(self.samples):
            axis[3].scatter(
                self.pca.components_[i, 3], self.pca.components_[i, 4],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        axis[3].set_xlabel("PC4 - {0}% variance".format(variance[3]))
        axis[3].set_ylabel("PC5 - {0}% variance".format(variance[4]))

        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plot_path = os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA_{0}.pdf".format(suffix))
        fig.savefig(plot_path, bbox_inches="tight")

        # 3 components
        fig = plt.figure()
        # plot each point
        ax = fig.add_subplot(111, projection='3d')
        for i, sample in enumerate(self.samples):
            ax.scatter(
                self.pca.components_[i, 0], self.pca.components_[i, 1], self.pca.components_[i, 2],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        ax.set_xlabel("PC1 - {0}% variance".format(variance[0]))
        ax.set_ylabel("PC2 - {0}% variance".format(variance[1]))
        ax.set_zlabel("PC3 - {0}% variance".format(variance[2]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plot_path = os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA_{0}_3comp.pdf".format(suffix))
        fig.savefig(plot_path, bbox_inches="tight")

    def plot_mds(self, n, suffix=""):
        # get unique colors
        colors = samples_to_color(samples)

        # plot
        fig, axis = plt.subplots(1)
        for i, sample in enumerate(self.samples):
            axis.scatter(
                self.mds_fit[i, 0], self.mds_fit[i, 1],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        axis.set_title("MDS on %i most variable regions" % n)
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plot_path = os.path.join(self.plots_dir, "all_sample_peaks.concatenated.MDS_{0}.pdf".format(suffix))
        fig.savefig(plot_path, bbox_inches="tight")


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


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_sample_id(name):
    return name.split("_")[3:4][0]


def samples_to_color(samples, method="mutation"):
    # dependent on igvh status
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
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            # get sample id
            _id = name_to_sample_id(sample.name)
            # get corresponding series from "clinical"
            sample_c = clinical[clinical['sample_id'] == _id].squeeze()
            # add timepoint
            sample.timepoint = sample_c.timepoint
            # add time since diagnosis
            sample.time_since_diagnosis = pd.to_datetime(sample_c['treatment_%i_date' % sample.timepoint]) - pd.to_datetime(sample_c['diagnosis_date'])
            if sample_c['treated'] == "Y":
                sample.treatment_active = True
                # if sample is treated, find out which treatment based on timepoint
                sample.treatment_type = sample_c['treatment_%i_regimen' % t]
                # CR, GR, PR, NR in this order of 'goodness'
                sample.treatment_response = sample_c['treatment_%i_response' % t]
            elif sample_c['treated'] == "N":
                sample.treatment_active = False
                sample.treatment_type = None
                sample.treatment_response = None
            elif pd.isnull(sample_c.treated):
                sample.treatment_active = None
                sample.treatment_type = None
                sample.treatment_response = None
        else:
            sample.timepoint = None
            sample.treatment_active = None
            sample.treatment_type = None
            sample.treatment_response = None
            sample.time_since_diagnosis = None
        new_samples.append(sample)
    return new_samples


def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)


def get_cluster_classes(den, label='ivl'):
    from collections import defaultdict

    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][j] for j in l]
        cluster_classes[c] = i_l

    return cluster_classes


def run_lola(bed_files, universe_file, output_folder):
    import rpy2.robjects as robj

    lola = robj.r("""
        function(bedFiles, universeFile, outputFolder) {
            library("bedr")
            library("LOLA")

            userUniverse  <- bedr::bed_to_granges(universeFile)

            dbPath = "/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/"

            dbPath = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/hg19/"
            regionDB = loadRegionDB(dbPath)

            if (typeof(bedFiles) == "character") {
                userSet <- bedr::bed_to_granges(bedFiles)
                lolaResults = runLOLA(userSet, userUniverse, regionDB, cores=12)
                lolaResults[order(support, decreasing=TRUE), ]
                writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
            } else if (typeof(bedFiles) == "double") {
                for (bedFile in bedFiles) {
                    userSet <- bedr::bed_to_granges(bedFile)
                    lolaResults = runLOLA(userSet, userUniverse, regionDB, cores=12)
                    lolaResults[order(support, decreasing=TRUE), ]
                    writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
                }
            }
        }
    """)

    # convert the pandas dataframe to an R dataframe
    lola(bed_files, universe_file, output_folder)


# Should we regenerate the data?
generate = False

# Get path configuration
data_dir = os.path.join('.', "data")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")

# Get clinical info
clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))
attributes = [
    "patient_id", "sample_id", "timepoint",
    "igvh_mutation_status", "patient_gender", "sample_viability",
    "patient_birth_date", "diagnosis_date"
]
clinical = clinical[attributes].drop_duplicates()

# Start project
# prj = pickle.load(open("prj.pickle", 'rb'))
prj = Project("cll-patients")
prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

# Annotate with igvh mutated status
# add "mutated" attribute to sample depending on IGVH mutation status
# Select ATAC-seq samples
prj.samples = annotate_igvh_mutations(prj.samples, clinical)

# Annotate with clinical data
prj.samples = annotate_treatments(prj.samples, clinical)

# Start analysis object
# only with ATAC-seq samples
analysis = Analysis(
    data_dir, plots_dir, [sample for sample in prj.samples if sample.technique == "ATAC-seq"],
    pickle_file=os.path.join(data_dir, "analysis.pickle")
)
analysis.prj = prj
analysis.clinical = clinical


# GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
# Get consensus peak set from all samples
if generate:
    analysis.get_consensus_sites()
else:
    analysis.sites = pybedtools.BedTool(os.path.join(data_dir, "all_sample_peaks.concatenated.bed"))

# Calculate peak support
if generate:
    analysis.calculate_peak_support()
else:
    analysis.support = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.support.csv"))

# Annotate peaks with closest gene
if generate:
    analysis.get_peak_gene_annotation()
else:
    analysis.gene_annotation = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.gene_annotation.csv"))

# Annotate peaks with genomic regions
if generate:
    analysis.get_peak_genomic_location()
else:
    analysis.region_annotation = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.region_annotation.csv"))

# Annotate peaks with ChromHMM state from CD19 cells
if generate:
    analysis.get_peak_chromatin_state()
else:
    analysis.chrom_state_annotation = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.chromatin_state.csv"))


# WORK WITH "OPENNESS"
# Get RPKM values for each peak in each sample
if generate:
    analysis.measure_chromatin_openness()
else:
    analysis.rpkm = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t")

# Annotate peaks with closest gene, chromatin state,
# genomic location, mean and variance measurements across samples
if generate:
    analysis.annotate_rpkm()
else:
    analysis.rpkm_annotated = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.rpkm.annotated.tsv"), sep="\t")


# plot general peak set features
if generate:
    analysis.plot_peak_characteristics()

# Plot rpkm features across peaks/samples
if generate:
    analysis.plot_rpkm()
    analysis.plot_variance()
    analysis.plot_sample_correlations()

# Observe exponential fit to the coeficient of variation
if generate:
    analysis.plot_qv2_fit()


# Try to separate samples in 2D space
if generate:
    # PCA
    # Decide on low-end threshold based on the elbow method
    t = 2
    analysis.filter_rpkm(t, method="rpkm")

    n = 3
    analysis.filter_rpkm(3, method="support")

    data = pd.merge(analysis.rpkm_filtered, analysis.chrom_states, on=['chrom', 'start', 'end'])

    analysis.pca_analysis(data[[sample.name for sample in analysis.samples]])
    analysis.plot_pca(suffix="mean>%i" % t)

    # Filter peaks based on sd
    n = 100
    atacseq_samples = [sample for sample in analysis.samples if sample.technique == "ATAC-seq"]
    most_variable_data = data.sort(["dispersion"], ascending=False).head(n)
    analysis.pca_analysis(most_variable_data[[sample.name for sample in analysis.samples]])
    analysis.plot_pca(suffix="mostdispersion")

    most_variable_data = data.sort(["qv2"], ascending=False).head(n)
    analysis.pca_analysis(most_variable_data[[sample.name for sample in analysis.samples]])
    analysis.plot_pca(suffix="mostqv2")

    # Use only Promoters
    promoter_data = data.ix[data['chromatin_state'].str.contains("Tss")]
    analysis.pca_analysis(promoter_data[[sample.name for sample in analysis.samples]])
    analysis.plot_pca(suffix="promoters")

    # Use only Enhancers
    enhancer_data = data.ix[data['chromatin_state'].str.contains("Enh")]
    analysis.pca_analysis(enhancer_data[[sample.name for sample in analysis.samples]])
    analysis.plot_pca(suffix="enhancers")

    # MDS
    # X most variable sites
    analysis.mds_analysis(most_variable_data[[sample.name for sample in analysis.samples]])
    analysis.plot_mds(n)

# COMPARISON mCLL - uCLL
# test if sites come from same population based on rpkm values
pvalues = analysis.rpkm.apply(
    lambda x: mannwhitneyu(
        x.loc[[sample.name for sample in samples if sample.mutated]],
        x.loc[[sample.name for sample in samples if not sample.mutated]])[1],
    axis=1
)

# correct for multiple testing
# qvalues = pd.Series(multipletests(pvalues)[1])

# get differential sites
significant = analysis.rpkm.ix[pvalues[pvalues < 0.0001].index][[sample.name for sample in analysis.samples]]
significant.to_csv(os.path.join(data_dir, "dors.mutated_vs_unmutated.tsv"), sep="\t", index=False)

# correlate samples on significantly different sites
sns.clustermap(
    significant.corr(),
    method="complete",
    square=True, annot=False,
    figsize=(20, 16)
)
plt.savefig(os.path.join(plots_dir, "dors.correlation_clustering.pdf"), bbox_inches="tight")

# cluster samples and sites
# plot heatmap of differentialy open sites
clustermap = sns.clustermap(
    significant,
    cmap=plt.get_cmap('YlGn'),
    square=False, annot=False,
)
plt.savefig(os.path.join(plots_dir, "dors.clustering.pdf"), bbox_inches="tight")

# get cluster assignments from linkage matrix
lm = clustermap.dendrogram_col.linkage

# plot dendrogram
# determine height to separate clusters
dendr = dendrogram(lm, labels=significant.columns, color_threshold=65)
plt.savefig(os.path.join(plots_dir, "dors.dendrogram.pdf"), bbox_inches="tight")

# plt.show()

# assign each sample to one cluster
clusters = get_cluster_classes(dendr)

# concatenate clusters on the unmutated side
unmutated_cluster = ['b', 'c', 'm', 'y']

# annotate samples with cluster
new_samples = list()
for cluster, sample_names in clusters.items():
    for sample_name in sample_names:
        s = [sample for sample in analysis.samples if sample.name == sample_name][0]
        s.cluster = cluster if cluster not in unmutated_cluster else 'b'
        new_samples.append(s)
samples = new_samples

# Repeat again independence test and
# get all differential sites (from 3 comparisons)
all_pvalues = pd.DataFrame()
for g1, g2 in itertools.combinations(['r', 'g', 'b'], 2):
    pvalues = pd.DataFrame(analysis.rpkm.apply(
        lambda x: mannwhitneyu(
            x.loc[[sample.name for sample in samples if sample.cluster == "r"]],
            x.loc[[sample.name for sample in samples if sample.cluster == "g"]])[1],
        axis=1
    ))
    pvalues['comparison'] = "-".join([g1, g2])
    all_pvalues = pd.concat([all_pvalues, pvalues])

all_pvalues.columns = ['p', 'comparison']
all_pvalues.to_csv(os.path.join("data", "dors.3_populations.pvalues.tsv"), sep="\t", index=False)

all_significant = analysis.rpkm.ix[all_pvalues[all_pvalues['p'] < 0.000001].index].drop_duplicates()

all_significant.to_csv(os.path.join("data", "dors.3_populations.significant.tsv"), sep="\t", index=False)

# correlate samples on significantly different sites
sns.clustermap(
    all_significant.corr(),
    method="complete",
    square=True, annot=False,
    figsize=(20, 16)
)
plt.savefig(os.path.join(plots_dir, "dors.3_populations.correlation_clustering.pdf"), bbox_inches="tight")

# heatmap all significat sites
clustermap = sns.clustermap(
    all_significant[[sample.name for sample in samples]],
    cmap=plt.get_cmap('YlGn'),
    square=False, annot=False
)
plt.savefig(os.path.join(plots_dir, "dors.3_populations.clustering.pdf"), bbox_inches="tight")

# export dors regions
for g1, g2 in itertools.combinations(['r', 'g', 'b'], 2):
    comparison = "-".join([g1, g2])
    specific = analysis.rpkm.loc[
        all_pvalues[
            (all_pvalues['p'] < 0.000001) &
            (all_pvalues['comparison'] == comparison)
        ].index,
        ['chrom', 'start', 'end']
    ].drop_duplicates()
    specific.to_csv(os.path.join(data_dir, "dors.{0}.bed".format(comparison)), sep="\t", index=False, header=None)

# run lola
# use all cll sites as universe
universe_file = os.path.join(data_dir, "all_sample_peaks.concatenated.bed")

for g1, g2 in itertools.combinations(['r', 'g', 'b'], 2):
    comparison = "-".join([g1, g2])
    bed_file = os.path.join(data_dir, "dors.{0}.bed".format(comparison))
    output_folder = os.path.join(data_dir, "lola", "dors_{0}".format(comparison))
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # run
    run_lola(bed_file, universe_file, output_folder)

# get closest gene: GO, KEGG, OMIM, mSigDB
# de novo motif finding - enrichment

# Subset data in Enhancers/Promoters
# stratify again

# INTER-SAMPLE VARIABILITY ANALYSIS
# Subsample peaks or reads and see the minimum required to form the clusters previously
# See which regions explain most of variability for each cluster

# CLASSIFICATION
# Stratify patients on:
# - treated vs untreated
# - ...
# Train classifiers on groups
# Predict for all samples
# Assess: ROC, AUC
