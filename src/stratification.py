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
import pickle


sns.set_style("whitegrid")
sns.set_context("paper")


# decorator for some methods of Analysis class
def pickle_me(function):
    def wrapper(obj):
        function(obj)
        pickle.dump(obj, open(obj.pickle_file, 'wb'))
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

    @pickle_me
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

    @pickle_me
    def annotate_peak_gene(self):
        # create bedtool with hg19 TSS positions
        hg19_ensembl_tss = pybedtools.BedTool(os.path.join(self.data_dir, "GRCh37_hg19_ensembl_genes.tss.bed"))
        # get closest TSS of each cll peak
        closest = self.sites.closest(hg19_ensembl_tss).to_dataframe()[['chrom', 'start', 'end', 'thickStart', 'blockCount']]
        # aggregate annotation per peak, concatenate various genes (comma-separated)
        closest = closest.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set(x))).reset_index()
        closest.columns = ['chrom', 'start', 'end', 'ensembl_gene_id', 'gene_name']
        # save to disk
        closest.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.closest_gene.csv"), index=False)

        self.closest_gene = closest

    @pickle_me
    def annotate_peak_class(self):
        # create bedtool with CD19 chromatin states
        states_cd19 = pybedtools.BedTool(os.path.join(self.data_dir, "E032_15_coreMarks_mnemonics.bed"))
        # intersect with cll peaks, to create annotation, get original peaks
        annotation = self.sites.intersect(states_cd19, wa=True, wb=True).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        # aggregate annotation per peak, concatenate various annotations (comma-separated)
        chrom_states = annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        chrom_states.columns = ['chrom', 'start', 'end', 'chromatin_state']
        # save to disk
        chrom_states.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.chromatin_state.csv"), index=False)

        self.chrom_states = chrom_states

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
        rpkm = pd.concat([rpk[["chrom", "start", "end"]], rpkm], axis=1)
        rpkm.to_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t", index=False)
        rpkm = pd.read_csv(os.path.join(self.data_dir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t")

        self.rpkm = rpkm

    @pickle_me
    def log_rpkm(self):
        # Log2 transform
        self.rpkm_log = self.rpkm
        self.rpkm_log[[sample.name for sample in self.samples]] = np.log2(1 + self.rpkm[[sample.name for sample in self.samples]])

    @pickle_me
    def rpkm_variance(self):
        # add support to rpkm
        self.rpkm = pd.merge(self.rpkm, self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        self.rpkm_log = pd.merge(self.rpkm_log, self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        # mean rpkm
        self.rpkm['mean'] = self.rpkm[[sample.name for sample in self.samples]].apply(lambda x: np.mean(x), axis=1)
        self.rpkm_log['mean'] = self.rpkm_log[[sample.name for sample in self.samples]].apply(lambda x: np.mean(x), axis=1)
        # dispersion (variance / mean)
        self.rpkm['dispersion'] = self.rpkm[[sample.name for sample in self.samples]].apply(lambda x: np.var(x) / np.mean(x), axis=1)
        self.rpkm_log['dispersion'] = self.rpkm_log[[sample.name for sample in self.samples]].apply(lambda x: np.var(x) / np.mean(x), axis=1)
        # qv2 vs mean rpkm
        self.rpkm['qv2'] = self.rpkm[[sample.name for sample in self.samples]].apply(lambda x: (np.std(x) / np.mean(x)) ** 2, axis=1)
        self.rpkm_log['qv2'] = self.rpkm_log[[sample.name for sample in self.samples]].apply(lambda x: (np.std(x) / np.mean(x)) ** 2, axis=1)

    @pickle_me
    def filter_rpkm(self, x):
        self.rpkm_filtered = self.rpkm[self.rpkm['mean'] > x]
        self.rpkm_filtered_log = self.rpkm_log[self.rpkm_log['mean'] > x]

    @pickle_me
    def pca(self):
        # PCA
        # normalize
        x = normalize(self.rpkm_log[[sample.name for sample in self.samples]])

        # random PCA
        pca = RandomizedPCA()
        self.pca_fit = pca.fit(x).transform(x)

    @pickle_me
    def mds(self, n=1000):
        # normalize, get *n* most variable sites
        x = normalize(self.rpkm_log.ix[self.rpkm_log[[sample.name for sample in self.samples]].apply(np.var, axis=1).order(ascending=False).index].head(n)[[sample.name for sample in self.samples]])

        # convert two components as we're plotting points in a two-dimensional plane
        # "precomputed" because we provide a distance matrix
        # we will also specify `random_state` so the plot is reproducible.
        mds = MDS()  # n_components=2, dissimilarity="precomputed", random_state=1)

        self.mds_fir = mds.fit_transform(x)

    def plot_peak_characteristics(self):
        # Plot cumulative number of peaks
        plt.plot(self.peak_count.keys(), self.peak_count.values(), 'o')
        plt.ylim(0, max(self.peak_count.values()) + 5000)
        plt.title("Cumulative peaks per sample")
        plt.xlabel("Number of samples")
        plt.ylabel("Total number of peaks")
        plt.savefig(os.path.join(self.plots_dir, "total_peak_count.per_patient.pdf"), bbox_inches="tight")
        plt.close('all')

        # Loop at summary statistics:
        # interval lengths
        sns.distplot([interval.length for interval in self.sites], bins=300, kde=False)
        plt.xlabel("peak width (bp)")
        plt.ylabel("frequency")
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.lengths.pdf"), bbox_inches="tight")
        plt.close('all')

        # plot support
        sns.distplot(self.support["support"], bins=40)
        plt.ylabel("frequency")
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.support.pdf"), bbox_inches="tight")
        plt.close('all')

    def plot_rpkm(self):
        # Plot
        rpkm_melted = pd.melt(self.rpkm, id_vars=["chrom", "start", "end"], var_name="sample", value_name="rpkm")

        # rpkm density
        # all in one plot
        for sample in self.samples:
            sns.distplot(self.rpkm[[sample.name]], hist=False, label=sample.name)
        # plt.legend()
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.all.pdf"), bbox_inches="tight")
        plt.close()

        # separately in one grid
        g = sns.FacetGrid(rpkm_melted, col="sample", aspect=2, col_wrap=4)
        g.map(sns.distplot, "rpkm", hist=False)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.distplot.pdf"), bbox_inches="tight")
        plt.close()

        # boxplot rpkm per sample
        # Plot the orbital period with horizontal boxes
        sns.boxplot(x="rpkm", y="sample", data=rpkm_melted)
        plt.xlim(0, 15)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.boxplot.pdf"), bbox_inches="tight")
        plt.close()

    def plot_variance(self):
        sns.jointplot('mean', "dispersion", data=self.rpkm_log)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.dispersion.pdf"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('mean', "qv2", data=self.rpkm_log)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.qv2_vs_mean.pdf"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('support', "qv2", data=self.rpkm_log)
        plt.savefig(os.path.join(self.plots_dir, "rpkm_per_sample.support_vs_qv2.pdf"), bbox_inches="tight")
        plt.close('all')

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.rpkm_log[self.rpkm_log[[sample.name for sample in self.samples]].apply(max, axis=1) > 3]

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

        X = np.array(self.rpkm_log['mean'])
        Y = np.array(self.rpkm_log['qv2'])

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

        correlations = self.rpkm_log[[sample.name for sample in self.samples]].corr()
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

    def plot_pca(self):
        from collections import Counter
        # get variance explained by each component
        variance = [np.round(i * 100, 0) for i in self.pca_fit.explained_variance_ratio_]

        # get colors
        # rainbow (unique color per sample)
        colors = cm.Paired(np.linspace(0, 1, len(self.samples)))

        # per patient
        patients = Counter([sample.patientID for sample in analysis.samples]).keys()
        color_dict = cm.Paired(np.linspace(0, 1, len(patients)))
        color_dict = dict(zip(patients, colors))
        colors = [color_dict[sample.patientID] for sample in self.samples]

        # dependent on igvh status
        colors = ['yellow' if sample.mutated else 'green' for sample in self.samples]

        # 2 components
        fig = plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.pca_fit.components_[i, 0], self.pca_fit.components_[i, 1],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        fig.axes[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
        fig.axes[0].set_ylabel("PC2 - {0}% variance".format(variance[1]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-2comp.1vs2.pdf"), bbox_inches="tight")

        # 2vs3 components
        fig = plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.pca_fit.components_[i, 1], self.pca_fit.components_[i, 2],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        fig.axes[0].set_xlabel("PC2 - {0}% variance".format(variance[1]))
        fig.axes[0].set_ylabel("PC3 - {0}% variance".format(variance[2]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-2comp.2vs3.pdf"), bbox_inches="tight")

        # 3vs4 components
        fig = plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.pca_fit.components_[i, 2], self.pca_fit.components_[i, 3],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        fig.axes[0].set_xlabel("PC3 - {0}% variance".format(variance[2]))
        fig.axes[0].set_ylabel("PC4 - {0}% variance".format(variance[3]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-2comp.3vs4.pdf"), bbox_inches="tight")

        # 4vs5 components
        fig = plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.pca_fit.components_[i, 3], self.pca_fit.components_[i, 4],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        fig.axes[0].set_xlabel("PC4 - {0}% variance".format(variance[3]))
        fig.axes[0].set_ylabel("PC5 - {0}% variance".format(variance[4]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-2comp.4vs5.pdf"), bbox_inches="tight")

        # 1vs3 components
        fig = plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.pca_fit.components_[i, 0], self.pca_fit.components_[i, 2],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        fig.axes[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
        fig.axes[0].set_ylabel("PC3 - {0}% variance".format(variance[2]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-2comp.1vs3.pdf"), bbox_inches="tight")

        # 3 components
        fig = plt.figure()
        # plot each point
        ax = fig.add_subplot(111, projection='3d')
        for i, sample in enumerate(self.samples):
            ax.scatter(
                self.pca_fit.components_[i, 0], self.pca_fit.components_[i, 1], self.pca_fit.components_[i, 2],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=100
            )
        ax.set_xlabel("PC1 - {0}% variance".format(variance[0]))
        ax.set_ylabel("PC2 - {0}% variance".format(variance[1]))
        ax.set_zlabel("PC3 - {0}% variance".format(variance[2]))
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.PCA-3comp.pdf"), bbox_inches="tight")

    def plot_mds(self, n=1000):
        # get unique colors
        colors = cm.Paired(np.linspace(0, 1, len(self.samples)))

        # plot
        plt.figure()
        for i, sample in enumerate(self.samples):
            plt.scatter(
                self.mds_fit[i, 0], self.mds_fit[i, 1],
                label=name_to_repr(sample.name),
                color=colors[i],
                s=50
            )
        plt.title("MDS on %i most variable genes" % n)
        plt.legend(loc='center left', ncol=3, bbox_to_anchor=(1, 0.5))
        plt.savefig(os.path.join(self.plots_dir, "all_sample_peaks.concatenated.MDS.pdf"), bbox_inches="tight")
        plt.close('all')


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


def annotate_igvh_mutations(samples, clinical):
    new_samples = list()

    for sample in samples:
        _id = name_to_sample_id(sample.name)
        if clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 1:
            sample.mutated = True
        elif clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 2:
            sample.mutated = False
        else:
            sample.mutated = None
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

# Select ATAC-seq samples
samples = [s for s in prj.samples if type(s) == ATACseqSample if s.cellLine == "CLL"]

# Annotate with igvh mutated status
# add "mutated" attribute to sample depending on IGVH mutation status
samples = annotate_igvh_mutations(samples, clinical)

# Start analysis object
analysis = Analysis(
    data_dir, plots_dir, samples,
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
    analysis.annotate_peak_gene()
else:
    analysis.closest_gene = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.closest_gene.csv"))

# Annotate peaks with ChromHMM state from CD19 cells
if generate:
    analysis.annotate_peak_class()
else:
    analysis.closest_gene = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.chromatin_state.csv"))

# plot general peak set features
if generate:
    analysis.plot_peak_characteristics()


# WORK WITH "OPENNESS"
# Get RPKM values for each peak in each sample
if generate:
    analysis.measure_chromatin_openness()
else:
    analysis.rpkm = pd.read_csv(os.path.join(data_dir, "all_sample_peaks.concatenated.rpkm.tsv"), sep="\t")

# Compute log2
if generate:
    analysis.log_rpkm()

# Get variance measurements across samples
if generate:
    analysis.rpkm_variance()

# Plot rpkm features across peaks/samples
if generate:
    analysis.plot_rpkm()
    analysis.plot_variance()
    analysis.plot_sample_correlations()

# Observe exponential fit to the coeficient of variation
if generate:
    analysis.plot_qv2_fit()

# Decide on low-end cut-off based on the elbow method
if generate:
    analysis.filter_rpkm(1)

# Try to separate samples in 2D space
if generate:
    analysis.pca()
    analysis.plot_pca()
    analysis.mds()
    analysis.plot_mds()

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
significant = analysis.rpkm.ix[pvalues[pvalues < 0.0001].index][[sample.name for sample in samples]]
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
        s = [sample for sample in samples if sample.name == sample_name][0]
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
