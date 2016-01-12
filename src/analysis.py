#!/usr/bin/env python

"""
This is the main script of the cll-patients project.
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import recipy
from argparse import ArgumentParser
import os
import sys
from pipelines.models import Project, ATACseqSample
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
sns.set_style("whitegrid")
sns.set_context("paper")
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('text', usetex='false')


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

    def estimate_peak_saturation(self, n=1000):
        """
        Randomizes sample order n times and measures the number of unique peaks after adding each sample after the other.
        Plots this.
        """
        import random
        from scipy import stats

        # GET CONSENSUS SITES ACROSS CLL ATAC-SEQ SAMPLES
        samples = [sample for sample in self.samples if sample.cellLine == "CLL" and sample.technique == "ATAC-seq"]

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
        samples = [sample for sample in self.samples if sample.cellLine == "CLL" and sample.technique == "ATAC-seq"]

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
        # add support to coverage - this is added here because support was calculated prior to coverage
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'])
        # calculate mean coverage
        self.coverage_qnorm_annotated['mean'] = self.coverage_qnorm_annotated[[sample.name for sample in atacseq_samples]].apply(lambda x: np.mean(x), axis=1)
        # calculate coverage variance
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

    def inspect_variability(self):
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
        ighv_u = [s.name for s in self.samples if not s.mutated]
        ighv_m = [s.name for s in self.samples if s.mutated]

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

    def gene_oppeness_across_samples(self):
        """
        Annotates peaks with closest gene.
        Needs files downloaded by prepare_external_files.py
        """
        from collections import OrderedDict
        from scipy.stats import ks_2samp

        sns.set(style="whitegrid", palette="pastel", color_codes=True)

        names = [s.name for s in self.samples]

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
        gene_values = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)][[sample.name for sample in self.samples]].T
        gene_values.columns = promoters.ix[gene_values.columns]['gene_name']
        promoter_data = pd.melt(gene_values, var_name="gene", value_name="openness")
        promoter_data['region'] = 'promoter'

        gene_values = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)][[sample.name for sample in self.samples]].T
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
        plt.close('all')

        sns.jointplot('mean', "qv2", data=self.coverage_qnorm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.qv2_vs_mean.svg"), bbox_inches="tight")
        plt.close('all')

        sns.jointplot('support', "qv2", data=self.coverage_qnorm_annotated)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.support_vs_qv2.svg"), bbox_inches="tight")
        plt.close('all')

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated[[sample.name for sample in self.samples]].apply(max, axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.dispersion.filtered.svg"), bbox_inches="tight")
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.plots_dir, "norm_counts_per_sample.support_vs_qv2.filtered.svg"), bbox_inches="tight")

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
                if sample.cellLine == "CLL":
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
                if sample.cellLine == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
                return colors
    # unique color per patient
    elif method == "unique":
        # per patient
        patients = set(sample.patientID for sample in samples)
        color_dict = cm.Paired(np.linspace(0, 1, len(patients)))
        color_dict = dict(zip(patients, color_dict))
        return [color_dict[sample.patientID] for sample in samples]
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
                if sample.cellLine == "CLL":
                    colors.append('gray')
                else:
                    colors.append('black')
        return colors
    elif method == "treatment":
        drugs = ['Alemtuz', "Ibrutinib"]
        colors = list()
        for sample in samples:
            if sample.treatment_active is False and sample.relapse is False:
                colors.append('peru')  # untreated samples
            elif sample.treatment_active is True and sample.treatment_type not in drugs:
                colors.append('black')  # chemo
            elif sample.treatment_active is True and sample.treatment_type in drugs:
                colors.append('green')  # antibodies
            elif sample.treatment_active is True and sample.treatment_type == "":
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
        patients = set(sample.patientID for sample in samples)
        symbol_dict = [c.next()[0] for _ in range(len(patients))]
        symbol_dict = dict(zip(patients, symbol_dict))
        return [symbol_dict[sample.patientID] for sample in samples]
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


def annotate_IGHV(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            h = clinical.loc[clinical['sample_id'] == _id, 'igvh_homology']
            if len(h) != 1:
                sample.ighv_homology = None
                sample.mutated = None
            else:
                sample.ighv_homology = h.tolist()[0]
                if clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 1:
                    sample.mutated = True
                elif clinical.loc[clinical['sample_id'] == _id, 'igvh_mutation_status'].tolist()[0] == 2:
                    sample.mutated = False
                else:
                    sample.mutated = None
        else:
            sample.ighv_homology = None
            sample.mutated = None
        new_samples.append(sample)
    return new_samples


def annotate_CD38(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            c = clinical.loc[clinical['sample_id'] == _id, 'CD38_cells_percentage']

            if len(c) != 1:
                sample.CD38_percentage = None
                sample.CD38 = None
            else:
                sample.CD38_percentage = c.tolist()[0]
                if clinical.loc[clinical['sample_id'] == _id, 'CD38_positive'].tolist()[0] == 2:
                    sample.CD38 = True
                elif clinical.loc[clinical['sample_id'] == _id, 'CD38_positive'].tolist()[0] == 1:
                    sample.CD38 = False
                else:
                    sample.CD38 = None
        else:
            sample.CD38_percentage = None
            sample.CD38 = None
        new_samples.append(sample)
    return new_samples


def annotate_ZAP70(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)

            z = clinical.loc[clinical['sample_id'] == _id, 'ZAP70_cells_percentage']

            if len(z) != 1:
                sample.ZAP70_percentage = None
                sample.ZAP70 = None
            else:
                sample.ZAP70_percentage = z.tolist()[0]

                # ZAP70 expression
                if clinical.loc[clinical['sample_id'] == _id, 'ZAP70_positive'].tolist()[0] == 2:
                    sample.ZAP70 = True
                elif clinical.loc[clinical['sample_id'] == _id, 'ZAP70_positive'].tolist()[0] == 1:
                    sample.ZAP70 = False
                else:
                    sample.ZAP70 = None

            # ZAP70 mono-allelic methylation/expression
            zm = clinical.loc[clinical['sample_id'] == _id, 'ZAP70_monoallelic_methylation']
            if len(zm) != 1:
                sample.ZAP70_monoallelic = False
            else:
                if zm.tolist()[0] == "Y":
                    sample.ZAP70_monoallelic = True
                else:
                    sample.ZAP70_monoallelic = False
        else:
            sample.ZAP70_percentage = None
            sample.ZAP70 = None
            sample.ZAP70_monoallelic = None
        new_samples.append(sample)
    return new_samples


def annotate_disease_treatments(samples, clinical):
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

            # Get patient timepoint
            sample.timepoint = sample_c['timepoint']
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
                        # if sample is collected after treatment end, mark as not under treatment
                        if treatment_date < treatment_end_date < sample.collection_date:
                            sample.treatment_active = False
                        # if sample is collected before treatment end, mark as under treatment
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
                setattr(sample, attr, pd.np.nan)

        # Append sample
        new_samples.append(sample)
    return new_samples


def annotate_mutations(samples, clinical):
    new_samples = list()

    clinical2 = clinical.replace(to_replace=["q", "loss", "\?"], value=[""] * 3, regex=True)

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            if clinical2.loc[clinical2['sample_id'] == _id, 'mutations'].empty:
                sample.mutations = []
            else:
                m = clinical2[clinical2['sample_id'] == _id]['mutations']
                if len(m) == 1:
                    if type(m.tolist()[0]) is str:
                        sample.mutations = m.tolist()[0].split(",")
                    else:
                        sample.mutations = None
                else:
                    sample.mutations = None
        else:
            sample.mutations = None
        new_samples.append(sample)
    return new_samples


def annotate_gender(samples, clinical):
    new_samples = list()

    for sample in samples:
        if sample.cellLine == "CLL" and sample.technique == "ATAC-seq":
            _id = name_to_sample_id(sample.name)
            g = clinical.loc[clinical['sample_id'] == _id, 'patient_gender']

            if len(g) != 1:
                sample.patient_gender = None
            else:
                if g.tolist()[0] == "F":
                    sample.patient_gender = "F"
                elif clinical.loc[clinical['sample_id'] == _id, 'patient_gender'].tolist()[0] == "M":
                    sample.patient_gender = "M"
                else:
                    sample.patient_gender = None
        else:
            sample.patient_gender = None
        new_samples.append(sample)
    return new_samples


def annotate_samples(samples, clinical):
    samples = annotate_IGHV(samples, clinical)
    samples = annotate_CD38(samples, clinical)
    samples = annotate_ZAP70(samples, clinical)
    samples = annotate_disease_treatments(samples, clinical)
    samples = annotate_gender(samples, clinical)
    samples = annotate_mutations(samples, clinical)
    return samples


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


def characterize_regions_composition(df, prefix, universe_df=None, plots_dir="results/plots"):
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


def classify_samples(analysis, sel_samples, labels, comparison, rerun=False):
    """
    Use a machine learning approach for sample classification based on known sample attributes.
    Extract features most important to separate samples and investigate those.
    """
    from tqdm import tqdm

    print("Trait:%s" % comparison)
    print("%i samples with trait annotated" % len(sel_samples))
    print(Counter(labels))

    dataframe_file = os.path.join(
        analysis.prj.dirs.data,
        "trait_specific",
        "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % comparison)
    if os.path.exists(dataframe_file) and not rerun:  # Load up
        dataframe = pd.read_csv(dataframe_file, sep="\t")
    else:  # Run analysis
        # Get all CLL ATAC-seq samples for validation
        all_samples = [s for s in analysis.samples if s.cellLine == "CLL" and s.technique == "ATAC-seq"]

        # Get colors depending on this feature label (comparison) (True == green; False == redish)
        palette = sns.color_palette("colorblind")
        comparison_colors = [palette[1] if l else palette[2] for l in labels]
        all_samples_colors = [comparison_colors[sel_samples.index(s)] if s in sel_samples else "gray" for s in all_samples]
        cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)  # for later usage in heatmaps

        # ALL CLL OPEN CHROMATIN REGIONS
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]]

        # BINARY CLASSIFICATION
        # get features and labels
        X = normalize(matrix).T
        y = np.array(labels)

        loo = cross_validation.LeaveOneOut(len(X))

        for i, (train_index, test_index) in tqdm(enumerate(loo)):
            # Skip validations on samples from the same patient
            # # get patient:sample mappings
            # p_s_mapping = dict(pd.DataFrame([pd.Series(s.__dict__) for s in sel_samples]).groupby('patientID').groups.items())
            # group = p_s_mapping.index([x for x in p_s_mapping.values() if train_index in x][0])
            # if test_index in p_s_mapping[group][:].pop(group):
            #     continue

            # print(i)
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
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.ROC_PRC.svg" % comparison), bbox_inches="tight")

        # Display training and prediction of pre-labeled samples of most informative features:
        # average feature importance across iterations
        mean_importance = importance.mean(axis=0)

        # visualize feature importance
        fig, axis = plt.subplots(1)
        sns.distplot(mean_importance, ax=axis)
        fig.savefig(os.path.join(
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.mean_importance.svg" % comparison), bbox_inches="tight")
        plt.close('all')

        # get important features
        # n = 500; x = matrix.loc[np.argsort(mean_importance)[-n:], :] # get n top features
        # or
        # Get most informative features
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]]
        x = matrix.loc[[i for i, j in enumerate(mean_importance > 0) if j == True], :]  # get features on the tail of the importance distribution
        sites_cluster = sns.clustermap(
            x,
            cmap=cmap,
            standard_scale=0,
            col_colors=comparison_colors,
            yticklabels=False)
        plt.savefig(os.path.join(
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.svg" % comparison), bbox_inches="tight")
        plt.close('all')

        # SEE ALL SAMPLES
        # Display most informative features for ALL samples:
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in all_samples]]
        # get important features
        x = matrix.loc[[i for i, j in enumerate(mean_importance > 1e-4) if j == True], :]  # get features on the tail of the importance distribution

        # get colors for each cluster
        group_number = 4 if comparison == 'IGHV' else 2  # 4 IGHV groups, everything else 2
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
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sample_correlation.svg" % comparison), bbox_inches="tight")
        plt.close('all')

        # pca on these regions
        pca_r(x, colors, os.path.join(
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.pca.sample_labels.svg" % comparison))

        # REGIONS
        # Characterize regions as a whole
        # region's composition
        characterize_regions_composition(df=analysis.coverage_qnorm_annotated.loc[x.index, :], prefix=comparison)
        # region's function
        output_dir = os.path.join(analysis.data_dir, "trait_specific", "%s_peaks" % comparison)
        characterize_regions_function(df=analysis.coverage_qnorm_annotated.loc[x.index, :], output_dir=output_dir, prefix=comparison)

        # Split in major groups
        if comparison == 'gender':
            gr = 2
        elif comparison == 'IGHV':
            gr = 4
        else:
            gr = 5
        Z = sites_cluster.dendrogram_row.linkage
        clusters = fcluster(Z, gr, criterion="maxclust")
        # visualize  cluster site attribution
        # get cluster colors
        cluster_colors = dict(zip(np.unique(clusters), sns.color_palette("colorblind")))
        colors = [cluster_colors[c] for c in clusters]
        sns.clustermap(
            x,
            cmap=cmap,
            standard_scale=0,
            col_colors=comparison_colors,
            row_colors=colors,
            yticklabels=False)
        plt.savefig(os.path.join(
            analysis.prj.dirs.plots,
            "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sites_labels.svg" % comparison), bbox_inches="tight")
        plt.close('all')

        # get cluster labels for sites
        # add cluster number to df
        dataframe = analysis.coverage_qnorm_annotated.loc[x.index, :]
        dataframe['importance'] = mean_importance[x.index]
        dataframe['cluster'] = clusters
        # Save whole dataframe as csv
        dataframe_file = os.path.join(
            analysis.prj.dirs.data,
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % comparison)
        dataframe.to_csv(dataframe_file, sep="\t", index=False)

        # Save as bed
        bed_file = os.path.join(
            analysis.prj.dirs.data,
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.sites.bed" % comparison)
        dataframe[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)

    # Region characterization
    # plot chromosome distribution of regions
    chrom_count = Counter(dataframe[['chrom', 'start', 'end']].chrom)
    chrom_count = np.array(sorted(chrom_count.items(), key=lambda x: x[1], reverse=True))
    fig, axis = plt.subplots(1, figsize=(10, 5))
    sns.barplot(chrom_count[:, 0], chrom_count[:, 1].astype(int), ax=axis)
    fig.savefig(os.path.join(
        analysis.prj.dirs.plots,
        "trait_specific", "cll_peaks.%s_significant.classification.random_forest.loocv.sites_location.svg" % comparison), bbox_inches="tight")

    for i, cluster in enumerate(np.unique(dataframe['cluster'])):
        # GET REGIONS FROM CLUSTER
        df = dataframe[dataframe['cluster'] == cluster]

        # ignore clusters with less than 20 regions
        if len(df) < 20:
            continue

        # output folder
        outdir = os.path.join("%s_cluster%i" % (comparison, cluster))
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # region's composition
        regions = characterize_regions_composition(df=df, prefix="%s_cluster%i" % (comparison, cluster), universe_df=analysis.coverage_qnorm_annotated)
        regions['cluster'] = cluster

        if i == 0:
            df3 = regions
        else:
            df3 = df3.append(regions)

        # region's function
        output_dir = os.path.join(analysis.data_dir, "trait_specific", "%s_peaks_cluster%i" % (comparison, cluster))
        characterize_regions_function(df=df, output_dir=output_dir, prefix="%s_cluster%i" % (comparison, cluster))

        # # parse meme-ame output
        # motifs = pd.DataFrame(parse_ame(os.path.join(output_dir, "meme")), columns=['motifs', 'q_values'])
        # motifs['cluster'] = cluster
        # if i == 0:
        #     df2 = motifs
        # else:
        #     df2 = pd.concat([df2, motifs])

    # # Plot region enrichment
    # df3['region'] = df3.index
    # df3 = df3.sort(['region'])
    # df3 = df3.replace(np.nan, 0)

    # g = sns.FacetGrid(col="region", data=df3[df3['variable'] == 'genomic_region'], col_wrap=3, sharey=True)
    # g.map(sns.barplot, "cluster", "value")
    # plt.savefig(os.path.join(
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.region_enrichment.svg" % comparison), bbox_inches="tight")

    # g = sns.FacetGrid(col="region", data=df3[df3['variable'] == 'chromatin_state'], col_wrap=3, sharey=True)
    # g.map(sns.barplot, "cluster", "value")
    # plt.savefig(os.path.join(
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.chromatin_enrichment.svg" % comparison), bbox_inches="tight")

    # # save data
    # df3.to_csv(os.path.join(
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.region_enrichment.csv" % comparison), index=False)

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
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.motif_enrichment.svg" % comparison), bbox_inches="tight")
    # plt.close('all')

    # sns.clustermap(df3[(df3.icol(0) > 30) & (df3.icol(1) > 30)], cmap=plt.cm.YlGn)
    # plt.savefig(os.path.join(
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.motif_enrichment.highest.svg" % comparison), bbox_inches="tight")
    # plt.close('all')

    # # save data
    # df3.to_csv(os.path.join(
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.motif_enrichment.csv" % comparison), index=False)

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
    #     analysis.prj.dirs.plots,
    #     "trait_specific", "%s_regions.motif_enrichment.difference.svg" % comparison), bbox_inches="tight")
    # plt.close('all')


def trait_analysis(analysis):
    # Gender
    sel_samples = [s for s in analysis.samples if type(s.patient_gender) is str]
    labels = np.array([1 if s.patient_gender == "M" else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="gender")

    # IGHV mutation status
    sel_samples = [s for s in analysis.samples if type(s.mutated) is bool]
    labels = np.array([1 if s.mutated else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="IGHV")

    # CD38 expression
    sel_samples = [s for s in analysis.samples if type(s.CD38) is bool]
    labels = np.array([1 if s.CD38 else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="CD38")

    # ZAP70 expression
    sel_samples = [s for s in analysis.samples if type(s.ZAP70) is bool]
    labels = np.array([1 if s.ZAP70 else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="ZAP70")

    # ZAP70 mono-allelic expression
    sel_samples = [s for s in analysis.samples if type(s.ZAP70_monoallelic) == bool]
    labels = np.array([1 if s.ZAP70_monoallelic else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="ZAP70_monoallelic")

    # Disease at Diagnosis - comparison in untreated samples
    # Primary vs Secondary CLL (progressed from MBL and SLL)
    sel_samples = [s for s in analysis.samples if type(s.diagnosis_disease) is str]
    labels = np.array([1 if s.diagnosis_disease == "CLL" else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="primary_CLL")

    # Treatment types
    # Under treatment
    sel_samples = [s for s in analysis.samples if type(s.treatment_active) is bool]
    labels = np.array([1 if s.treatment_active else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="treated")

    # Under treament but with different types
    # Chemotherapy
    chemo_drugs = ['Chlor', 'Chlor R', 'B Of', 'BR', 'CHOPR']
    sel_samples = [s for s in analysis.samples if s.treatment_active]
    labels = np.array([1 if s.treatment_type in chemo_drugs else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="chemo_treated")
    # B cell antibodies
    target_drugs = ['Alemtuz', 'Ibrutinib']
    sel_samples = [s for s in analysis.samples if s.treatment_active]
    labels = np.array([1 if s.treatment_type in target_drugs else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="target_treated")

    # Mutations / abnormalities
    muts = ['del13', 'del11', 'tri12']
    muts += ['SF3B1', 'ATM', 'NOTCH1', 'BIRC3', 'BCL2', 'TP53', 'MYD88', 'CHD2', 'NFKIE']
    # see mutations:
    # Counter([x for s in prj.samples if s.cellLine == "CLL" and s.technique == "ATAC-seq" for x in s.mutations])
    for mut in muts:
        # later add as well:
        # IGHVun +/- del; IGHVmut +/- del
        sel_samples = [s for s in analysis.samples if type(s.mutations) is list]
        labels = np.array([1 if mut in s.mutations else 0 for s in sel_samples])
        classify_samples(analysis, sel_samples, labels, comparison=mut)

    # Relapse
    # "relapse", ("True", "False"), # relapse or before relapse
    sel_samples = [s for s in analysis.samples if type(s.relapse) is bool]
    labels = np.array([1 if s.relapse else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="relapse")

    # Early vs Late
    # "progression", ("True", "False"), # early vs late timepoint
    sel_samples = [s for s in analysis.samples if type(s.diagnosis_collection) is bool]
    labels = np.array([1 if s.diagnosis_collection else 0 for s in sel_samples])
    classify_samples(analysis, sel_samples, labels, comparison="early_diagnosis")


def create_clinical_epigenomic_space(analysis):
    """
    """
    from sklearn.decomposition import PCA
    import cPickle as pickle
    import itertools
    from matplotlib.offsetbox import AnchoredText

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
            analysis.prj.dirs.plots, "trait_specific", "cll_peaks.medical_epigenomics_space.svg")
        fig.savefig(output_pdf, bbox_inches='tight')

    def sum_cartesian_vectors(vectors):
        """
        Given A = (xi, yi, zi, ...) and  B = (xj, yj, zj, ...);
        sum A + B = (xi + xj, yi + yj, zi + zj, ...)
        """
        return vectors.sum(axis=0)

    def plot_space(x, y, scale=1.0, axis=None):
        def pairwise(iterable):
            from itertools import tee, izip
            "s -> (s0,s1), (s1,s2), (s2, s3), ..."
            a, b = tee(iterable)
            next(b, None)
            return izip(a, b)

        # Get weighted variance to scale vectors
        # lam = (np.array(pca.explained_variance_[:2]) * np.sqrt(len(x))) ** scale
        # # divide observations by variance
        # xx = pd.DataFrame(x[0:2, :].T / lam)  # samples

        # or

        # Z-score variables to get them into same space
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)
        yy = y.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        # Plot samples and variables in same space (biplot)
        if axis is None:
            fig, axis = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)

        sample_colors = samples_to_color(analysis.samples, "unique")
        sample_symbols = samples_to_symbol(analysis.samples, "unique")

        var_colors = sns.color_palette("Paired", yy.shape[0])

        # Plot observations (samples)
        samples = pd.DataFrame([pd.Series(sample.__dict__) for sample in analysis.samples])
        for i in range(len(xx)):
            axis.scatter(xx.ix[i][0], xx.ix[i][1], s=50, color=sample_colors[i], marker=sample_symbols[i])
            # axis.annotate(
            #     "{0} - {1}".format(int(samples.ix[i]['patientID']), int(samples.ix[i]['timepoint'])),
            #     (xx.ix[i][0], xx.ix[i][1]),
            #     fontsize=8)  # add text with patient ID and timepoint

        # Plot variables (clinical features)
        for i, trait in enumerate(yy.index):
            axis.plot((0, yy.ix[trait][0]), (0, yy.ix[trait][1]), '-o', color=var_colors[i], label=trait)

        # add dashed line between patient's timepoints
        for patient, indexes in samples.groupby('patientID').groups.items():
            print(list(sorted(pairwise(samples.ix[indexes]["timepoint"]))))
            for t1, t2 in pairwise(samples.ix[indexes]["timepoint"]):
                tt1 = samples[(samples["timepoint"] == sorted([t1, t2])[0]) & (samples["patientID"] == patient)].index
                tt2 = samples[(samples["timepoint"] == sorted([t1, t2])[1]) & (samples["patientID"] == patient)].index
                # axis.plot((xx.ix[t1][0], xx.ix[t2][0]), (xx.ix[t1][1], xx.ix[t2][1]), "--", alpha=0.8, color="black")
                axis.annotate(
                    "",  # samples.ix[t1]["patientID"],
                    xy=(xx.ix[tt1][0], xx.ix[tt1][1]), xycoords='data',
                    xytext=(xx.ix[tt2][0], xx.ix[tt2][1]), textcoords='data',
                    arrowprops=dict(arrowstyle="fancy", color="0.5", shrinkB=5, connectionstyle="arc3,rad=0.3",),
                )
        axis.legend()

        if axis is None:
            output_pdf = os.path.join(
                analysis.prj.dirs.plots, "trait_specific", "cll_peaks.medical_epigenomics_space.biplot.%s.svg" % "-".join(y.index))
            fig.savefig(output_pdf, bbox_inches='tight')
        else:
            return axis

    #

    traits = ["IGHV", "CD38", "ZAP70", "primary_CLL", "treated", "chemo_treated", "target_treated"]
    muts = ['del11', 'tri12']  # abnormalities
    muts += ['TP53']  # mutations
    traits += muts

    # Here I am removing traits which had poor classification performance or too few patients associated
    # because they exist in either only one patient

    # Put trait-specific chromatin regions in one matrix
    features = pd.DataFrame()
    for trait in traits:
        file_name = os.path.join(
            analysis.prj.dirs.data,
            "trait_specific",
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
        try:
            df = pd.read_csv(file_name, sep="\t")
            df['trait'] = trait
            features = features.append(df, ignore_index=True)
        except IOError:
            print("Trait %s did not generate any associated regions" % trait)

    # # Save whole dataframe as csv
    # features.to_csv(os.path.join(analysis.prj.dirs.data, "trait_specific", "cll_peaks.medical_epigenomics_space.csv"), sep="\t", index=False)
    # features = pd.read_csv(os.path.join(analysis.prj.dirs.data, "trait_specific", "cll_peaks.medical_epigenomics_space.csv"), sep="\t")

    # Plot matrix of overlap between sets
    m = pd.DataFrame()  # build matrix of common features
    for t1, t2 in itertools.combinations(traits, 2):
        a = set(features[features['trait'] == t1][['start', 'end']].apply(sum, axis=1))
        b = set(features[features['trait'] == t2][['start', 'end']].apply(sum, axis=1))
        m.loc[t1, t2] = len(a.intersection(b)) / float(len(features[features['trait'] == t1].drop_duplicates()))
        m.loc[t2, t1] = len(a.intersection(b)) / float(len(features[features['trait'] == t2].drop_duplicates()))
    m.replace(np.nan, 1, inplace=True)
    sns.clustermap(np.log2(m.sort(axis=0).sort(axis=1)))
    output_pdf = os.path.join(
        analysis.prj.dirs.plots, "trait_specific", "cll_peaks.medical_epigenomics_space.common_trait_regions.svg")
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close('all')

    # # Investigate feature importance
    # # g = sns.FacetGrid(features[features['importance'] > 0.001], col="trait", col_wrap=5)
    # g = sns.FacetGrid(features, col="trait", col_wrap=5)
    # g.map(sns.distplot, "importance")

    # #
    # z_imp = pd.DataFrame()
    # for t in features.trait.unique():
    #     x = features[features['trait'] == t]['importance']
    #     z_imp = z_imp.append((x - x.mean()) / x.std())
    # [sns.distplot(z_imp.ix[i].dropna()) for i in range(len(z_imp))]

    # [sns.distplot(features[features['trait'] == t]['importance']) for t in traits]

    # # see fraction of features from total which are explained with clinical associations
    # features.shape[0]
    # features[['chrom', 'start', 'end']].drop_duplicates().shape[0]  # unique chromatin features
    fig, axis = plt.subplots(nrows=2, ncols=5, sharex=True, sharey=True, figsize=(60, 20))
    for i in range(1, len(traits) + 1):
        # PCA
        # get a matrix of unique features for all samples
        x = features[features['trait'].str.contains("|".join(traits[:i]))].drop_duplicates(['chrom', 'start', 'end'])[[s.name for s in analysis.samples] + ["trait"]]
        # save csvs for pca
        pca = PCA()
        x_new = pca.fit_transform(x.T.drop("trait"))

        # pickle.dump(x_new, open(os.path.join(analysis.prj.dirs.data, "trait_specific", "cll_peaks.medical_epigenomics_space.pca.pickle"), "w"))
        # x_new = pickle.load(open(os.path.join(analysis.prj.dirs.data, "trait_specific", "cll_peaks.medical_epigenomics_space.pca.pickle"), "r"))

        # # plot % explained variance per PC
        # fig, axis = plt.subplots(1)
        # axis.plot(range(1, len(pca.explained_variance_) + 1), pca.explained_variance_, 'o-')
        # axis.set_xlabel("PC")
        # axis.set_ylabel("% variance")
        # fig.savefig(os.path.join(analysis.prj.dirs.plots, "trait_specific", "cll_peaks.medical_epigenomics_space.pca_variance.svg"), bbox_inches='tight')

        # plot PCA
        # plot_pca(x_new)

        # Variable space
        # get one vector loading per feature (sum vectors of each chromatin feature)
        samples = features.columns[features.columns.str.contains("CLL")]
        loadings = pd.DataFrame(pca.components_, columns=x.index).T
        loadings['trait'] = x['trait']

        plot_space(x=pd.DataFrame(x_new), y=loadings.groupby('trait').sum(), axis= axis.flatten()[i - 1])

    # Save fig
    output_pdf = os.path.join(
        analysis.prj.dirs.plots, "trait_specific", "cll_peaks.medical_epigenomics_space.biplot.traits.svg")
    fig.savefig(output_pdf, bbox_inches='tight')

    # Get survival predictions


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

    # prj = Project("cll-patients")
    # prj.addSampleSheet("sequencing_sample_annotation_submission.csv")

    # Annotate with clinical data
    prj.samples = annotate_samples(prj.samples, clinical)

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
        'CLL_ATAC-seq_4621_1-5-36904_ATAC16-2_hg19']
    samples = [sample for sample in prj.samples if sample.technique == "ATAC-seq" and sample.cellLine == "CLL" and sample.name not in samples_to_exclude]

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
        analysis.estimate_peak_saturation(n=100)
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
        analysis.plot_coverage()
        analysis.plot_variance()
        analysis.plot_sample_correlations()
        # Observe exponential fit to the coeficient of variation
        analysis.plot_qv2_fit()

        # Characterize all CLL regions as a whole
        # region's composition
        characterize_regions_composition(df=analysis.coverage_qnorm_annotated, prefix="all_regions")
        # region's function
        output_dir = os.path.join(analysis.data_dir, "%s_peaks" % "all_regions")
        characterize_regions_function(df=analysis.coverage_qnorm_annotated, output_dir=output_dir, prefix="all_regions")
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
        analysis.coverage_qnorm = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.coverage_qnorm.log2.tsv"), sep="\t")
        analysis.coverage_qnorm_annotated = pd.read_csv(os.path.join(analysis.prj.dirs.data, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")

    # TRAIT-SPECIFIC ANALYSIS
    trait_analysis(analysis)
    # Clinical epigenomic space
    create_clinical_epigenomic_space(analysis)
    # Characterize patients in space

    # Add survival layer

    #

    #

    # Compare go terms/LOLA results from regions
    # compare_go_terms(
    #     [1, 2],
    #     ['data/mutation_peaks_cluster1/mutation_cluster1_regions.seq2pathway.csv', 'data/mutation_peaks_cluster2/mutation_cluster2_regions.seq2pathway.csv'],
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


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
