#!/usr/bin/env python

"""
This is the main script of the cll-patients project.
"""

# %logstart  # log ipython session

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['backend'] = "Agg"
# import recipy
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
try:  # stupid bug, importing it twice works
    from sklearn import cross_validation
except AttributeError:
    from sklearn import cross_validation
from sklearn.preprocessing import normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
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


def copy(obj):
    def copy(self):
        """
        Copy self to a new object.
        """
        from copy import deepcopy

        return deepcopy(self)
    obj.copy = copy
    return obj


@copy
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
    def normalize_coverage_quantiles(self, samples, samples2=None):
        # Normalize by quantiles
        to_norm = self.coverage[[s.name for s in samples]]
        self.coverage_qnorm = pd.DataFrame(
            normalize_quantiles_r(np.array(to_norm)),
            index=to_norm.index,
            columns=to_norm.columns
        )
        # Log2 transform
        self.coverage_qnorm = np.log2(1 + self.coverage_qnorm)

        if samples2 is not None:
            # Normalize by quantiles
            to_norm2 = self.coverage[[s.name for s in samples2]]
            qnorm2 = pd.DataFrame(
                normalize_quantiles_r(np.array(to_norm2)),
                index=to_norm2.index,
                columns=to_norm2.columns
            )
            # Log2 transform
            qnorm2 = np.log2(1 + qnorm2)

            # Join with other matrix
            self.coverage_qnorm = self.coverage_qnorm.join(qnorm2)

        # Add coordinates
        self.coverage_qnorm = self.coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        # Save
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
        # read list in
        sel_genes = pd.read_csv(os.path.join("metadata", "gene_lists", "bcell_cll_genelist.tsv"), sep="\t")
        sel_genes = OrderedDict(zip(sel_genes["gene_name"], sel_genes["ensembl_gene_id"]))

        # get distance to gene and ensembl gene id annotation in whole matrix
        df = pd.merge(self.coverage_qnorm_annotated, self.gene_annotation, on=['chrom', 'start', 'end', 'gene_name'])

        # GET 1(!) element per gene
        # get peaks around promoter (+/- 1kb)
        df2 = df[df['distance'] <= 2500]
        # promoters
        promoter_index = df2.groupby(["ensembl_gene_id"]).apply(lambda x: np.argmin((x['distance'])))
        promoters = df2.ix[promoter_index]

        # get peaks away from promoters (> 1kb)
        # enhancers = df[df['distance'] > 2500]
        # enhancers = df2.groupby(["ensembl_gene_id"])[["variance", "std_deviation", "dispersion", "qv2"]].apply(np.mean).reset_index()
        df2 = df[df['distance'] > 2500]
        # enhancers
        enhancer_index = df2.groupby(["ensembl_gene_id"]).apply(lambda x: np.argmin((x['distance'])))
        enhancers = df2.ix[enhancer_index]

        # Figure 2a - variability
        # 81 genes on top of all genes
        genes_str = "|".join(sel_genes.values())
        cll_p = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)]
        cll_e = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)]

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.distplot(np.log2(1 + promoters['std_deviation']), ax=axis[0], bins=100)
        sns.distplot(np.log2(1 + enhancers['std_deviation']), ax=axis[1], bins=100)
        sns.distplot(np.log2(1 + cll_p['std_deviation']), ax=axis[0], rug=True, bins=20)
        sns.distplot(np.log2(1 + cll_e['std_deviation']), ax=axis[1], rug=True, bins=20)
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.std_deviation.log2.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.distplot(np.log2(1 + promoters['variance']), ax=axis[0], bins=100)
        sns.distplot(np.log2(1 + enhancers['variance']), ax=axis[1], bins=100)
        sns.distplot(np.log2(1 + cll_p['variance']), ax=axis[0], rug=True, bins=20)
        sns.distplot(np.log2(1 + cll_e['variance']), ax=axis[1], rug=True, bins=20)
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.variance.log2.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.distplot(np.log2(1 + promoters['dispersion']), ax=axis[0], bins=100)
        sns.distplot(np.log2(1 + enhancers['dispersion']), ax=axis[1], bins=100)
        sns.distplot(np.log2(1 + cll_p['dispersion']), ax=axis[0], rug=True, bins=20)
        sns.distplot(np.log2(1 + cll_e['dispersion']), ax=axis[1], rug=True, bins=20)
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.dispersion.log2.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(2, sharex=True, sharey=True)
        sns.distplot(np.log2(1 + promoters['qv2']), ax=axis[0], bins=100)
        sns.distplot(np.log2(1 + enhancers['qv2']), ax=axis[1], bins=100)
        sns.distplot(np.log2(1 + cll_p['qv2']), ax=axis[0], rug=True, bins=20)
        sns.distplot(np.log2(1 + cll_e['qv2']), ax=axis[1], rug=True, bins=20)
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.qv2.log2.svg"), bbox_inches="tight")

        # test difference in distributions
        D, p = ks_2samp(promoters['std_deviation'], cll_p["std_deviation"])
        D, p = ks_2samp(enhancers['std_deviation'], cll_e["std_deviation"])
        D, p = ks_2samp(promoters['variance'], cll_p["variance"])
        D, p = ks_2samp(enhancers['variance'], cll_e["variance"])
        D, p = ks_2samp(promoters['dispersion'], cll_p["dispersion"])
        D, p = ks_2samp(enhancers['dispersion'], cll_e["dispersion"])
        D, p = ks_2samp(promoters['qv2'], cll_p["qv2"])
        D, p = ks_2samp(enhancers['qv2'], cll_e["qv2"])

        # random subset of genes
        n = 1000
        r_p_values = list()
        for i in range(n):
            r_p = promoters.ix[np.random.choice(promoters.index, len(cll_p), replace=False)]
            # r_e = enhancers.ix[np.random.choice(enhancers.index, len(cll_e), replace=False)]
            # plt.scatter(r_p["mean"], r_p["variance"], alpha=0.3)
            r_p_values.append(ks_2samp(r_p['variance'], cll_p["variance"])[1])

        # Number of elements
        # random subset of genes with similar number of reg elements
        elements_per_gene = df.groupby("ensembl_gene_id").apply(len)
        cll_elements_per_gene = elements_per_gene.ix[cll_p["ensembl_gene_id"]]

        fig, axis = plt.subplots(1)
        sns.distplot(np.log2(1 + elements_per_gene), ax=axis, bins=100)
        sns.distplot(np.log2(1 + cll_elements_per_gene), ax=axis, bins=20)
        fig.savefig(os.path.join(self.plots_dir, "genes.number_of elements.log2.svg"), bbox_inches="tight")
        # test
        D, p = ks_2samp(elements_per_gene, cll_elements_per_gene)

        # now chose one gene with the same number of elements
        n = 1000
        rm_p_values = list()
        for i in range(n):
            genes = list()
            for gene in cll_elements_per_gene:
                genes.append(np.random.choice(elements_per_gene[elements_per_gene == gene].index, 1, replace=False)[0])
            rm_p = promoters[promoters["ensembl_gene_id"].str.contains("|".join(genes))]
            rm_p_values.append(ks_2samp(rm_p['variance'], cll_p["variance"])[1])

        # random subset of genes with similar mean accessibility

        # housekeeping genes
        hk_genes = pd.read_csv(os.path.join("metadata", "gene_lists", "housekeeping_genes.tsv"), sep="\t")
        hk_genes = OrderedDict(zip(hk_genes["gene_name"], hk_genes["ensembl_gene_id"]))
        genes_str = "|".join(hk_genes.values())
        hk_p = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)]
        hk_e = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)]

        # test
        ks_2samp(promoters['variance'], cll_p["variance"])[1]
        ks_2samp(hk_p['variance'], cll_p["variance"])[1]
        ks_2samp(r_p['variance'], cll_p["variance"])[1]
        ks_2samp(rm_p['variance'], cll_p["variance"])[1]

        # plot
        plt.scatter(promoters["mean"], promoters["variance"], alpha=0.3, color="b")
        plt.scatter(hk_p["mean"], hk_p["variance"], alpha=1, color="g")
        plt.scatter(cll_p["mean"], cll_p["variance"], alpha=1, color="r")

        sns.jointplot(promoters["mean"], promoters["variance"])
        sns.jointplot(cll_p["mean"], cll_p["variance"])
        sns.jointplot(r_p["mean"], r_p["variance"])
        sns.jointplot(hk_p["mean"], hk_p["variance"])

        # Plot distributions of amplitude (fold_change)
        fig, axis = plt.subplots(1, figsize=(15, 10))
        sns.distplot(promoters['amplitude'], color="b", ax=axis)
        sns.distplot(enhancers['amplitude'], color="y", ax=axis)
        fig.savefig(os.path.join(self.plots_dir, "all_genes.distplot.svg"), bbox_inches='tight')

        # Inspect influence of CpG Islands
        # divide promoters and enhancers in CGI/non-CGI containing
        # intersect with CGI data
        cgi = pybedtools.BedTool("data/external/cpgIsland.hg19.bed")
        peaks_cgi = self.sites.intersect(cgi).to_dataframe()

        # split into CGI/non-CGI
        promoters_cgi = promoters.reset_index().merge(peaks_cgi, on=["chrom", "start", "end"], how="inner").set_index('index')
        promoters_non_cgi = promoters.ix[~promoters.index.isin(promoters_cgi.index)]
        enhancers_cgi = enhancers.reset_index().merge(peaks_cgi, on=["chrom", "start", "end"], how="inner").set_index('index')
        enhancers_non_cgi = enhancers.ix[~enhancers.index.isin(enhancers_cgi.index)]

        # plot
        fig, axis = plt.subplots(2)
        p = ks_2samp(promoters_non_cgi['variance'], promoters_cgi["variance"])[1]
        # sns.distplot(np.log2(1 + promoters["variance"]), bins=50, ax=axis[0])
        sns.distplot(np.log2(1 + promoters_non_cgi["variance"]), bins=50, ax=axis[0], label="non_CGI")
        sns.distplot(np.log2(1 + promoters_cgi["variance"]), bins=50, ax=axis[0], label="CGI\np-value: {:.2e}".format(p))
        # sns.distplot(np.log2(1 + enhancers["variance"]), bins=50, ax=axis[1])
        p = ks_2samp(enhancers_non_cgi['variance'], enhancers_cgi["variance"])[1]
        sns.distplot(np.log2(1 + enhancers_non_cgi["variance"]), bins=50, ax=axis[1], label="non_CGI")
        sns.distplot(np.log2(1 + enhancers_cgi["variance"]), bins=50, ax=axis[1], label="CGI\np-value: {:.2e}".format(p))
        axis[0].set_title("promoters")
        axis[1].set_title("distal elements")
        axis[0].set_xlabel("")
        axis[1].set_xlabel("log2(1 + variance)")
        axis[0].legend()
        axis[1].legend()
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.cgislands.variance.log2.svg"), bbox_inches="tight")

        # Now compare ALL vs CLL split into CGI/non-CGI
        genes_str = "|".join(sel_genes.values())
        cll_p = promoters[promoters['ensembl_gene_id'].str.contains(genes_str)]
        cll_e = enhancers[enhancers['ensembl_gene_id'].str.contains(genes_str)]

        cll_p_cgi = cll_p.reset_index().merge(peaks_cgi, on=["chrom", "start", "end"], how="inner").set_index('index')
        cll_p_non_cgi = cll_p.ix[~cll_p.index.isin(cll_p_cgi.index)]
        cll_e_cgi = cll_e.reset_index().merge(peaks_cgi, on=["chrom", "start", "end"], how="inner").set_index('index')
        cll_e_non_cgi = cll_e.ix[~cll_e.index.isin(cll_e_cgi.index)]

        # Test overrepresentation of CpG Islands
        from scipy.stats import fisher_exact
        fisher_exact([[cll_p_cgi.shape[0], promoters_cgi.shape[0]], [cll_p_non_cgi.shape[0], promoters_non_cgi.shape[0]]])

        # plot
        fig, axis = plt.subplots(2)
        sns.distplot(np.log2(1 + promoters_non_cgi["variance"]), bins=50, ax=axis[0], label="non_CGI")
        sns.distplot(np.log2(1 + promoters_cgi["variance"]), bins=50, ax=axis[0], label="CGI")
        sns.distplot(np.log2(1 + cll_p_non_cgi["variance"]), bins=20, ax=axis[0], label="non_CGI CLL")
        sns.distplot(np.log2(1 + cll_p_cgi["variance"]), bins=20, ax=axis[0], label="CGI CLL")

        sns.distplot(np.log2(1 + enhancers_non_cgi["variance"]), bins=50, ax=axis[1], label="non_CGI")
        sns.distplot(np.log2(1 + enhancers_cgi["variance"]), bins=50, ax=axis[1], label="CGI")
        sns.distplot(np.log2(1 + cll_e_non_cgi["variance"]), bins=20, ax=axis[1], label="non_CGI CLL")
        sns.distplot(np.log2(1 + cll_e_cgi["variance"]), bins=20, ax=axis[1], label="CGI CLL")
        axis[0].set_title("promoters")
        axis[1].set_title("distal elements")
        axis[0].set_xlabel("")
        axis[1].set_xlabel("log2(1 + variance)")
        axis[0].legend()
        axis[1].legend()
        fig.savefig(os.path.join(self.plots_dir, "prom-enh.cll-specfic.cgislands.variance.log2.svg"), bbox_inches="tight")

        #

        # Specific gene analysis
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
        fig.savefig(os.path.join(self.plots_dir, "relevant_genes.full.violinplot.svg"), bbox_inches='tight')

        # sort by predefined order (intensity/functional classes)
        sorterIndex = dict(zip(sel_genes.keys(), range(len(sel_genes.keys()))))
        boxplot_data['order'] = boxplot_data['gene'].map(sorterIndex)
        boxplot_data.sort('order', ascending=False, inplace=True)

        fig, axis = plt.subplots(1, figsize=(45, 10))
        sns.violinplot(data=boxplot_data, x="gene", y="openness", hue="region", split=True, inner="quart", palette={"promoter": "b", "enhancer": "y"}, jitter=True, ax=axis)
        fig.savefig(os.path.join(self.plots_dir, "relevant_genes.full.violinplot.funcorder.svg"), bbox_inches='tight')

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
        plt.hexbin(self.coverage_qnorm_annotated["mean"], self.coverage_qnorm_annotated["qv2"], bins='log')
        plt.savefig(os.path.join(self.plots_dir, "mean_qv2.svg"), bbox_inches="tight")
        # plot mean vs dispersion
        plt.hexbin(self.coverage_qnorm_annotated["mean"], self.coverage_qnorm_annotated["dispersion"], bins='log')
        plt.savefig(os.path.join(self.plots_dir, "mean_dispersion.svg"), bbox_inches="tight")

        # divide samples per IGHV status
        ighv_u = [s.name for s in samples if not s.ighv_mutation_status]
        ighv_m = [s.name for s in samples if s.ighv_mutation_status]

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

        # significantly variable regions
        # by threshold
        indexes = [i for i, p in p_values.items() if -np.log10(p) > 1.3]
        # filter out regions with mean across all samples lower than 1
        indexes = [i for i in indexes if self.coverage_qnorm_annotated.ix[i]["mean"] > 1]

        # plot uncorrected p_values
        sns.distplot(-np.log10(np.array(p_values.values())), hist=False)
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.p_values.svg"), bbox_inches="tight")
        plt.close("all")
        # volcano plot
        plt.scatter(np.log2(u['dispersion'] / m['dispersion']), -np.log10(np.array(p_values.values())))
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.volcano.svg"), bbox_inches="tight")
        plt.close("all")

        # correct p-values
        p_values = dict(zip(p_values.keys(), multipletests(p_values.values(), method="fdr_bh")[1]))
        # plot p-value distribution
        sns.distplot(-np.log10(np.array(p_values.values())), hist=False)
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.p_values-corrected.svg"), bbox_inches="tight")
        plt.close("all")
        # volcano plot
        plt.scatter(np.log2(u['dispersion'] / m['dispersion']), -np.log10(np.array(p_values.values())))
        plt.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.volcano-corrected.svg"), bbox_inches="tight")
        plt.close("all")

        # hexbin plots
        g = sns.FacetGrid(df, col="variable", sharex=False, sharey=False)
        for i, variable in enumerate(df["variable"].unique()):
            a = df[
                (df["variable"] == variable) &
                (df["subtype"] == "u")
            ]["value"].reset_index(drop=True)
            b = df[
                (df["variable"] == variable) &
                (df["subtype"] == "m")
            ]["value"].reset_index(drop=True)
            if i > 1:
                a = np.log2(1 + a)
                b = np.log2(1 + b)

            # distribution as hexbin
            g.axes[0][i].hexbin(a, b, bins="log")
            g.axes[0][i].set_title(variable)

            # get original indexes of significant regions
            aa = a.ix[indexes]
            bb = b.ix[indexes]

            # colors
            colors = (aa > bb).replace(True, "#d55e00").replace(False, "#0072b2")

            # significant as scatter
            g.axes[0][i].scatter(aa, bb, color=colors)

            # x=y line
            lims = [
                np.min([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # min of both axes
                np.max([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # max of both axes
            ]
            g.axes[0][i].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
            g.axes[0][i].set_aspect('equal')

        g.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.hexbin.significant.svg"), bbox_inches="tight")
        plt.close("all")

        # hexbin plots with color of highest variance to mean ratio ~= dispersion
        # get colors according to dispersion
        a = df[(df["variable"] == "dispersion") & (df["subtype"] == "u")]["value"].reset_index(drop=True).ix[indexes]
        b = df[(df["variable"] == "dispersion") & (df["subtype"] == "m")]["value"].reset_index(drop=True).ix[indexes]
        colors = (a > b).replace(True, "#d55e00").replace(False, "#0072b2")

        g = sns.FacetGrid(df, col="variable", sharex=False, sharey=False)
        for i, variable in enumerate(df["variable"].unique()):
            a = df[
                (df["variable"] == variable) &
                (df["subtype"] == "u")
            ]["value"].reset_index(drop=True)
            b = df[
                (df["variable"] == variable) &
                (df["subtype"] == "m")
            ]["value"].reset_index(drop=True)
            if i > 1:
                a = np.log2(1 + a)
                b = np.log2(1 + b)

            # distribution as hexbin
            g.axes[0][i].hexbin(a, b, bins="log")
            g.axes[0][i].set_title(variable)

            # get original indexes of significant regions
            aa = a.ix[indexes]
            bb = b.ix[indexes]

            # significant as scatter
            g.axes[0][i].scatter(aa, bb, color=colors)

            # x=y line
            lims = [
                np.min([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # min of both axes
                np.max([g.axes[0][i].get_xlim(), g.axes[0][i].get_ylim()]),  # max of both axes
            ]
            g.axes[0][i].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
            g.axes[0][i].set_aspect('equal')

        g.savefig(os.path.join(self.plots_dir, "ighv_mutation_variable_comparison.hexbin.significant.dispersion_colored.svg"), bbox_inches="tight")
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
        self.coverage_qnorm_annotated.ix[v_u][['chrom', 'start', 'end']].to_csv(
            os.path.join(self.data_dir, "ighv_mutation_variable_comparison", "ighv_mutation_variable_comparison.significant-uCLL.bed"),
            index=False, sep="\t", header=False)
        self.coverage_qnorm_annotated.ix[v_m][['chrom', 'start', 'end']].to_csv(
            os.path.join(self.data_dir, "ighv_mutation_variable_comparison", "ighv_mutation_variable_comparison.significant-mCLL.bed"),
            index=False, sep="\t", header=False)

        universe = os.path.join(self.data_dir, "cll_peaks.bed")

        # run lola
        lola(
            "ighv_mutation_variable_comparison.significant-uCLL.bed",
            universe,
            os.path.join(self.data_dir, "ighv_mutation_variable_comparison/uCLL"))
        lola(
            "ighv_mutation_variable_comparison.significant-mCLL.bed",
            universe,
            os.path.join(self.data_dir, "ighv_mutation_variable_comparison/mCLL"))

        # get gene names
        out = dict()
        out['u'] = self.coverage_qnorm_annotated.ix[v_u]['gene_name'].dropna().unique()
        out['m'] = self.coverage_qnorm_annotated.ix[v_m]['gene_name'].dropna().unique()
        # write gene names to file
        for c, genes in out.items():
            with open(os.path.join(self.data_dir, "ighv_mutation_variable_comparison/%sCLL.variable_genes.txt" % c), 'w') as handle:
                handle.writelines("\n".join(genes.tolist()))

    def correlate_expression(self):
        """
        Investigate degree of correlation between accessibility and gene expression using several approaches.
        """
        import itertools
        from scipy.stats import pearsonr

        def bin_distance(x):
            ranges = np.arange(0, 100000, 1000)
            return ranges.searchsorted(x)

        def hexbin(x, y, color, **kwargs):
            cmap = sns.light_palette(color, as_cmap=True)
            plt.hexbin(x, y, gridsize=15, cmap=cmap, **kwargs)

        # Get RNA
        # get samples with matched ATAC-seq and RNA-seq
        atac_ids = [s.sample_id for s in self.samples if s.library == "ATAC-seq"]
        rna_ids = [s.sample_id for s in self.samples if s.library == "RNA-seq"]
        atac_samples = [s for s in self.samples if s.library == "ATAC-seq"]
        rna_samples = [s for s in self.samples if s.library == "RNA-seq"]

        # get transcript/gene ids
        expression = pd.read_csv(
            os.path.join(rna_samples[0].paths.sample_root, "bowtie1_hg19_cdna", "bitSeq", rna_samples[0].name + ".tr"),
            sep=" ", skiprows=1, header=None)[[0, 1]]

        names = list()
        for sample in rna_samples:
            # if not os.path.exists(os.path.join(sample.paths.sample_root, "bowtie1_hg19_cdna", "bitSeq", sample.name + ".counts")):
            #     print sample.name
            #     continue
            # read in counts
            expr = pd.read_csv(
                os.path.join(sample.paths.sample_root, "bowtie1_hg19_cdna", "bitSeq", sample.name + ".counts"),
                skiprows=1, header=None, names=[sample.name])
            # append
            expression = pd.concat([expression, expr], ignore_index=True, axis=1)
            names.append(sample.name)
        expression.columns = ["ensembl_gene_id", "ensembl_transcript_id"] + names

        # get gene-level quantification by having transcript with highest value
        expression_genes = expression.groupby("ensembl_gene_id").apply(max)

        # save raw counts
        expression_genes.to_csv(os.path.join(self.data_dir, "cll_expression_matrix.csv"))

        # log expression
        expression_genes = np.log2(expression_genes.drop(["ensembl_gene_id", "ensembl_transcript_id"], axis=1))

        # average across all samples
        expression_genes["median"] = expression_genes.apply(np.median, axis=1)

        # save expression matrix
        expression_genes.to_csv(os.path.join(self.data_dir, "cll_expression_matrix.log2.csv"))
        # expression_genes = pd.read_csv(os.path.join(self.data_dir, "cll_expression_matrix.log2.csv")).set_index("ensembl_gene_id")

        # Explore expression
        fig, axis = plt.subplots(1)
        sns.distplot(expression_genes['median'], bins=100, ax=axis)
        fig.savefig(os.path.join(self.plots_dir, "expression_all_samples.median.svg"), bbox_inches="tight")

        # a few samples vs samples
        fig, axis = plt.subplots(2, 2)
        axis = axis.flatten()
        for i in range(4):
            axis[i].scatter(
                expression_genes[expression_genes.columns[i]],
                expression_genes[expression_genes.columns[i + 1]])
        fig.savefig(os.path.join(self.plots_dir, "expression_samples_vs_samples.median.svg"), bbox_inches="tight")

        # By IGHV group
        fig, axis = plt.subplots(3, figsize=(10, 30))
        for i, (group1, group2) in enumerate(itertools.combinations(["uCLL", "iCLL", "mCLL"], 2)):
            m1 = expression_genes[[s.name for s in rna_samples if (s.ighv_group == group1) and (s.name in names)]]
            m2 = expression_genes[[s.name for s in rna_samples if (s.ighv_group == group2) and (s.name in names)]]
            axis[i].scatter(m1.apply(np.median, axis=1), m2.apply(np.median, axis=1))
            axis[i].set_xlabel(group1 + " - %i samples" % m1.shape[1])
            axis[i].set_ylabel(group2 + " - %i samples" % m2.shape[1])
        fig.savefig(os.path.join(self.plots_dir, "expression_ighv-groups.median.svg"), bbox_inches="tight")

        # Get openness
        openness = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated["chrom"].str.contains("chr[^X|Y]")]
        # get closest gene info with ensembl ids
        g = pd.read_csv(os.path.join(self.data_dir, "cll_peaks.gene_annotation.csv"))
        openness = pd.merge(openness, g)
        # add ensembl id as index
        openness.index = openness["ensembl_gene_id"]

        # Make dataframe of matched accessibility and expression per patient sample including distance
        matched = pd.DataFrame()
        for sample_id in set(atac_ids).intersection(set(rna_ids)):
            a = openness[[s.name for s in atac_samples if s.sample_id == sample_id] + ["distance"]]
            a.columns = ["atac", "distance"]
            r = expression_genes[[s.name for s in rna_samples if (s.sample_id == sample_id) and (s.name in names)]]
            if r.shape[1] != 1:
                continue
            r.columns = ["rna"]

            # join both
            df = a.join(r).dropna()
            df["sample_id"] = sample_id
            matched = matched.append(df)
        # Put gene-element pairs in bins dependent on distance
        ranges = np.arange(0, 100000, 1000)
        matched["distance_bin"] = ranges.searchsorted(matched["distance"])

        # save
        matched.drop("ensembl_gene_id", axis=1).to_csv(os.path.join(self.data_dir, "cll_expression_accessibility_matrix.csv"))
        matched = pd.read_csv(os.path.join(self.data_dir, "cll_expression_accessibility_matrix.csv"))

        # filtering
        matched[
            (matched["atac"] > 2) &
            (matched["rna"] > 2.5)
        ]

        # Brute force
        # all genes with associated elements vs median of all patients
        x = matched[['atac', 'sample_id']].groupby(level=0)['atac'].apply(np.median)
        y = matched[['rna', 'sample_id']].groupby(level=0)['rna'].apply(np.median)
        j = sns.jointplot(x, y, kind="hex", color="#4CB391", joint_kws={'bins': 'log'})
        sns.regplot(x, y, ax=j.ax_joint, scatter=False)
        j.fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.brute_force_median.hex.svg"), bbox_inches="tight")

        # Brute force, patient-specific
        # all genes with associated elements within each matched patient
        j = sns.jointplot(matched["atac"], matched["rna"], kind="hex", color="#4CB391", joint_kws={'bins': 'log'})
        sns.regplot(matched["atac"], matched["rna"], ax=j.ax_joint, scatter=False)
        j.fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.hex.svg"), bbox_inches="tight")

        # plotted individually
        g = sns.FacetGrid(matched, col="sample_id", col_wrap=3)
        g.map(plt.hexbin, "atac", "rna", bins='log', gridsize=50, marginals=False, edgecolors=None)
        g.fig.savefig(os.path.join(self.plots_dir, "norm_counts.mean.per_genomic_region.hex.svg"), bbox_inches="tight")
        plt.close()

        # Promoters only, patient-specific
        # only promoters with associated elements within each matched patient
        # get only genes within 5kb
        proms = matched[matched["distance"] < 2500]
        j = sns.jointplot(proms["atac"], proms["rna"], kind="hex", color="#4CB391", joint_kws={'bins': 'log'})
        sns.regplot(proms["atac"], proms["rna"], ax=j.ax_joint, scatter=False)
        j.fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.promoter.hex.svg"), bbox_inches="tight")

        # Distal elemnts only, patient-specific
        # only promoters with associated elements within each matched patient
        # get only genes within 5kb
        distal = matched[matched["distance"] > 2500]
        j = sns.jointplot(distal["atac"], distal["rna"], kind="hex", color="#4CB391", joint_kws={'bins': 'log'})
        sns.regplot(distal["atac"], distal["rna"], ax=j.ax_joint, scatter=False)
        j.fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.distal.hex.svg"), bbox_inches="tight")

        # Promoters only vs distal elements only, patient-specific
        # only promoters and only distal elements with associated elements within each matched patient
        fig, axis = plt.subplots(2)
        axis[0].hexbin(proms["atac"], proms["rna"], color="#4CB391", bins='log')
        sns.regplot(distal["atac"], distal["rna"], ax=axis[0], scatter=False)
        axis[1].hexbin(distal["atac"], distal["rna"], color="#4CB391", bins='log')
        sns.regplot(distal["atac"], distal["rna"], ax=axis[1], scatter=False)
        fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.promoter_vs_distal.hex.svg"), bbox_inches="tight")

        # Promoters only vs distal elements only, patient-specific, distance dependent
        # only promoters and only distal elements with associated elements within each matched patient
        # ranges = dict(zip(range(100), [(i * 1000, (i * 1000) + 1000) for i in range(100)]))
        cor = dict()
        cor_filtered = dict()
        for i in range(100):
            subset = matched[
                matched["distance_bin"] == i
            ]
            # add correlation
            cor[i] = (pearsonr(subset["atac"], subset["rna"]), subset.shape[0])

            subset = matched[
                (matched["distance_bin"] == i) &
                (matched["rna"] > 2)
            ]
            # add correlation
            cor_filtered[i] = (pearsonr(subset["atac"], subset["rna"]), subset.shape[0])

        # plot correlation vs distance
        fig, axis = plt.subplots(2)
        axis[0].plot(cor.keys(), [k[0][0] for k in cor.values()])
        axis[1].plot(cor.keys(), -np.log10([k[0][1] for k in cor.values()]))
        axis[1].twinx().plot(cor_filtered.keys(), np.log10([k[1] for k in cor.values()]), color="green")
        axis[0].set_title("ATAC-seq/RNA-seq concordance")
        axis[0].set_ylabel("Pearson correlation")
        axis[1].set_ylabel("correlation p-value (-log10)")
        axis[1].set_xlabel("distance to TSS (kb)")
        fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.correlation_distance.svg"), bbox_inches="tight")

        # plot correlation vs distance
        fig, axis = plt.subplots(2)
        axis[0].plot(cor_filtered.keys(), [k[0][0] for k in cor_filtered.values()])
        axis[1].plot(cor_filtered.keys(), -np.log10([k[0][1] for k in cor_filtered.values()]))
        axis[1].twinx().plot(cor_filtered.keys(), np.log10([k[1] for k in cor_filtered.values()]), color="green")
        axis[0].set_title("ATAC-seq/RNA-seq concordance")
        axis[0].set_ylabel("Pearson correlation")
        axis[1].set_ylabel("correlation p-value (-log10)")
        axis[1].set_xlabel("distance to TSS (kb)")
        fig.savefig(os.path.join(self.plots_dir, "expression_oppenness_correlation.patient_matched.expressed_genes.correlation_distance.svg"), bbox_inches="tight")

        """
        IGHV differential expression relation to accessibility.
        """
        from scipy.stats import pearsonr

        samples = [sample for sample in self.samples if sample.library == "ATAC-seq"]

        # IGHV dependent colors
        group_colors = samples_to_color(samples, "IGHV")
        homology_colors = samples_to_color(samples, "ighv_homology")
        # colorbar
        cmap = plt.get_cmap('summer')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=86, vmax=100))
        sm._A = []

        # set genes as index
        matched = matched.set_index("ensembl_gene_id")
        openness = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated["chrom"].str.contains("chr[^X|Y]")]
        # get closest gene info with ensembl ids
        g = pd.read_csv(os.path.join(self.data_dir, "cll_peaks.gene_annotation.csv"))
        openness = pd.merge(openness, g)

        group = "{}_{}".format("uCLL", "mCLL")
        # Get differentially expressed genes
        diff_expr = pd.read_csv(os.path.join(self.data_dir, "gene_expression.differential.{}-{}_results.degs.csv".format("uCLL", "mCLL")), index_col=0)

        # # Separate genes in up/down regulated
        up = diff_expr[diff_expr["log2FoldChange"] > 0]
        down = diff_expr[diff_expr["log2FoldChange"] < 0]
        # Get genes of each group
        up = openness[openness["ensembl_gene_id"].isin(up.index)]
        down = openness[openness["ensembl_gene_id"].isin(down.index)]
        # Get samples of each group
        up = up[[s.name for s in samples]].mean()
        down = down[[s.name for s in samples]].mean()

        # plot mean vs rank
        fig, axis = plt.subplots(2, 2, figsize=(10, 8))
        # Plot Rank vs mean
        axis[0][0].scatter(up.rank(ascending=False), up, label=group, color=group_colors)
        axis[0][1].scatter(up.rank(ascending=False), up, label=group, color=homology_colors)
        axis[1][0].scatter(down.rank(ascending=False), down, label=group, color=group_colors)
        axis[1][1].scatter(down.rank(ascending=False), down, label=group, color=homology_colors)
        plt.legend()
        axis[0][0].set_ylabel("Mean accessibility")
        axis[1][0].set_ylabel("Mean accessibility")
        axis[0][0].set_title("uCLL differential regions")
        axis[0][1].set_title("mCLL differential regions")
        axis[1][0].set_xlabel("Sample rank in mean accessibility")
        axis[1][1].set_xlabel("Sample rank in mean accessibility")
        fig.savefig(os.path.join(self.plots_dir, "gene_expression.differential_regions.accessibility_rank.svg"), bbox_inches="tight")

        # plot mean uCLL vs mean mCLL
        fig2, axis2 = plt.subplots(2, 1, figsize=(10, 8))
        r, p = pearsonr(up, down)
        # sns.regplot(up, down, scatter=False, fit_reg=True, ax=axis2[0])
        # sns.regplot(up, down, scatter=False, fit_reg=True, ax=axis2[1])
        axis2[0].scatter(up, down, label="r={}\np={}".format(r, p), color=group_colors)
        axis2[1].scatter(up, down, label="r={}\np={}".format(r, p), color=homology_colors)
        plt.legend()
        plt.colorbar(sm, orientation="horizontal", aspect=5.)
        axis2[0].set_xlabel("Mean accessibility in uCLL differential regions")
        axis2[0].set_ylabel("Mean accessibility in mCLL differential regions")
        axis2[1].set_xlabel("Mean accessibility in uCLL differential regions")
        axis2[1].set_ylabel("Mean accessibility in mCLL differential regions")
        fig2.savefig(os.path.join(self.plots_dir, "gene_expression.differential_regions.accessibility.svg"), bbox_inches="tight")

    def correlate_expression_spanish_cohort(self):
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

    def correlate_methylation(self, samples):
        """
        See how accessibility varies in windows with TF motifs
        shown to have different methylation levels along B cell development
        (and IGHV status of CLL)
        """
        from scipy.stats import pearsonr

        # IGHV dependent colors
        group_colors = samples_to_color(samples, "IGHV")
        homology_colors = samples_to_color(samples, "ighv_homology")
        # colorbar
        cmap = plt.get_cmap('summer')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=86, vmax=100))
        sm._A = []

        # Kulis 2012 dataset
        acc_meth = pd.DataFrame()

        # plot mean vs rank
        fig, axis = plt.subplots(2, 2, figsize=(10, 8))
        for i, group in enumerate(["U-M", "M-U"]):
            # Intersect TF binding site windows with CLL peaks
            diff_meth = pybedtools.BedTool(os.path.join("data", "external", "Kulis_DNAme.%s.hg19.bed" % group))
            overlap = self.sites.intersect(diff_meth, wa=True).to_dataframe()
            overlap.columns = ["chrom", "start", "end"]
            peaks_overlap = pd.merge(self.coverage_qnorm_annotated, overlap, how="right", on=["chrom", "start", "end"])
            mean_tf = peaks_overlap[[s.name for s in samples]].mean()
            acc_meth[group] = mean_tf
            # sns.regplot(mean_tf.rank(ascending=False), mean_tf, label=group, fit_reg=False)
            axis[0][i].scatter(mean_tf.rank(ascending=False), mean_tf, label=group, color=group_colors)
            axis[1][i].scatter(mean_tf.rank(ascending=False), mean_tf, label=group, color=homology_colors)
        axis[0][0].set_ylabel("Mean accessibility")
        axis[1][0].set_ylabel("Mean accessibility")
        axis[0][0].set_title("uCLL hypermethylated regions")
        axis[0][1].set_title("mCLL hypermethylated regions")
        axis[1][0].set_xlabel("Sample rank in mean accessibility")
        axis[1][1].set_xlabel("Sample rank in mean accessibility")
        fig.savefig(os.path.join(self.plots_dir, "hypermethylated_regions.accessibility_rank.svg"), bbox_inches="tight")

        # plot mean uCLL vs mean mCLL
        r, p = pearsonr(acc_meth["U-M"], acc_meth["M-U"])
        fig, axis = plt.subplots(2, figsize=(10, 15))
        axis[0].scatter(acc_meth["U-M"], acc_meth["M-U"], color=group_colors, label="r={}\np={}".format(r, p))
        axis[1].scatter(acc_meth["U-M"], acc_meth["M-U"], color=homology_colors, label="r={}\np={}".format(r, p))
        plt.legend()
        plt.colorbar(sm, orientation="horizontal", aspect=5.)
        axis[0].set_xlabel("Mean accessibility in uCLL hypermethylated regions")
        axis[0].set_ylabel("Mean accessibility in mCLL hypermethylated regions")
        axis[1].set_xlabel("Mean accessibility in uCLL hypermethylated regions")
        axis[1].set_ylabel("Mean accessibility in mCLL hypermethylated regions")
        fig.savefig(os.path.join(self.plots_dir, "hypermethylated_regions.accessibility.svg"), bbox_inches="tight")

        # Oakes 2016 dataset
        # Intersect TF binding site windows with CLL peaks
        all_me_tfbs = pybedtools.BedTool(os.path.join("data", "external", "TF_DNAme_windows.hg19.bed"))
        overlap = self.sites.intersect(all_me_tfbs, wa=True).to_dataframe()
        overlap.columns = ["chrom", "start", "end"]
        peaks_overlap = pd.merge(self.coverage_qnorm_annotated, overlap, how="right", on=["chrom", "start", "end"])

        # Get mean accessibility per sample in those peaks
        mean_acc = peaks_overlap[[s.name for s in samples]].mean()

        # Plot rank vs mean accessibility
        fig, axis = plt.subplots(1, figsize=(8, 8))
        sns.regplot(mean_acc.rank(ascending=False), mean_acc, label="all", fit_reg=False, color="black", ax=axis)

        # Split by TF
        tfs = ["AP-1", "EBF1", "RUNX3", "IRF4", "OCT2", "NFkB"]
        for tf in tfs:
            tfbs = pybedtools.BedTool(os.path.join("data", "external", "TF_DNAme_windows.%s.hg19.bed" % tf))
            overlap = self.sites.intersect(tfbs, wa=True).to_dataframe()
            overlap.columns = ["chrom", "start", "end"]
            tfbs_overlap = pd.merge(self.coverage_qnorm_annotated, overlap, how="right", on=["chrom", "start", "end"])
            mean_tf = tfbs_overlap[[s.name for s in samples]].mean()
            sns.regplot(mean_tf.rank(ascending=False), mean_tf, label=tf, fit_reg=False, ax=axis)
        axis.legend()
        fig.savefig(os.path.join(self.plots_dir, "methylated_regions_TFBS.accessibility_rank.svg"), bbox_inches="tight")

        # Plot rank vs mean accessibility
        fig, axis = plt.subplots(4, 2, figsize=(6, 10))
        axis = axis.flatten()
        axis[0].scatter(mean_acc.rank(ascending=False), mean_acc, label="all", color=homology_colors)
        axis[0].set_title("all")
        axis[0].set_xlim((-5, 93))

        # Split by TF
        tfs = ["AP-1", "EBF1", "RUNX3", "IRF4", "OCT2", "NFkB"]
        for i, tf in enumerate(tfs):
            tfbs = pybedtools.BedTool(os.path.join("data", "external", "TF_DNAme_windows.%s.hg19.bed" % tf))
            overlap = self.sites.intersect(tfbs, wa=True).to_dataframe()
            overlap.columns = ["chrom", "start", "end"]
            tfbs_overlap = pd.merge(self.coverage_qnorm_annotated, overlap, how="right", on=["chrom", "start", "end"])
            mean_tf = tfbs_overlap[[s.name for s in samples]].mean()
            axis[i + 1].scatter(mean_tf.rank(ascending=False), mean_tf, label=tf, color=homology_colors)
            axis[i + 1].set_title(tf)
            axis[i + 1].set_xlim((-5, 93))
        axis[i + 1].set_xlabel("sample rank")
        axis[i + 2].set_xlabel("sample rank")
        axis[i + 2].set_visible(False)
        [axis[i].set_ylabel("mean accessibility") for i in range(0, 8, 2)]
        plt.colorbar(sm, orientation="horizontal", aspect=5.)
        fig.savefig(os.path.join(self.plots_dir, "methylated_regions_TFBS.accessibility_rank.svg"), bbox_inches="tight")

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

        # this is loaded now
        df = pd.read_csv(os.path.join(self.data_dir, "cll_peaks.support.csv"))
        sns.distplot(df["support"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.support.distplot.svg"), bbox_inches="tight")
        plt.close()

    def plot_coverage(self):
        data = self.coverage_qnorm_annotated.copy()

        # distribution of count attributes
        sns.distplot(data["mean"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.mean.distplot.svg"), bbox_inches="tight")
        plt.close()
        sns.distplot(data["qv2"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.qv2.distplot.svg"), bbox_inches="tight")
        plt.close()
        sns.distplot(data["dispersion"], rug=False)
        plt.savefig(os.path.join(self.plots_dir, "cll_peaks.dispersion.distplot.svg"), bbox_inches="tight")
        plt.close()

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


def normalize_variation_r(array):
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("limma")')
    normq = robjects.r('normalizeVSN')
    return np.array(normq(array + 1))


def name_to_repr(name):
    return "_".join([name.split("_")[0]] + [name.split("_")[2]] + name.split("_")[3:4])


def name_to_id(name):
    """This returns joined patient and sample IDs"""
    return "_".join([name.split("_")[2]] + name.split("_")[3:4])


def name_to_patient_id(name):
    return name.split("_")[2]


def name_to_sample_id(name):
    return name.split("_")[3]


def samples_to_color(samples, trait="IGHV"):
    # unique color per patient
    if trait == "patient":
        patients = set([sample.patient_id for sample in samples])
        color_dict = cm.Paired(np.linspace(0, 1, len(patients)))
        color_dict = dict(zip(patients, color_dict))
        return [color_dict[sample.patient_id] for sample in samples]
    # rainbow (unique color per sample)
    elif trait == "unique_sample":
        return cm.Paired(np.linspace(0, 1, len(samples)))
        # gender
    elif trait == "gender":
        colors = list()
        for sample in samples:
            if sample.patient_gender == "F":
                colors.append('red')
            elif sample.patient_gender == "M":
                colors.append('blue')
            else:
                colors.append('gray')
        return colors
    # disease at diagnosis time
    elif trait == "disease":
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
    # dependent on trait threshold
    elif trait in ["IGHV"]:
        # This uses sns colorblind pallete
        colors = list()
        for sample in samples:
            if getattr(sample, trait) == 1:
                colors.append(sns.color_palette("colorblind")[0])  # blue #0072b2
            elif getattr(sample, trait) == 0:
                colors.append(sns.color_palette("colorblind")[2])  # vermillion #d55e00
            else:
                colors.append('gray')
        return colors
    # IGHV homology color scale from min to max
    if trait in ["ighv_homology"]:
        vmin = min([getattr(s, trait) for s in samples])
        # This uses sns summer colormap
        cmap = plt.get_cmap('summer')
        # scale colormap to min and max ighv homology
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=100)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        # get colors acordingly
        colors = list()
        for sample in samples:
            if vmin <= getattr(sample, trait) <= 100:
                colors.append(m.to_rgba(getattr(sample, trait)))
            else:
                colors.append('gray')
        return colors
    else:
        raise ValueError("trait %s is not valid" % trait)


def all_sample_colors(samples, order=""):
    return [
        samples_to_color(samples, "patient"),
        samples_to_color(samples, "gender"),
        samples_to_color(samples, "disease"),
        samples_to_color(samples, "IGHV"),
        samples_to_color(samples, "ighv_homology")
    ]


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
            "cll_peaks.%s_significant.classification.random_forest.loocv.ROC_PRC.svg" % trait), bbox_inches="tight")

        # Display training and prediction of pre-labeled samples of most informative features:
        # average feature importance across iterations
        mean_importance = importance.mean(axis=0)

        # visualize feature importance
        fig, axis = plt.subplots(1)
        sns.distplot(mean_importance, ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.mean_importance.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # SEE VALUES OF ALL SAMPLES IN IMPORTANT FEATURES
        # Display most informative features for ALL samples:
        matrix = analysis.coverage_qnorm_annotated[[s.name for s in all_samples]]
        # get important features
        x = matrix.loc[[i for i, j in enumerate(mean_importance > 1e-4) if j == True], :]  # get features on the tail of the importance distribution

        # Add info
        dataframe = analysis.coverage_qnorm_annotated.loc[x.index, :]
        # add difference of standardized openness between positive and negative groups used in classification
        df2 = dataframe[[s.name for s in sel_samples]].apply(lambda j: (j - j.min()) / (j.max() - j.min()), axis=0)
        dataframe["change"] = df2.icol([i for i, l in enumerate(labels) if l == 1]).mean(axis=1) - df2.icol([i for i, l in enumerate(labels) if l == 0]).mean(axis=1)
        # add direction of chromatin feature association with trait
        dataframe['direction'] = dataframe['change'].apply(lambda x: 1 if x > 0 else -1)
        # add feature importance
        dataframe["importance"] = mean_importance[x.index]

        # Save whole dataframe as csv
        dataframe_file = os.path.join(
            analysis.data_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
        dataframe.to_csv(dataframe_file, index=False)

        # Save as bed
        bed_file = os.path.join(
            analysis.data_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.sites.bed" % trait)
        dataframe[["chrom", "start", "end"]].to_csv(bed_file, sep="\t", header=False, index=False)

        # sample correlation dendrogram
        dend = sns.clustermap(
            x.corr(),
            col_colors=all_samples_colors,  # all_sample_colors(all_samples),
            row_colors=all_samples_colors)
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sample_correlation.svg" % trait), bbox_inches="tight")
        plt.close("all")

        # extract linkage matrix
        from scipy.cluster.hierarchy import fcluster
        lm = dend.dendrogram_col.linkage

        mclust = 4  # this is determined with supervision
        clusters = fcluster(lm, mclust, criterion='maxclust')
        col_dict = dict(zip(np.unique(clusters), sns.color_palette("colorblind")[1:mclust + 1]))
        # get colors according to cluster
        cluster_colors = [col_dict[c] for c in clusters]

        # save sample assignment
        assignments = pd.DataFrame([[s.name for s in all_samples], list(clusters)]).T
        assignments.columns = ["sample", "cluster"]
        assignments.to_csv(os.path.join(
            analysis.data_dir, "cll_peaks.%s_significant.classification.random_forest.loocv.pca.sample_assignments.csv" % trait), index=False)

        # pca on these regions
        pca_r(x, cluster_colors, os.path.join(
            analysis.plots_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.pca.sample_labels.svg" % trait))

        # region-sample heatmap with colors of the direction each region is associated to
        region_colors = dict(zip([1, -1], sns.color_palette("colorblind")[1:3]))
        colors = [region_colors[c] for c in dataframe['direction']]
        sns.clustermap(
            x,
            cmap=cmap,
            standard_scale=0,
            col_colors=cluster_colors,
            row_colors=colors,
            yticklabels=False)
        plt.savefig(os.path.join(
            analysis.plots_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.clustering_sites.sites_labels.svg" % trait), bbox_inches="tight")
        plt.close("all")


def classification_random(analysis, sel_samples, labels, trait, n=100):
    """
    Use a machine learning approach for sample classification based on known sample attributes.
    Extract features most important to separate samples and investigate those.
    """
    print("Trait:%s" % trait)
    print("%i samples with trait annotated" % len(sel_samples))
    print(Counter(labels))

    # BINARY CLASSIFICATION ON ALL CLL OPEN CHROMATIN REGIONS
    # get features and labels
    X = normalize(analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]].T)

    fig, axis = plt.subplots(1, 2, figsize=(12, 5))
    metrics = pd.DataFrame(index=["fpr", "tpr", "roc_auc", "TNR", "TPR", "recall", "precision", "aps"])
    for i in range(n):
        # randomize labels
        y = np.array(np.random.permutation(labels))

        loo = cross_validation.LeaveOneOut(len(X))

        for j, (train_index, test_index) in enumerate(loo):
            print(i, j)
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

            if j == 0:
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
        tn = len([1 for k in range(len(binary_scores)) if (binary_labels[k] == 1) and (binary_scores[k] == 1)])
        n = len([1 for k in range(len(binary_scores)) if binary_labels[k] == 1])
        TNR = tn / float(n)
        # Sensitivity (TP / P)
        tp = len([1 for k in range(len(binary_scores)) if (binary_labels[k] == 0) and (binary_scores[k] == 0)])
        p = len([1 for k in range(len(binary_scores)) if binary_labels[k] == 0])
        TPR = tp / float(p)
        # FPR (FP / P)
        fn = len([1 for k in range(len(binary_scores)) if (binary_labels[k] == 0) and (binary_scores[k] == 1)])
        p = len([1 for k in range(len(binary_scores)) if binary_labels[k] == 0])
        FPR = fn / float(p)
        # FNR (FN / P)
        fn = len([1 for k in range(len(binary_scores)) if (binary_labels[k] == 1) and (binary_scores[k] == 0)])
        p = len([1 for k in range(len(binary_scores)) if binary_labels[k] == 1])
        FNR = fn / float(p)

        # Compute ROC curve and ROC area for each class
        fpr, tpr, _ = roc_curve(y_all_test, y_all_scores[:, 1], pos_label=1)
        roc_auc = auc(fpr, tpr, reorder=True)
        # Compute Precision-Recall and average precision
        precision, recall, _ = precision_recall_curve(y_all_test, y_all_scores[:, 1], pos_label=1)
        binary_labels = [0 if x == classifier.classes_[0] else 1 for x in y_all_test]
        aps = average_precision_score(binary_labels, y_all_scores[:, 1])

        # Plot ROC and PRC curves
        axis[0].plot(fpr, tpr, alpha=0.3)
        axis[1].plot(recall, precision, alpha=0.3)
        metrics[i] = [fpr, tpr, roc_auc, TNR, TPR, recall, precision, aps]
        metrics.T.to_csv(os.path.join(
            analysis.plots_dir,
            "cll_peaks.%s-random_significant.classification.random_forest.loocv.metrics.csv" % trait), index=False)

    # calculate mean metrics
    # values
    auc_m = np.median([x[2] for x in metrics])
    TNR_m = np.median([x[3] for x in metrics])
    TPR_m = np.median([x[4] for x in metrics])
    APS_m = np.median([x[-1] for x in metrics])
    # rates
    # here we are going to interpolate values between points in the ROC curves
    # across all iterations in a range (e.g. 0-100)
    fpr_m = list()
    tpr_m = list()
    recall_m = list()
    precision_m = list()
    for i in np.array(range(101)) / 100.:
        fpr_m.append(np.median([np.interp(i, x[0], np.array(range(len(x[0]))) / float(len(x[0]))) for x in metrics]))
        tpr_m.append(np.median([np.interp(i, x[1], np.array(range(len(x[1]))) / float(len(x[1]))) for x in metrics]))
        recall_m.append(np.median([np.interp(i, x[5], np.array(range(len(x[5]))) / float(len(x[5]))) for x in metrics]))
        precision_m.append(np.median([np.interp(i, x[6], np.array(range(len(x[6]))) / float(len(x[6]))) for x in metrics]))

    fig, axis = plt.subplots(1, 2, figsize=(12, 5))

    # plot mean curves
    axis[0].plot(fpr_m, tpr_m, label='ROC (median AUC = {0:0.2f}; median TNR = {1:0.2f}; median TPR = {2:0.2f})'.format(auc_m, TNR_m, TPR_m), color="black")
    axis[1].plot(recall_m, precision_m, label='PRC (median AUC = {0:0.2f})'.format(APS_m), color="black")

    # other stylistics
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
    fig.savefig(os.path.join(
        analysis.plots_dir,
        "cll_peaks.%s-random_significant.classification.random_forest.loocv.ROC_PRC.svg" % trait), bbox_inches="tight")


def unsupervised(analysis, samples):
    """
    Run trait classification (with independent validation if possible for that trait)
    on all samples with known annotation of the trait.
    """
    from sklearn.decomposition import PCA

    # Approach 1: all regions
    X = analysis.coverage_qnorm_annotated[[s.name for s in samples]]

    pca = PCA()
    x_new = pca.fit_transform(X.T)
    # transform again
    x = pd.DataFrame(x_new)
    xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

    # plot % explained variance per PC
    fig, axis = plt.subplots(1)
    axis.plot(
        range(1, len(pca.explained_variance_) + 1),  # all PCs
        (pca.explained_variance_ / pca.explained_variance_.sum()) * 100, 'o-')  # % of total variance
    axis.set_xlabel("PC")
    axis.set_ylabel("% variance")
    fig.savefig(os.path.join(analysis.plots_dir, "cll_peaks.all_sites.pca.explained_variance.svg"), bbox_inches='tight')

    # plot PCA
    # with samples colored for all traits
    # sample_colors = samples_to_color(samples, "patient")
    sample_symbols = samples_to_symbol(samples, "unique")

    all_colors = all_sample_colors(samples)
    traits = [
        "patient", "gender", "disease", "IGHV", "ighv_homology"]

    for pc in [0, 1, 2, 3, 4, 5]:
        fig, axis = plt.subplots(5, figsize=(3, 12))
        axis = axis.flatten()
        for i in range(len(traits)):
            for j in range(len(xx)):
                axis[i].scatter(xx.ix[j][pc], xx.ix[j][pc + 1], s=50, color=all_colors[i][j], marker=sample_symbols[j])
            axis[i].set_title(traits[i])
        fig.savefig(os.path.join(analysis.plots_dir, "cll_peaks.all_sites.pca.pc%i_vs_pc%i.svg" % (pc + 1, pc + 2)), bbox_inches="tight")


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


def characterize_regions_structure(df, prefix, output_dir, universe_df=None):
    # use all cll sites as universe
    if universe_df is None:
        universe_df = pd.read_csv(os.path.join("data", "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")

    # compare amplitude and support
    fig, axis = plt.subplots(2, 1, figsize=(5, 10))
    for i, var in enumerate(['fold_change', 'support']):
        sns.distplot(df[[var]], ax=axis[i])
        sns.distplot(universe_df[[var]], ax=axis[i])
        axis[i].set_title(var)
    fig.savefig(os.path.join(output_dir, "%s_regions.amplitude_support.svg" % prefix), bbox_inches="tight")

    # compare 4 measurements of variability
    fig, axis = plt.subplots(2, 2, figsize=(10, 10))
    axis = axis.flatten()
    for i, var in enumerate(['variance', 'std_deviation', 'dispersion', 'qv2']):
        sns.distplot(df[[var]], ax=axis[i])
        sns.distplot(universe_df[[var]], ax=axis[i])
        axis[i].set_title(var)
    fig.savefig(os.path.join(output_dir, "%s_regions.variability.svg" % prefix), bbox_inches="tight")

    # compare genomic regions and chromatin_states
    enrichments = pd.DataFrame()
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

        # g = sns.FacetGrid(col="region", data=data, col_wrap=3, sharey=True)
        # g.map(sns.barplot, "set", "value")
        # plt.savefig(os.path.join(output_dir, "%s_regions.%s.svg" % (prefix, var)), bbox_inches="tight")

        fc = pd.DataFrame(np.log2(both['subset'] / both['all']), columns=['value'])
        fc['variable'] = var

        # append
        enrichments = enrichments.append(fc)

    # save
    enrichments.to_csv(os.path.join(output_dir, "%s_regions.region_enrichment.csv" % prefix), index=True)


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
        for gene in df['gene_name'].dropna():
            handle.write(gene + "\n")
    # write gene names to file
    genes_file = os.path.join(output_dir, "%s_regions.closest_genes.txt" % prefix)
    with open(genes_file, 'w') as handle:
        for gene in df['gene_name'].dropna():
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

    meme_ame(fasta_file, meme_output)


def characterize_regions(analysis, traits, nmin=100):
    """
    Characterize structural-, functionally and in the chromatin regions trait-specific regions.
    """
    all_chrom_counts = Counter(analysis.coverage_qnorm_annotated['chrom'])

    for trait in traits:
        print(trait)
        # Load dataframe with trait-specific regions
        dataframe = pd.read_csv(
            os.path.join(
                analysis.data_dir,
                "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait))

        # Output of all plots/exports will be the same for the same trait
        output_dir = os.path.join(analysis.plots_dir, "%s_regions" % trait)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # ALL REGIONS
        # region's structure
        characterize_regions_structure(df=dataframe, prefix=trait, output_dir=output_dir)

        # plot chromosome distribution of regions
        abs_counts = Counter(dataframe[['chrom', 'start', 'end']].chrom)
        rel_counts = {chrom: count / float(all_chrom_counts[chrom]) for chrom, count in abs_counts.items()}
        abs_counts = np.array(sorted(abs_counts.items(), key=lambda x: x[1], reverse=True))
        rel_counts = np.array(sorted(rel_counts.items(), key=lambda x: x[1], reverse=True))

        fig, axis = plt.subplots(2, figsize=(10, 5))
        sns.barplot(abs_counts[:, 0], abs_counts[:, 1].astype(int), ax=axis[0])
        sns.barplot(rel_counts[:, 0], rel_counts[:, 1].astype(float), ax=axis[1])
        fig.savefig(os.path.join(output_dir, "cll_peaks.%s_regions.chromossomes.svg" % trait), bbox_inches="tight")

        # region's function
        characterize_regions_function(df=dataframe, output_dir=output_dir, prefix=trait)

        # GROUPS OF REGIONS
        # Region characterization of each cluster of regions
        for i, direction in enumerate(np.unique(dataframe['direction'])):
            print(direction)
            # GET REGIONS FROM direction
            df = dataframe[dataframe['direction'] == direction]

            # ignore groups of regions with less than nmin regions
            if df.shape[0] < nmin:
                continue

            # output folder
            outdir = os.path.join(output_dir, "%s_direction%i" % (trait, direction))
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            # region's structure
            characterize_regions_structure(df=df, prefix="%s_direction%i" % (trait, direction), universe_df=analysis.coverage_qnorm_annotated, output_dir=outdir)

            # region's function
            characterize_regions_function(df=df, output_dir=outdir, prefix="%s_direction%i" % (trait, direction))

    # Gather data
    structures = pd.DataFrame()
    lolas = pd.DataFrame()
    motifs = pd.DataFrame()
    gos = pd.DataFrame()
    paths = pd.DataFrame()
    diseases = pd.DataFrame()
    gos2 = pd.DataFrame()
    for trait in traits:
        for direction in ["all", 1, -1]:
            print(trait, direction)
            if direction == "all":
                output_dir = os.path.join("results/plots", "%s_regions" % trait)
            else:
                output_dir = os.path.join("results/plots", "%s_regions" % trait, "%s_direction%i" % (trait, direction))
            # # Gather structure
            results = "%s_regions.region_enrichment.csv" % trait if direction == "all" else "%s_direction%i_regions.region_enrichment.csv" % (trait, direction)
            structure = pd.read_csv(os.path.join(output_dir, results))
            structure["trait"] = trait
            structure["direction"] = direction
            structures = structures.append(structure)

            # Gather LOLA
            lola = pd.read_csv(os.path.join(output_dir, "lola/allEnrichments.txt"), sep="\t")
            lola["trait"] = trait
            lola["direction"] = direction
            lolas = lolas.append(lola)

            # Gather TF motifs
            # parse meme-ame output
            motif = pd.DataFrame(parse_ame(os.path.join(output_dir, "meme")), columns=['motifs', 'q_values'])
            motif['trait'] = trait
            motif['direction'] = direction
            motifs = motifs.append(motif)

            # Gather GO terms
            results = "%s_regions.seq2pathway.csv" % trait if direction == "all" else "%s_direction%i_regions.seq2pathway.csv" % (trait, direction)
            go = pd.read_csv(os.path.join(output_dir, results)).drop_duplicates()
            go["trait"] = trait
            go["direction"] = direction
            gos = gos.append(go)

            # Gather pathways
            results = "%s_regions.pathways.csv" % trait if direction == "all" else "%s_direction%i_regions.pathways.csv" % (trait, direction)
            path = pd.read_csv(os.path.join(output_dir, results), sep="\t")
            path["trait"] = trait
            path["direction"] = direction
            paths = paths.append(path)

            # Gather diseases
            results = "%s_regions.disease.csv" % trait if direction == "all" else "%s_direction%i_regions.disease.csv" % (trait, direction)
            disease = pd.read_csv(os.path.join(output_dir, results), sep="\t")
            disease["trait"] = trait
            disease["direction"] = direction
            diseases = diseases.append(disease)

            # Gather goterms from kobas
            results = "%s_regions.goterms.csv" % trait if direction == "all" else "%s_direction%i_regions.goterms.csv" % (trait, direction)
            go = pd.read_csv(os.path.join(output_dir, results), sep="\t")
            go["trait"] = trait
            go["direction"] = direction
            gos2 = gos2.append(go)

    # Save data
    output_dir = os.path.join("results/plots", "trait_specific")
    structures.to_csv(os.path.join(output_dir, "trait_specific_regions.region_enrichment.csv"), index=False)
    lolas.to_csv(os.path.join(output_dir, "trait_specific_regions.lola_enrichment.csv"), index=False)
    gos.to_csv(os.path.join(output_dir, "trait_specific_regions.goterm_enrichment.csv"), index=False)
    motifs.to_csv(os.path.join(output_dir, "trait_specific_regions.motif_enrichment.csv"), index=False)
    paths.to_csv(os.path.join(output_dir, "trait_specific_regions.pathway_enrichment.csv"), index=False)
    diseases.to_csv(os.path.join(output_dir, "trait_specific_regions.disease_enrichment.csv"), index=False)
    gos2.to_csv(os.path.join(output_dir, "trait_specific_regions.goterm_enrichment-kobas.csv"), index=False)

    # Plot
    # Structure
    structures.columns = ["feature", "value", "type", "trait", "direction"]
    structures["direction"] = structures["direction"].astype(str)
    structures["value"] = structures["value"].astype(float)
    structures["trait"] = structures["trait"] + " - " + structures["direction"]
    structures = pd.melt(structures, ["feature", "trait", "direction", "type"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(structures[structures['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix)
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.region_enrichment.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(structures[structures['direction'] != "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix)
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.region_enrichment.direction.svg"))
    plt.close("all")

    # LOLA
    # create readable name
    cols = ['description', u'cellType', u'tissue', u'antibody', u'treatment', u'dataSource']
    lolas['feature'] = lolas[cols].apply(lambda x: " ".join([i for i in x if i is not pd.np.nan]), axis=1)
    lolas = lolas.rename(columns={"pValueLog": "value"})
    lolas = lolas[["feature", "value", "trait", "direction"]]
    lolas["direction"] = lolas["direction"].astype(str)
    lolas["value"] = lolas["value"].astype(float)
    lolas["trait"] = lolas["trait"] + " - " + lolas["direction"]
    lolas = pd.melt(lolas, ["feature", "trait", "direction", "type"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(lolas[lolas['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.lola_enrichment.all.svg"))
    plt.close("all")
    sns.clustermap(matrix, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.lola_enrichment.all.std_score.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(lolas[lolas['direction'] != "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.lola_enrichment.direction.svg"))
    plt.close("all")
    matrix = pd.pivot_table(lolas[lolas['direction'] != "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.lola_enrichment.direction.std_score.svg"))
    plt.close("all")

    # Motifs
    motifs = motifs.rename(columns={"motifs": "feature", "q_values": "value"})
    motifs["direction"] = motifs["direction"].astype(str)
    motifs["value"] = motifs["value"].astype(float)
    motifs["trait"] = motifs["trait"] + " - " + motifs["direction"]
    motifs = pd.melt(motifs, ["feature", "trait", "direction"])
    # plot only all peaks from trait
    matrix = pd.pivot_table(motifs[motifs['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(-np.log10(matrix))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.motif_enrichment.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(motifs[motifs['direction'] != "all"], values="value", index="feature", columns="trait")
    sns.clustermap(-np.log10(matrix))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.motif_enrichment.direction.svg"))
    plt.close("all")

    # GO terms
    gos = gos.rename(columns={"Name": "feature", "FDR": "value"})
    gos = gos[["feature", "value", "trait", "direction"]]
    gos["direction"] = gos["direction"].astype(str)
    gos["value"] = -np.log10(gos["value"].astype(float))
    gos["trait"] = gos["trait"] + " - " + gos["direction"]
    gos = pd.melt(gos, ["feature", "trait", "direction"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(gos[gos['direction'] == "all"], values="value", index="feature", columns="trait")
    # prune results to get less terms
    matrix_filtered = matrix[matrix.apply(lambda x: all([i > 0.05 for i in x]), axis=1)]
    sns.clustermap(matrix_filtered, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.goterm_enrichment.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(gos[gos['direction'] != "all"], values="value", index="feature", columns="trait")
    # prune results to get less terms
    matrix_filtered = matrix[matrix.apply(lambda x: all([i > 0.05 for i in x]), axis=1)]
    sns.clustermap(matrix_filtered, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.goterm_enrichment.direction.svg"))
    plt.close("all")

    # Pathways
    paths = paths.rename(columns={"#Term": "feature", "P-Value": "value"})
    paths["feature"] = paths["Database"] + " - " + paths["feature"]
    paths = paths[["feature", "value", "trait", "direction"]]
    paths["direction"] = paths["direction"].astype(str)
    paths["value"] = -np.log10(paths["value"].astype(float))
    paths["trait"] = paths["trait"] + " - " + paths["direction"]
    paths = pd.melt(paths, ["feature", "trait", "direction"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(paths[paths['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.pathway_enrichment.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(paths[paths['direction'] != "all"], values="value", index="feature", columns="trait")
    matrix_filtered = matrix[matrix.apply(lambda x: max(x) > 1.5, axis=1)]
    sns.clustermap(matrix_filtered, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.pathway_enrichment.all.filtered.svg"))
    plt.close("all")

    # Diseases
    diseases = diseases.rename(columns={"#Term": "feature", "P-Value": "value"})
    diseases["feature"] = diseases["Database"] + " - " + diseases["feature"]
    diseases = diseases[["feature", "value", "trait", "direction"]]
    diseases["direction"] = diseases["direction"].astype(str)
    diseases["value"] = -np.log10(diseases["value"].astype(float))
    diseases["trait"] = diseases["trait"] + " - " + diseases["direction"]
    diseases = pd.melt(diseases, ["feature", "trait", "direction"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(diseases[diseases['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.disease_enrichment.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(diseases[diseases['direction'] == "all"], values="value", index="feature", columns="trait")
    matrix_filtered = matrix[matrix.apply(lambda x: max(x) > 1.5, axis=1)]
    sns.clustermap(matrix_filtered, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.disease_enrichment.all.fitered.svg"))
    plt.close("all")

    # Go terms 2
    gos2 = gos2.rename(columns={"#Term": "feature", "P-Value": "value"})
    gos2["feature"] = gos2["Database"] + " - " + gos2["feature"]
    gos2 = gos2[["feature", "value", "trait", "direction"]]
    gos2["direction"] = gos2["direction"].astype(str)
    gos2["value"] = -np.log10(gos2["value"].astype(float))
    gos2["trait"] = gos2["trait"] + " - " + gos2["direction"]
    gos2 = pd.melt(gos2, ["feature", "trait", "direction"])
    # only genomic regions
    # plot only all peaks from trait
    matrix = pd.pivot_table(gos2[gos2['direction'] == "all"], values="value", index="feature", columns="trait")
    sns.clustermap(matrix, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.goterms_enrichment-kobas.all.svg"))
    plt.close("all")
    # plot only two directions
    matrix = pd.pivot_table(gos2[gos2['direction'] != "all"], values="value", index="feature", columns="trait")
    matrix_filtered = matrix[matrix.apply(lambda x: max(x) > 1.5, axis=1)]
    sns.clustermap(matrix_filtered, standard_scale=1, figsize=(20, 20))
    plt.savefig(os.path.join(output_dir, "trait_specific_regions.goterms_enrichment-kobas.all.fitered.svg"))
    plt.close("all")


def export_matrices(analysis, sel_samples, labels, trait, validation=""):
    """
    Save matrices used for classification.
    """
    print("Trait:%s" % trait)
    print("%i samples with trait annotated" % len(sel_samples))
    print(Counter(labels))

    # normalize
    X = pd.DataFrame(normalize(analysis.coverage_qnorm_annotated[[s.name for s in sel_samples]].T), index=[s.name for s in sel_samples])
    # add labels
    X['target'] = labels

    # write out with sample name
    outfile_file = os.path.join(analysis.data_dir, "ml_matrix.%s.%scsv" % (trait, validation))
    X.T.to_csv(outfile_file, index=True)


def trait_analysis(analysis, samples, traits):
    """
    Run trait classification (with independent validation if possible for that trait)
    on all samples with known annotation of the trait.
    """
    for trait in traits:
        # cross-validated classification
        sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]
        labels = np.array([getattr(s, trait) for s in sel_samples])
        classify_samples(analysis, sel_samples, labels, trait, rerun=True)
        export_matrices(analysis, sel_samples, labels, trait)

        # random classifiers
        sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]
        labels = np.array([getattr(s, trait) for s in sel_samples])
        classification_random(analysis, sel_samples, labels, trait, n=100)


def join_trait_specific_regions(analysis, traits):
    # Put trait-specific chromatin regions in one matrix
    features = pd.DataFrame()
    for trait in traits:
        file_name = os.path.join(
            analysis.data_dir,
            "cll_peaks.%s_significant.classification.random_forest.loocv.dataframe.csv" % trait)
        try:
            df = pd.read_csv(file_name)
            # assert df.shape == df.dropna().shape  # check there are no nans - in fact there can be (e.g. gene_name column)
            df['trait'] = trait

            df = pd.merge(
                df[["chrom", "start", "end", "change", "direction", "importance", "trait"]],
                analysis.coverage_qnorm_annotated,
                on=["chrom", "start", "end"])

            features = features.append(df, ignore_index=True)
        except IOError:
            print("Trait %s did not generate any associated regions" % trait)

    # # Save whole dataframe as csv
    features.to_csv(os.path.join(analysis.data_dir, "cll.trait-specific_regions.csv"), index=False)


def characterize_regions_chromatin(analysis, traits, extend=False):
    """
    For each trait-associated region, get ratios of active/repressed and poised/repressed chromatin,
    across all patients or groups of patients with same IGHV mutation status.

    Visualize histone mark ratios dependent on the positive- or negative-association of regions with
    accessibility signal.
    """
    def coverage_around(sites, samples, diameter=1000):
        sites2 = sites.slop(b=diameter, genome="hg19")
        sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites2]
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

        # Add interval description to df
        # get original strings of intervals (before slop)
        original_sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites]
        # get indexes of intervals which we got coverage from
        original_sites_str_covered = [original_sites_str[i] for i, s in enumerate(coverage.index) if s in sites_str]
        assert len(original_sites_str_covered) == len(coverage.index)
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            original_sites_str_covered
        )
        coverage["chrom"] = [x[0] for x in ints]
        coverage["start"] = [int(x[1]) for x in ints]
        coverage["end"] = [int(x[2]) for x in ints]

        # save to disk
        coverage.to_csv(os.path.join("data", "cll_peaks.raw_coverage.histmods.%i_extended.tsv" % diameter), sep="\t", index=True)
        # coverage = pd.read_csv(os.path.join("data", "cll_peaks.raw_coverage.histmods.%i_extended.tsv" % diameter), sep="\t", index_col=0)

        # normalize and log2
        to_norm = coverage.iloc[:, :len(samples)]
        coverage_qnorm = pd.DataFrame(
            normalize_quantiles_r(np.array(to_norm)),
            index=to_norm.index,
            columns=to_norm.columns
        )
        # Log2 transform
        coverage_qnorm = np.log2(1 + coverage_qnorm)

        coverage_qnorm = coverage_qnorm.join(coverage[['chrom', 'start', 'end']])
        coverage_qnorm.to_csv(os.path.join("data", "cll_peaks.raw_coverage.histmods.%i_extended.qnorm.log2.tsv"), sep="\t", index=False)

        return coverage_qnorm

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

    # samples
    samples = analysis.samples
    # read in dataframe with counts over CLL ATAC-seq peaks
    features = pd.read_csv(os.path.join(analysis.data_dir, "cll.trait-specific_regions.csv"))
    if extend:
        # alternatively, replace ChIPmentation counts with new read count around peaks (e.g. 1kb)
        chip_samples = [s for s in analysis.samples if s.library == "ChIPmentation"]
        sites = pybedtools.BedTool(os.path.join(analysis.data_dir, "cll_peaks.bed"))
        features_extended = coverage_around(sites, chip_samples, diameter=1000)
        for sample in chip_samples:
            features[sample.name] = features_extended[sample.name].ix[features.index].reset_index(drop=True)

        # or
        # to_norm = features[features.columns[features.columns.str.contains("CLL")]]
        # to_append = features[features.columns[~features.columns.str.contains("CLL")]]
        # features.columns[features.columns.str.contains("CLL")]
        # n = pd.DataFrame(
        #     normalize_quantiles_r(np.array(to_norm)),
        #     index=to_norm.index,
        #     columns=to_norm.columns
        # )
        # features_extended = pd.concat([to_append, n], axis=1)

    # get intensity of chromatin marks
    groups = ["all", "uCLL", "iCLL", "mCLL"]
    marks = ["H3K27ac", "H3K4me1", "H3K27me3"]
    # across all patients
    for group in groups:
        for mark in ["H3K27ac", "H3K4me1", "H3K27me3"]:
            features["%s_intensity_%s" % (group, mark)] = get_intensities(features, group, mark)

    # compute ratios of chromatin marks
    # across all patients and for each each IGHV class
    for group in marks:
        for ratio in ["A:R", "P:R"]:
            features["%s_ratio_%s" % (group, ratio)] = calculate_ratio(features, group, ratio)

    # save dataframe with intensities/ratios of chromatin marks per peak
    features.to_csv(os.path.join(
        analysis.data_dir,
        "cll.trait-specific_regions.%shistone_intensities_ratios.csv" % ("extended." if extend else "")), index=False)

    # extend = True
    # traits = ['IGHV']
    # features = pd.read_csv(os.path.join(
    #     analysis.data_dir,
    #     "trait_specific",
    #     "cll.trait-specific_regions.%shistone_intensities_ratios.csv" % "extended." if extend else ""))

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
        sample_colors = [trait_colors[s.ighv_mutation_status] if not pd.isnull(s.ighv_mutation_status) else "gray" for s in samples]

        # colormap for heatmap values
        cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)

        # heatmap
        dend = sns.clustermap(
            df[[s.name for s in samples]],
            cmap=cmap,
            standard_scale=0,
            col_colors=sample_colors,
            yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s.%sclustermap.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # order rows by ATAC-seq dendrogram
        chrom = df[df.columns[df.columns.str.contains("intensity")]]
        order = dend.dendrogram_row.dendrogram['leaves']
        chrom = chrom.ix[order]

        # order columns by histone mark
        chrom = chrom[["_".join([g, "intensity", mark]) for mark in marks for g in groups]]
        sns.heatmap(chrom, yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_mark.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # order columns by group
        chrom = chrom[["_".join([g, "intensity", mark]) for g in groups for mark in marks]]
        sns.heatmap(chrom, yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_group.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # Z-scored, order columns by histone mark
        chrom = chrom[["_".join([g, "intensity", mark]) for mark in marks for g in groups]]
        sns.heatmap(chrom.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_mark.zscore_rows.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # Z-scored, order columns by group
        chrom = chrom[["_".join([g, "intensity", mark]) for g in groups for mark in marks]]
        sns.heatmap(chrom.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_group.zscore_rows.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        chrom_specfic = chrom[chrom.columns[~chrom.columns.str.contains("all")]]

        # Z-scored, order columns by histone mark
        chrom_specfic = chrom_specfic[["_".join([g, "intensity", mark]) for mark in marks for g in groups]]
        sns.heatmap(chrom_specfic.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_mark.specific.zscore_rows.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # Z-scored, order columns by group
        chrom_specfic = chrom_specfic[["_".join([g, "intensity", mark]) for g in groups for mark in marks]]
        sns.heatmap(chrom_specfic.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
        plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_group.specific.zscore_rows.svg" % (trait, "extended." if extend else "")), bbox_inches="tight")
        plt.close("all")

        # Z score separately each mark
        for mark in ["H3K27ac", "H3K4me1", "H3K27me3"]:
            df = chrom_specfic[chrom_specfic.columns[chrom_specfic.columns.str.contains(mark)]]

            # Z-scored, order columns by histone mark
            df = df[["_".join([g, "intensity", mark]) for g in ["uCLL", "iCLL", "mCLL"]]]
            sns.heatmap(df.apply(lambda x: (x - x.mean()) / x.std(), axis=0), yticklabels=False)
            plt.savefig(os.path.join(analysis.plots_dir, "%s_histones.%sordered_mark.specific.zscore_rows.%s_only.separate.svg" % (trait, "extended." if extend else "", mark)), bbox_inches="tight")
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
    for trait in ["IGHV"]:
        p = features_slim[features_slim['trait'] == trait].dropna()

        g = sns.FacetGrid(p, col="group", row="direction", legend_out=True, row_order=[-1, 1], margin_titles=True)
        g.map(sns.violinplot, "mark", "intensity", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_intensity.%smark_centric.svg" % (trait, "extended." if extend else "")),
            bbox_inches='tight')

        g = sns.FacetGrid(p, col="group", row="mark", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "intensity", split=True, order=[-1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_intensity.%sdirection_centric.svg" % (trait, "extended." if extend else "")),
            bbox_inches='tight')

        g = sns.FacetGrid(p, row="mark", col="direction", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "group", "intensity", split=True, order=["uCLL", "iCLL", "mCLL"])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_intensity.%sgroup_centric.svg" % (trait, "extended." if extend else "")),
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
    for trait in ["IGHV"]:
        p = features_slim[features_slim['trait'] == trait].dropna()

        g = sns.FacetGrid(p, col="group", row="direction", legend_out=True)
        g.map(sns.violinplot, "type", "ratio", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_ratios.%smark_centric.svg" % (trait, "extended." if extend else "")),
            bbox_inches='tight')

        g = sns.FacetGrid(p, col="group", row="type", legend_out=True)
        g.map(sns.violinplot, "direction", "ratio", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_ratios.%sdirection_centric.svg" % (trait, "extended." if extend else "")),
            bbox_inches='tight')

        g = sns.FacetGrid(p, row="type", col="direction", legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "group", "ratio", split=True, order=["uCLL", "iCLL", "mCLL"])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.chromatin_ratios.%sgroup_centric.svg" % (trait, "extended." if extend else "")),
            bbox_inches='tight')


def characterize_regions_expression(analysis, traits, extend=False):
    """
    For each trait-associated region, get of genes in those regions,
    across all patients or groups of patients with same IGHV mutation status.

    Visualize dependent on the positive- or negative-association of regions with
    accessibility signal.
    """
    import itertools

    def get_mean_expression(df, group):
        """
        Get mean expression of genes across a subset of samples (uCLL or mCLL).
        """
        samples = [s for s in analysis.samples if s.library == "RNA-seq"]
        # subset samples according to requested group (all, uCLL, iCLL or mCLL)
        if group == "all":
            g = samples
        else:
            g = [s for s in samples if s.ighv_group == group]
        # mean of that mark across the subset of samples
        return pd.DataFrame(df[[s.name for s in g]].mean(axis=1), columns=["value"]).reset_index()

    def calculate_fold_change(df, group1, group2):
        """
        Get mean fold-change of gene expression in two groups of samples.
        """
        samples = [s for s in analysis.samples if s.library == "RNA-seq"]
        # subset samples according to requested group (all, uCLL, iCLL or mCLL)
        if group1 == "all":
            u = samples
        else:
            u = [s for s in samples if s.ighv_group == group1]
        if group2 == "all":
            m = samples
        else:
            m = [s for s in samples if s.ighv_group == group2]

        # get ratio of the requested histone marks
        u_s = df[[s.name for s in u]]
        m_s = df[[s.name for s in m]]
        return pd.DataFrame(np.log2(((1 + u_s.mean(1)) / (1 + m_s.mean(1)))), columns=["value"]).reset_index()

    # read in trait-specific regions
    features = pd.read_csv(os.path.join(analysis.data_dir, "cll.trait-specific_regions.csv"))
    # Get ensembl_gene ids mapping
    features = pd.merge(
        features,
        analysis.gene_annotation, on=["chrom", "start", "end"])

    # Get gene expression matrix
    expression_genes = pd.read_csv(os.path.join(analysis.data_dir, "cll_expression_matrix.log2.csv")).set_index("ensembl_gene_id")

    # Quantile normalization
    to_norm = expression_genes[[s.name for s in analysis.samples if s.library == "RNA-seq"]]
    expression_genes_qnorm = pd.DataFrame(
        normalize_quantiles_r(np.array(to_norm)),
        index=to_norm.index,
        columns=to_norm.columns
    )

    # Get expression intensities and ratios on each group of trait-specific regions
    for trait in traits:
        # subset trait
        feature_specific = features[features["trait"] == trait]

        # get intensity of gene expression
        exp_int = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        # across all patients
        for group in groups:
            for direction in ["all", 1, -1]:
                # Get genes within regions of trait, direction-specific
                if direction == "all":
                    index = feature_specific["ensembl_gene_id"].dropna().drop_duplicates()
                else:
                    index = feature_specific[feature_specific["direction"] == direction]["ensembl_gene_id"].dropna().drop_duplicates()
                expression_genes_regions = expression_genes_qnorm.ix[index]
                df = get_mean_expression(expression_genes_regions, group)
                df["kind"] = "intensity"
                df["group"] = group
                df["direction"] = direction
                exp_int = exp_int.append(df, ignore_index=True)

        # compute ratios of gene expression
        exp_ratio = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        for group1, group2 in itertools.permutations(groups, 2):
            for direction in ["all", 1, -1]:
                # Get genes within regions of trait, direction-specific
                if direction == "all":
                    index = feature_specific["ensembl_gene_id"].dropna().drop_duplicates()
                else:
                    index = feature_specific[feature_specific["direction"] == direction]["ensembl_gene_id"].dropna().drop_duplicates()
                expression_genes_regions = expression_genes_qnorm.ix[index]
                df = calculate_fold_change(expression_genes_regions, group1, group2)
                df["kind"] = "ratio"
                df["group"] = "{}_{}".format(group1, group2)
                df["direction"] = direction
                exp_ratio = exp_ratio.append(df, ignore_index=True)

        # save dataframe with intensities/ratios of chromatin marks per peak
        pd.concat([exp_int, exp_ratio]).to_csv(os.path.join(
            analysis.data_dir,
            "cll.trait-specific_regions.%s.expression_intensities_ratios.csv" % trait), index=False)

        g = sns.FacetGrid(exp_int, col="direction", col_wrap=3, legend_out=True, col_order=["all", -1, 1], margin_titles=True)
        g.map(sns.violinplot, "group", "value", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_intensity.group_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_int, col="group", col_wrap=3, legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "value", split=True, order=["all", -1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_intensity.direction_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_ratio, col="direction", col_wrap=3, legend_out=True, col_order=["all", -1, 1], margin_titles=True)
        g.map(sns.violinplot, "group", "value", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.group_centric.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_ratio, col="group", col_wrap=3, legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "value", split=True, order=["all", -1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.svg" % trait),
            bbox_inches='tight')

        # Filter out non-expressed genes across all samples
        d = exp_int[(exp_int['group'] == 'all') & (exp_int['direction'] == 'all')]
        expressed_genes = d[d["value"] >= 1]["ensembl_gene_id"]

        # new dataframes
        exp_int2 = exp_int[exp_int["ensembl_gene_id"].isin(expressed_genes)]
        exp_ratio2 = exp_ratio[exp_ratio["ensembl_gene_id"].isin(expressed_genes)]

        g = sns.FacetGrid(exp_int2, col="direction", col_wrap=3, legend_out=True, col_order=["all", -1, 1], margin_titles=True)
        g.map(sns.violinplot, "group", "value", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_intensity.group_centric.expressed.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_int2, col="group", col_wrap=3, legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "value", split=True, order=["all", -1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_intensity.direction_centric.expressed.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_ratio2, col="direction", col_wrap=3, legend_out=True, col_order=["all", -1, 1], margin_titles=True)
        g.map(sns.violinplot, "group", "value", split=True)
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.group_centric.expressed.svg" % trait),
            bbox_inches='tight')

        g = sns.FacetGrid(exp_ratio2, col="group", col_wrap=3, legend_out=True, margin_titles=True)
        g.map(sns.violinplot, "direction", "value", split=True, order=["all", -1, 1])
        sns.despine()
        plt.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.expressed.svg" % trait),
            bbox_inches='tight')

        # Get fold-change over the other regions within each group
        exp_int3 = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        for group in groups:
            df = pd.Series()
            a = exp_int[
                (exp_int["group"] == group) &
                (exp_int["direction"] == 1)]["value"].mean()
            b = exp_int[
                (exp_int["group"] == group) &
                (exp_int["direction"] == -1)]["value"].mean()
            df["value"] = np.log2(a / b)
            df["direction"] = 1
            exp_int3[group] = df

        fig, axis = plt.subplots(1)
        sns.barplot(groups, exp_int3.T["value"], ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.mean_fold_change_in_regions.svg" % trait),
            bbox_inches='tight')

        # Get fold-change over the other regions within each group and over all samples
        exp_int4 = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        for group in groups:
            df = pd.Series()
            df["value"] = np.log2(exp_int3[group]["value"] / exp_int3["all"]["value"])
            exp_int4[group] = df

        fig, axis = plt.subplots(1)
        sns.barplot(groups, exp_int4.T["value"], ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.mean_fold_change_in_regions.over_all_samples.svg" % trait),
            bbox_inches='tight')

        # Filter out non-expressed genes across all samples
        d = exp_int[(exp_int['group'] == 'all') & (exp_int['direction'] == 'all')]
        expressed_genes = d[d["value"] >= 1]["ensembl_gene_id"]
        exp_int2 = exp_int[exp_int["ensembl_gene_id"].isin(expressed_genes)]

        # Get fold-change over the other regions within each group
        exp_int3 = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        for group in groups:
            df = pd.Series()
            a = exp_int2[
                (exp_int2["group"] == group) &
                (exp_int2["direction"] == 1)]["value"].mean()
            b = exp_int2[
                (exp_int2["group"] == group) &
                (exp_int2["direction"] == -1)]["value"].mean()
            df["value"] = np.log2(a / b)
            df["kind"] = "intensity"
            df["group"] = group
            exp_int3[group] = df

        fig, axis = plt.subplots(1)
        sns.barplot(groups, exp_int3.T["value"], ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.expressed.mean_fold_change_in_regions.svg" % trait),
            bbox_inches='tight')

        # Get fold-change over the other regions within each group and over all samples
        exp_int4 = pd.DataFrame()
        groups = ["all", "uCLL", "iCLL", "mCLL"]
        for group in groups:
            df = pd.Series()
            df["value"] = np.log2(exp_int3[group]["value"] / exp_int3["all"]["value"])
            exp_int4[group] = df

        fig, axis = plt.subplots(1)
        sns.barplot(groups, exp_int4.T["value"], ax=axis)
        fig.savefig(os.path.join(
            analysis.plots_dir, "cll.%s-specific_regions.expression_ratio.direction_centric.expressed.mean_fold_change_in_regions.over_all_samples.svg" % trait),
            bbox_inches='tight')


def main():
    # Parse arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = add_args(parser)
    args = parser.parse_args()

    # Create Project and Analysis objects or load them
    if args.generate:
        # Start project
        prj = Project("metadata/project_config.yaml")
        prj.add_sample_sheet()

        # annotated samples with a few more things:
        prj.samples = annotate_samples(prj.samples, prj.sheet.df.columns.tolist())

        # Start analysis object
        analysis = Analysis(
            data_dir=os.path.join(".", "data_submission"),
            plots_dir=os.path.join(".", "results", "plots"),
            samples=prj.samples,
            pickle_file=os.path.join(".", "data", "analysis.pickle")
        )
        # pair analysis and Project
        analysis.prj = prj
    else:
        analysis = pickle.load(open(os.path.join(".", "data", "analysis.pickle"), 'rb'))

    # Create subsets of samples
    atac_seq_samples = [sample for sample in analysis.samples if sample.library == "ATAC-seq"]
    chipmentation_samples = [sample for sample in analysis.samples if sample.library == "ChIPmentation"]
    chromatin_samples = [sample for sample in analysis.samples if sample.library in ["ATAC-seq", "ChIPmentation"]]
    rna_samples = [sample for sample in analysis.samples if sample.library == "RNA-seq"]

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
        analysis.measure_coverage(chromatin_samples)
        # normalize coverage values
        analysis.normalize_coverage_quantiles(atac_seq_samples, chipmentation_samples)  # normalize them seperately
        # Annotate peaks with closest gene, chromatin state,
        # genomic location, mean and variance measurements across samples
        analysis.annotate(atac_seq_samples)  # samples on which peak metrics will be calculated on
    else:
        analysis.coverage = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.raw_coverage.tsv"), sep="\t", index_col=0)
        analysis.coverage_qnorm = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.coverage_qnorm.log2.tsv"), sep="\t")
        analysis.coverage_qnorm_annotated = pd.read_csv(os.path.join(analysis.data_dir, "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")

    # Characterize all CLL regions as a whole
    # region's structure
    characterize_regions_structure(df=analysis.coverage_qnorm_annotated, prefix="all_regions", output_dir=analysis.plots_dir)
    # region's function
    output_dir = os.path.join(analysis.data_dir, "%s_peaks" % "all_regions")
    characterize_regions_function(df=analysis.coverage_qnorm_annotated, output_dir=output_dir, prefix="all_regions")

    # Plot oppenness across peaks/samples
    analysis.plot_coverage()
    analysis.plot_variance()
    # Observe exponential fit to the coeficient of variation
    analysis.plot_qv2_fit()
    # Correlate with expression
    analysis.correlate_expression()
    analysis.correlate_methylation(atac_seq_samples)

    # INTER-PATIENT VARIATION (Figure 2)
    # cross-cohort variation at gene level
    analysis.gene_oppeness_across_samples(atac_seq_samples)
    # try seeing what most variable regions are with LOLA (no good results - not included)
    analysis.variability()
    # inter-group variation
    analysis.inspect_variability(atac_seq_samples)
    # unsupervised analysis
    unsupervised(analysis, atac_seq_samples)

    # TRAIT-SPECIFIC ANALYSIS (Figure 3)
    traits = ["IGHV"]
    trait_analysis(analysis, atac_seq_samples, traits)

    join_trait_specific_regions(analysis, traits)

    # characterize trait-specific regions
    # structurally, functionaly and in light of histone marks
    characterize_regions(analysis, traits)
    characterize_regions_chromatin(analysis, traits, extend=False)

    # # save "digested" clinical sheet to disk
    if args.generate:
        prj.sheet.asDataFrame().drop_duplicates().to_csv("clinical_annotation_digested.csv", index=False)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
