#!/usr/bin/env python

"""

"""

import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from pipelines.models import Project
from networkx.readwrite import json_graph
import json


sns.set(style="whitegrid")


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


def create_graph(network_file):
    # read in network file
    df = pd.read_csv(network_file, sep="\t")

    # create network
    G = nx.DiGraph()
    for i in df.index:
        G.add_edge(df.ix[i]['TF'], df.ix[i]['gene'], weight=df.ix[i]['interaction_score'])
    return G


def intersect_sum_weights(G1, G2):
    """
    Intersect two graphs, adding weights
    and keeping track of how many graphs
    have contributed to that edge's weight.
    """
    # the new graph
    G3 = nx.DiGraph()

    for u, v in G1.edges_iter():
        weight1 = G1[u][v]['weight']
        if (u, v) in G2.edges():
            weight2 = G2[u][v]['weight']
            # sum weights
            weight = weight2 + weight1

            # keep track of how many graphs have contributed to this edge's weight
            count = 0
            if "count" in G1[u][v].keys():
                count = G1[u][v]['count'] + count
            else:
                count = 1
            if "count" in G2[u][v].keys():
                count = G2[u][v]['count'] + count
            else:
                count = 2

            # add that value

            # add edge to new graph
            G3.add_edge(u, v, weight=weight, attr_dict={'count': count})
    return G3


def average_weights(G):
    """
    Averages the weight of each edge
    by the number of graphs which have
    contributed to the edge's weight
    (see intersect_sum_weights).
    """
    for u, v in G.edges_iter():
        weight = G[u][v]['weight']
        if 'count' in G[u][v].keys():
            weight /= G[u][v]['count']
        G.add_edge(u, v, weight=weight)
    return G


def describe_graph(G):
    """Graph description"""

    # GRAPH DESCRIPTION
    graph_desc = pd.Series()
    # n. nodes
    graph_desc["number_of_nodes"] = G.number_of_nodes()
    # n. edges
    graph_desc["number_of_edges"] = G.number_of_edges()
    # n. of selfloops
    graph_desc["number_of_selfloops"] = len(G.selfloop_edges())

    # density
    graph_desc["average_shortest_path_length"] = nx.average_shortest_path_length(G)
    # connectivity
    # graph_desc.append(pd.Series(nx.degree_assortativity_coefficient(G), name="degree_assortativity_coefficient"))
    graph_desc["degree_pearson_correlation_coefficient"] = nx.degree_pearson_correlation_coefficient(G)

    # NODE DESCRIPTION
    node_desc = list()
    # n. of neighbours
    node_desc.append(pd.Series(G.degree(), name="degree"))
    node_desc.append(pd.Series(nx.average_neighbor_degree(G), name="average_neighbor_degree"))
    # n. of outgoing
    outgoing = pd.Series(G.in_degree(), name="in_degree")
    node_desc.append(outgoing)
    # n. of incoming
    incoming = pd.Series(G.out_degree(), name="out_degree")
    node_desc.append(incoming)
    # fold change out/in
    ratio = np.log2(outgoing + 1) - np.log2(incoming + 1)
    node_desc.append(pd.Series(ratio, name="out_in_degree_fold_change"))

    # centrality
    # degree based
    node_desc.append(pd.Series(nx.degree_centrality(G), name="degree_centrality"))
    node_desc.append(pd.Series(nx.in_degree_centrality(G), name="in_degree_centrality"))
    node_desc.append(pd.Series(nx.out_degree_centrality(G), name="out_degree_centrality"))
    # closest-path based
    # node_desc.append(pd.Series(nx.closeness_centrality(G), name="closeness_centrality"))
    # node_desc.append(pd.Series(nx.betweenness_centrality(G), name="betweenness_centrality"))
    # # eigenvector-based
    # node_desc.append(pd.Series(nx.eigenvector_centrality(G), name="eigenvector_centrality"))
    # node_desc.append(pd.Series(nx.katz_centrality_numpy(G), name="katz_centrality"))
    # # load-based
    # node_desc.append(pd.Series(nx.load_centrality(G), name="load_centrality"))

    return (graph_desc, pd.DataFrame(node_desc).T)


def color_nodes_with_degree_ratio(G):
    # Node in/out degree ratio
    degree_ratio = pd.Series(G.in_degree()) - pd.Series(G.out_degree())

    colors = dict()
    for u in degree_ratio.index:
        if degree_ratio.ix[u] > 0:
            colors[u] = "#FF6666"
        elif degree_ratio.ix[u] < 0:
            colors[u] = "#CCCCFF"
        else:
            colors[u] = "#808080"

    nx.set_node_attributes(G, 'color', colors)
    return G


def graph_to_json_d3(G, json_file):
    # write network to disk
    # this can be visualized with D3.js (e.g. http://bl.ocks.org/mbostock/4062045#index.html)
    json_data = json_graph.node_link_data(G)
    with open(json_file, "w") as handle:
        json.dump(json_data, handle)


def json_to_graph(json_file):
    return json_graph.node_link_graph(json.load(open(json_file, "r")))


def samples_to_networks(samples):
    for i, sample in enumerate(samples):
        print(sample)
        # Create sample graph
        graph_file = os.path.join(data_dir, "footprints", sample.name + ".piq.TF-TF_interactions.tsv")
        G = create_graph(graph_file)

        # Describe network
        graph_desc, node_desc = describe_graph(G)
        graph_desc.name = sample.name
        node_desc["sample"] = sample.name
        if i == 0:
            graph = [graph_desc]
            nodes = node_desc
        else:
            graph.append(graph_desc)
            nodes = nodes.append(node_desc)

        # Intersect graphs, keeping weight from each
        if i == 0:
            master_graph = G
        else:
            master_graph = intersect_sum_weights(master_graph, G)
    graph = pd.DataFrame(graph)

    G = average_weights(master_graph)

    return G, graph, nodes


def plot_graph_attributes(graph, nodes, prefix):
    # INDIVIDUAL NETWORKS
    # Plot graph attributes
    graph['sample'] = graph.index
    # in stripplots
    g = sns.PairGrid(
        graph.sort("number_of_nodes", ascending=False),
        x_vars=graph.columns[:4], y_vars=["sample"],
        size=10, aspect=.25)
    g.map(sns.stripplot, size=10, orient="h",
          palette="Reds_r", edgecolor="gray")
    plt.savefig(prefix + ".svg", bbox_inches="tight")

    # Plot node attributes
    nodes['TF'] = nodes.index
    df = pd.melt(nodes, id_vars=['TF', 'sample'])
    # all values of all samples
    g = sns.FacetGrid(df, col="variable", sharey=False, col_wrap=3, margin_titles=True, size=4)
    g.map(sns.barplot, "TF", "value")
    plt.savefig(prefix + ".nodes.svg", bbox_inches="tight")

    # mean vs -std across samples, in a grid of variables
    # calculate mean across samples for each variable
    mean = df.groupby(["TF", "variable"]).aggregate(np.mean).reset_index()
    mean.columns = ["TF", "variable", "mean"]
    # calculate std across samples for each variable
    std = df.groupby(["TF", "variable"]).aggregate(np.std).reset_index()
    std.columns = ["TF", "variable", "std"]
    # merge both (mean, std)
    df2 = pd.merge(mean, std).dropna()
    # plot
    g = sns.FacetGrid(df2, col="variable", hue="TF", sharex=False, sharey=False, col_wrap=4, margin_titles=True)
    g.map(plt.scatter, "mean", "std")
    plt.savefig(prefix + ".nodes.stats.svg", bbox_inches="tight")


def subtract_networks_nodes(G1, G2):
    """
    Subtract two graphs, by removing weights from common edges.
    """
    # the new graph
    G3 = nx.DiGraph()

    for u in G1.nodes_iter():
        if u not in G2.nodes():
            G3.add_node(u)
    return G3


def subtract_networks_edges(G1, G2):
    """
    Subtract two graphs, by removing common edges.
    """
    # the new graph
    G3 = nx.DiGraph()

    for u, v in G1.edges_iter():
        weight1 = G1[u][v]['weight']
        if (u, v) not in G2.edges():
            G3.add_edge(u, v, weight=weight1)
    return G3


def subtract_networks_edges_weights(G1, G2):
    """
    Subtract two graphs, by removing weights from common edges.
    """
    # the new graph
    G3 = nx.DiGraph()

    for u, v in G1.edges_iter():
        weight1 = G1[u][v]['weight']
        if (u, v) in G2.edges():
            weight2 = G2[u][v]['weight']
            # sum weights
            weight = weight1 - weight2

            # add edge to new graph
            G3.add_edge(u, v, weight=weight)
        else:
            G3.add_edge(u, v, weight=weight1)
    return G3


def filter_networks_edges(G, threshold):
    """
    Subtract two graphs, by removing common edges.
    """
    # the new graph
    G2 = nx.DiGraph()

    for u, v in G.edges_iter():
        weight1 = G[u][v]['weight']
        if weight1 >= threshold:
            G2.add_edge(u, v, weight=weight1)
    return G2


# Get path configuration
data_dir = os.path.join('.', "data")
results_dir = os.path.join('.', "results")
plots_dir = os.path.join(results_dir, "plots")


# simpler ighv network comparison
new_graph_file = os.path.join("netx/merged-samples_all_all.piq.TF-gene_interactions.filtered.tsv")
G = create_graph(new_graph_file)

tfs = pd.read_table(new_graph_file)['TF'].tolist()

# degree of all nodes
d = pd.Series(G.degree())
d.sort()
fig, axis = plt.subplots(1)
colors = ["#d7191c" if i in tfs else "#669fb8" for i in d.index]
axis.scatter(d.rank(ascending=False, method='first'), np.log2(d), linewidth=0, color=colors)
fig.savefig("all.degree.svg", bbox_inches="tight")

# subnetwork of most connected nodes
# get subgraph
g = G.subgraph(d[d > 200].index)
nx.set_node_attributes(g, 'abs_degree', {n: int(d.to_dict()[n]) for n in g.nodes()})
nx.write_gexf(g, "merged-samples_all_all.piq.TF-gene_interactions.filtered.top.gexf")


# IGHV comparison
graph_file = os.path.join(data_dir, "merged-samples_mutated_False.piq.TF-gene_interactions.filtered.tsv")
uG = create_graph(graph_file)

graph_file = os.path.join(data_dir, "merged-samples_mutated_True.piq.TF-gene_interactions.filtered.tsv")
mG = create_graph(graph_file)

g = uG.subgraph(d[d > 200].index)
nx.set_node_attributes(g, 'abs_degree', {n: int(d.to_dict()[n]) for n in g.nodes()})
nx.write_gexf(g, "merged-samples_mutated_False.piq.TF-gene_interactions.filtered.top.gexf")

g = mG.subgraph(d[d > 200].index)
nx.set_node_attributes(g, 'abs_degree', {n: int(d.to_dict()[n]) for n in g.nodes()})
nx.write_gexf(g, "merged-samples_mutated_True.piq.TF-gene_interactions.filtered.top.gexf")


unode_degree = pd.Series(uG.degree(), name="degree")
mnode_degree = pd.Series(mG.degree(), name="degree")

# Compare nodes
df = pd.DataFrame([unode_degree, mnode_degree]).T
df.columns = ['u', 'm']
# normalize
df = df.apply(lambda x: (x - x.mean() / x.std()), axis=0)
#
df = np.log2(1 + df.fillna(0))
plt.scatter(df['m'], df['u'])
plt.plot([0, 4], [0, 4], '--')


# let's look at the difference now
diff = df['u'] - df['m']

# coverage_qnorm_annotated = pd.read_csv(os.path.join("data", "cll_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t")
# names = coverage_qnorm_annotated.columns[coverage_qnorm_annotated.columns.str.contains("CLL")]
# # average openness across samples, sum oppenness across different elements of same gene
# weights = np.log10(coverage_qnorm_annotated[names + ["gene_name"]].groupby(['gene_name']).apply(np.mean, axis=1).groupby(level=[0]).sum())

# weights = (weights - weights.min()) / (weights.max() - weights.min())

# score = diff * weights
# score = score.dropna()

fig, axis = plt.subplots(1)
sns.distplot(diff, ax=axis)
fig.savefig("all_fc.distplot.svg", bbox_inches="tight")

fig, axis = plt.subplots(1)
colors = list()
for i in diff:
    if i > 1:
        colors.append("#d55e00")
    elif i < -1:
        colors.append("#0072b2")
    else:
        colors.append("gray")

axis.scatter(diff.rank(ascending=False, method='first'), diff, linewidth=0, color=colors)
fig.savefig("all_fc.rank_fc.svg", bbox_inches="tight")

# cll genes' change
diffcll = diff.ix[cllgenes]
diffcll.sort(ascending=False)
fig, axis = plt.subplots(1)
colors = list()
for i in diffcll:
    if i > 1:
        colors.append("#d55e00")
    elif i < -1:
        colors.append("#0072b2")
    else:
        colors.append("gray")
axis.scatter(diffcll.rank(ascending=False), diffcll, linewidth=0, color=colors)
fig.savefig("all_fc.rank_fc.cllgenes.svg", bbox_inches="tight")


# Visualize single nodes in the two networks
examples = ["PAX9", "CD9", "CD38", "CD22"]

tfs = G.nodes()

for gene in examples:
    p = nx.spring_layout(G.subgraph([n for n in [i for i in nx.all_neighbors(uG, gene)] if n in tfs] + [gene]))
    nx.draw(
        uG.subgraph([n for n in [i for i in nx.all_neighbors(uG, gene)] if n in tfs] + [gene]),
        pos=p
    )
    nx.draw_networkx_labels(
        uG.subgraph([n for n in [i for i in nx.all_neighbors(uG, gene)] if n in tfs] + [gene]),
        pos=p
    )
    nx.write_gexf(uG.subgraph([n for n in [i for i in nx.all_neighbors(uG, gene)] if n in tfs] + [gene]), "graph_modules.%s.uCLL.gexf" % gene)

    #
    p = nx.spring_layout(G.subgraph([n for n in [i for i in nx.all_neighbors(mG, gene)] if n in tfs] + [gene]))
    nx.draw(
        mG.subgraph([n for n in [i for i in nx.all_neighbors(mG, gene)] if n in tfs] + [gene]),
        pos=p
    )
    nx.draw_networkx_labels(
        mG.subgraph([n for n in [i for i in nx.all_neighbors(mG, gene)] if n in tfs] + [gene]),
        pos=p
    )
    nx.write_gexf(mG.subgraph([n for n in [i for i in nx.all_neighbors(mG, gene)] if n in tfs] + [gene]), "graph_modules.%s.mCLL.gexf" % gene)


# Visualize network of changes in big network of all interactions
# make network with differential
diff_reg_genes = diff[(diff > 3) | (diff < -2)].index.tolist()

diff_reg_genes_p = diff[(diff > 3)]
diff_reg_genes_n = diff[(diff < -2)]
fc_std_p = (diff_reg_genes_p - diff_reg_genes_p.min()) / (diff_reg_genes_p.max() - diff_reg_genes_p.min()) + 1
fc_std_n = (diff_reg_genes_n - diff_reg_genes_n.min()) / (diff_reg_genes_n.max() - diff_reg_genes_n.min())
fc_std_n -= 2
fc_std = fc_std_n.append(fc_std_p)

# get neighbours
p = [item for sublist in [G.predecessors(i) for i in diff_reg_genes if i in G.nodes()] for item in sublist] + diff_reg_genes
# s = [item for sublist in [G.successors(i) for i in diff_reg_genes if i in G.nodes()] for item in sublist]

# get subgraph
g = G.subgraph(p)

# add fold-change attribute
fc_dict = {node: float(change) for node, change in fc_std.to_dict().items() if change < 10 and node in g.nodes()}
absfc_dict = {node: abs(float(change)) for node, change in fc_std.to_dict().items() if change < 10 and node in g.nodes()}

for node in g.nodes():
    if node not in fc_dict.keys():
        fc_dict[node] = 0
        absfc_dict[node] = 0

nx.set_node_attributes(g, 'fold_change', fc_dict)
nx.set_node_attributes(g, 'absfold_change', absfc_dict)
nx.write_gexf(g, "netx.fold_change.string.gexf")


# Test enrichment of oncogenes/tumour suppressors in lists
from scipy.stats import fisher_exact
# read in tumour suppresor list
ts = pd.read_table("tumour_suppressors.txt", header=None)[0].tolist()
onc = pd.read_table("oncogenes.txt", header=None)[0].tolist()

# store results
df2 = pd.DataFrame()
# read in tumour suppresor list
# uCLL
print("uCLL")
a = len([i for i in diff[diff > 1].index if i in onc])
b = len([i for i in diff[diff > 1].index if i not in onc])
c = len([i for i in diff[diff < 1].index if i in onc])
d = len([i for i in diff[diff < 1].index if i not in onc])
o, p = fisher_exact([[a, b], [c, d]])
df2 = df2.append(pd.Series([a / float(b), o, p, "onc", "u"]), ignore_index=True)

a = len([i for i in diff[diff > 1].index if i in ts])
b = len([i for i in diff[diff > 1].index if i not in ts])
c = len([i for i in diff[diff < 1].index if i in ts])
d = len([i for i in diff[diff < 1].index if i not in ts])
o, p = fisher_exact([[a, b], [c, d]])
df2 = df2.append(pd.Series([a / float(b), o, p, "ts", "u"]), ignore_index=True)

# mCLL
print("mCLL")
a = len([i for i in diff[diff < -1].index if i in onc])
b = len([i for i in diff[diff < -1].index if i not in onc])
c = len([i for i in diff[diff > -1].index if i in onc])
d = len([i for i in diff[diff > -1].index if i not in onc])
o, p = fisher_exact([[a, b], [c, d]])
df2 = df2.append(pd.Series([a / float(b), o, p, "onc", "m"]), ignore_index=True)

a = len([i for i in diff[diff < -1].index if i in ts])
b = len([i for i in diff[diff < -1].index if i not in ts])
c = len([i for i in diff[diff > -1].index if i in ts])
d = len([i for i in diff[diff > -1].index if i not in ts])
o, p = fisher_exact([[a, b], [c, d]])
df2 = df2.append(pd.Series([a / float(b), o, p, "ts", "m"]), ignore_index=True)
df2.columns = ["count", "odds", "pvalue", "comparison", "cluster"]

fig, axis = plt.subplots(1)
sns.barplot('cluster', 'count', hue='comparison', data=df2, ax=axis)
fig.savefig("oncogene_tumour-suppressors.enrichment.svg")


# Exclusive nodes
fc[fc.isnull()]

ex = fc[fc.isnull()].index.tolist()

orphans = dict()
for gene in ex:
    if gene in unode_degree.index:
        orphans[gene] = unode_degree.ix[gene]
    if gene in mnode_degree.index:
        orphans[gene] = -mnode_degree.ix[gene]

orphans = pd.Series(orphans)
orphans.sort()

sns.distplot(orphans)
plt.savefig("orphans.distplot.svg", bbox_inches="tight")

plt.scatter(orphans.rank(ascending=False, method='first'), orphans)
plt.savefig("orphans.rank_degree.svg", bbox_inches="tight")


# Exclusive interactions

# simpler ighv network comparison
graph_file = os.path.join(data_dir, "merged-samples_mutated_False.piq.TF-gene_interactions.filtered.tsv")
uG = create_graph(graph_file)

graph_file = os.path.join(data_dir, "merged-samples_mutated_True.piq.TF-gene_interactions.filtered.tsv")
mG = create_graph(graph_file)


uI = set(["-".join([i, j]) for i, j in uG.edges()])
mI = set(["-".join([i, j]) for i, j in mG.edges()])

uE = uI.difference(mI)
mE = mI.difference(uI)
