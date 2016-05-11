#!/usr/bin/env python

"""

"""

import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
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

# list of sex chromosome genes
sex_genes = pd.read_csv("sex_genes.csv")['name'].unique().tolist()

# get gene expression
expression_genes = pd.read_csv(os.path.join(data_dir, "cll_expression_matrix.log2.csv"))
# map ensembl to genesymbol
# add gene name and ensemble_gene_id
ensembl_gtn = pd.read_table(os.path.join(data_dir, "external", "ensemblToGeneName.txt"), header=None)
ensembl_gtn.columns = ['ensembl_transcript_id', 'gene_name']
ensembl_gtp = pd.read_table(os.path.join(data_dir, "external", "ensGtp.txt"), header=None)[[0, 1]]
ensembl_gtp.columns = ['ensembl_gene_id', 'ensembl_transcript_id']
gene_annotation = pd.merge(ensembl_gtn, ensembl_gtp, how="left")[["gene_name", "external", "ensembl_gene_id"]].drop_duplicates()

# simpler ighv network comparison
new_graph_file = os.path.join("merged-samples_all_all.piq.TF-gene_interactions.filtered.tsv")
tfs = pd.read_table(new_graph_file)['TF'].drop_duplicates()
G = create_graph(new_graph_file)

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
u = pd.read_csv("merged-samples_mutated_False.piq.TF-gene_interactions.tsv", sep="\t")
u = u[u['interaction_score'] > 1]
u.to_csv("merged-samples_mutated_False.piq.TF-gene_interactions.filtered.tsv", sep="\t", index=False)
m = pd.read_csv("merged-samples_mutated_True.piq.TF-gene_interactions.tsv", sep="\t")
m = m[m['interaction_score'] > 1]
m.to_csv("merged-samples_mutated_True.piq.TF-gene_interactions.filtered.tsv", sep="\t", index=False)

graph_file = os.path.join("merged-samples_mutated_False.piq.TF-gene_interactions.filtered.tsv")
uG = create_graph(graph_file)

graph_file = os.path.join("merged-samples_mutated_True.piq.TF-gene_interactions.filtered.tsv")
mG = create_graph(graph_file)

# remove sex genes
uG.remove_nodes_from(sex_genes)
mG.remove_nodes_from(sex_genes)

# remove outgoing nodes of non expressed TFs
# uG.remove_edges_from(uG.out_edges(expressed_tfs))
# mG.remove_edges_from(mG.out_edges(expressed_tfs))

d = pd.Series(uG.degree())
g = uG.subgraph(d[d > 20].index)
# g = uG
nx.set_node_attributes(g, 'abs_degree', {n: int(d.to_dict()[n]) for n in g.nodes()})
nx.write_gexf(g, "merged-samples_mutated_False.piq.TF-gene_interactions.filtered.filtered_gene_exp.top.gexf")
nx.write_gexf(uG, "merged-samples_mutated_False.piq.TF-gene_interactions.filtered.filtered_gene_exp.gexf")

d = pd.Series(mG.degree())
g = mG.subgraph(d[d > 20].index)
nx.set_node_attributes(g, 'abs_degree', {n: int(d.to_dict()[n]) for n in g.nodes()})
nx.write_gexf(g, "merged-samples_mutated_True.piq.TF-gene_interactions.filtered.filtered_gene_exp.top.gexf")
nx.write_gexf(mG, "merged-samples_mutated_True.piq.TF-gene_interactions.filtered.filtered_gene_exp.gexf")


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

# write differential results
diff.to_csv(os.path.join("merged-samples.TF-gene_interactions.IGHV_differential.csv"), index=True)


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
cllgenes = pd.read_csv(os.path.join("metadata", "gene_lists", "bcell_cll_genelist.tsv"), sep="\t")["gene_name"]
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


# CD19-CLL comparison
# Load up joint CLL network
new_graph_file = os.path.join("/home/arendeiro/cll-patients/netwx/merged-samples_all_all.piq.TF-gene_interactions.filtered.tsv")
G = create_graph(new_graph_file)

# Load up CD19 network
out_dir = os.path.abspath(os.path.join(data_dir, "external", "CD19_DNase"))
foots_dir = os.path.join(out_dir, "footprints")
label = "CD19_DNase"
graph_file = os.path.join(foots_dir, label + ".piq.TF-gene_interactions.filtered.tsv")
cd19G = create_graph(graph_file)

# get node degree
cll_node_degree = pd.Series(G.degree(), name="degree")
cd19_node_degree = pd.Series(cd19G.degree(), name="degree")

# Compare nodes
df = pd.DataFrame([cll_node_degree, cd19_node_degree]).T
df.columns = ['cll', 'cd19']
# normalize
df = df.apply(lambda x: (x - x.mean() / x.std()), axis=0)

# handle nodes not detected in both
# df = np.log2(1 + df.fillna(0))  # assume no binding
df = np.log2(1 + df.dropna())  # exclude genes

plt.scatter(df['cll'], df['cd19'])
plt.plot([0, 14], [0, 14], '--')

# let's look at the difference now
diff = df['cll'] - df['cd19']


# Describe structure of CD19 and CLL networks
descG = describe_graph(G)
descGcd19 = describe_graph(cd19G)
plot_graph_attributes()

descG[0]["group"] = "CLL"
descGcd19[0]["group"] = "CD19"

desc = pd.DataFrame([descG[0], descGcd19[0]])
melted = pd.melt(desc, ['group'])

g = sns.FacetGrid(melted, col="variable", col_wrap=4, legend_out=True, margin_titles=True, sharey=False)
g.map(sns.barplot, "group", "value")
sns.despine()
plt.savefig(os.path.join(
    plots_dir, "cll-cd19_network.structure_comparison.net.svg"),
    bbox_inches='tight')


# CD19+ degree of all nodes

import matplotlib
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white")
sns.set_style({'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica']})
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex=False)

tfs = pd.read_table(new_graph_file)['TF'].drop_duplicates().tolist()
d = pd.Series(cd19G.degree())
d.sort()
fig, axis = plt.subplots(1)
colors = ["#d7191c" if i in tfs else "#669fb8" for i in d.index]
axis.scatter(d.rank(ascending=False, method='first'), np.log2(d), linewidth=0, color=colors)
sns.despine()
fig.savefig(os.path.join(
    plots_dir, "cd19_network.node_degree.svg"), bbox_inches="tight")

d = pd.Series(G.degree())
d.sort()
fig, axis = plt.subplots(1)
colors = ["#d7191c" if i in tfs else "#669fb8" for i in d.index]
axis.scatter(d.rank(ascending=False, method='first'), np.log2(d), linewidth=0, color=colors)
sns.despine()
fig.savefig(os.path.join(
    plots_dir, "cll_network.node_degree.svg"), bbox_inches="tight")
