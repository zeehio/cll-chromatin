#!/usr/bin/env python

"""
This script makes plots to describe the cohort.
"""

import os
# from pipelines.models import Project
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


sns.set_style('whitegrid')
sns.set_context('paper')


# # start project
# prj = Project('cll-patients')
# prj.addSampleSheet("metadata/sequencing_sample_annotation.csv")

plotsDir = os.path.join('results', 'plots')

# Get clinical info
clinical = pd.read_csv(os.path.join('metadata', 'clinical_annotation.csv'))
clinical.index = clinical['patient_id']

# replace dubious/inconsistent
mut_counts = clinical['mutations'].apply(lambda x: x.replace("?", "") if type(x) is str else np.nan)
clinical['mutations'] = mut_counts.apply(lambda x: x.replace("TP53loss", "TP53") if type(x) is str else np.nan)

# Make categorical numeric
binars = ['patient_gender', 'diagnosis_disease', 'treated', 'sample_shippment_batch']
mapping = {'F': 0, 'M': 1, 'CLL': 0, 'MBL': 1, 'SLL': 2, 'Y': 0, 'N': 1, 'A': 0, 'B': 1, 'C': 2}
for bina in binars:
    g = list()
    for x in clinical[bina]:
        if x in mapping.keys():
            g.append(mapping[x])
        else:
            g.append(np.nan)
    clinical[bina] = g

clinical['mutations_notch1'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('NOTCH1')]
clinical['mutations_sf3b1'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('SF3B1')]
clinical['mutations_del13'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('del13')]
clinical['mutations_del11q'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('del11q')]
clinical['mutations_tp53'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('TP53')]
clinical['mutations_tri12'] = [1 if i else np.nan for i in clinical['mutations'].str.contains('tri12')]

# mutation count per patient
clinical['mutations_count'] = clinical['mutations'].apply(lambda x: len(x.split(',')) if type(x) is str else 0)

# Parse dates
dates = [
    'patient_birth_date', 'patient_death_date', 'diagnosis_date',
    'patient_last_checkup_date', 'sample_collection_date', 'sample_processing_date']
for col in dates:
    clinical[col] = pd.to_datetime(clinical[col])

# Get age at diagnosis
clinical['patient_age_at_diagnosis'] = (clinical['diagnosis_date'] - clinical['patient_birth_date']).astype('timedelta64[D]') / 365
# Get age at time of sample
clinical['patient_age_at_collection'] = (clinical['sample_collection_date'] - clinical['patient_birth_date']).astype('timedelta64[D]') / 365

# Dead or alive
clinical['patient_alive'] = [1 if type(x) is pd.tslib.Timestamp else 0 for x in clinical['patient_death_date']]


# HEATSTRIP PLOT
attrs = [
    'igvh_homoology', 'igvh_mutation_status', 'diagnosis_disease', 'patient_gender', 'patient_alive', 'treated', 'lymp_count',
    'patient_age_at_diagnosis', 'patient_age_at_collection',
    'mutations_count', 'mutations_notch1', 'mutations_sf3b1', 'mutations_del13', 'mutations_del11q', 'mutations_tp53', 'mutations_tri12',
    'sample_viability', 'sample_shippment_batch',
]

# get only attrs
df = clinical[attrs]
df = df.sort(columns=['igvh_mutation_status', 'diagnosis_disease', 'treated', 'patient_gender', 'igvh_homoology'])

f, axarr = plt.subplots(len(attrs), 1, figsize=(20, 12))
for i, attr in enumerate(attrs):
    if i != len(attrs) - 1:
        if attr in ['igvh_mutation_status', 'igvh_homoology']:
            sns.heatmap(pd.DataFrame(df[attr]).T, cmap=plt.get_cmap('summer'), ax=axarr[i], xticklabels=False)
        elif attr == 'diagnosis_disease':
            sns.heatmap(pd.DataFrame(df[attr]).T, cmap=plt.get_cmap('Paired'), ax=axarr[i], xticklabels=False)
        elif attr in ['treated', 'sample_shippment_batch', 'patient_alive']:
            sns.heatmap(pd.DataFrame(df[attr]).T, cmap=plt.get_cmap('YlGn_r'), ax=axarr[i], xticklabels=False)
        elif attr in ['patient_gender']:
            cmap = plt.get_cmap('RdYlBu')
            cmap.set_gamma(2)
            sns.heatmap(pd.DataFrame(df[attr]).T, cmap=cmap, ax=axarr[i], xticklabels=False)
        else:
            sns.heatmap(pd.DataFrame(df[attr]).T, cmap=plt.get_cmap('gray_r'), ax=axarr[i], xticklabels=False)
    else:
        sns.heatmap(pd.DataFrame(df[attr]).T, cmap=plt.get_cmap('gray_r'), ax=axarr[i], xticklabels=True)
        axarr[i].set_xticklabels(pd.DataFrame(df[attr]).T.columns, rotation=45)
    axarr[i].set_yticklabels([attr], rotation=0)

plt.savefig(os.path.join(plotsDir, 'cohort.pdf'), bbox_inches='tight')

# save plot data
df.to_csv(os.path.join("data", "cohort.csv"), index=False)


# EVALUATE MUTATIONS IN COHORT
# Plot mutations
# count
mut_counts = Counter([i for x in clinical['mutations'].apply(lambda x: x.split(",") if type(x) == str else []) for i in x])
mut_counts.pop(" ")
mut_counts = pd.DataFrame([mut_counts.keys(), mut_counts.values()]).T
mut_counts.sort(columns=[1], inplace=True, ascending=False)

# plot
fig, axis = plt.subplots(figsize=(8, 6))
sns.barplot(data=mut_counts, x=0, y=1, ax=axis)
axis.set_title("Mutated genes and genomic aberrations in cohort")
axis.set_xlabel("genes")
axis.set_ylabel("frequency")
plt.savefig(os.path.join(plotsDir, 'cohort_mutations.pdf'), bbox_inches='tight')

# save plot data
df.to_csv(os.path.join("data", "cohort_mutations.csv"), index=False)


# Kapplar-Mayer plot
# drop patients without dob
df = clinical[clinical['patient_birth_date'].notnull()]

df3 = pd.DataFrame(columns=['years', 'survival', 'igvh_mutation_status'])

for status in [1, 2]:
    df2 = df[df['igvh_mutation_status'] == status]

    # get death time after diagnosis
    death_diagnosis = df2['patient_death_date'] - df2['diagnosis_date']
    death_diagnosis = death_diagnosis.drop_duplicates().dropna()
    death_diagnosis = np.array([x.days for x in death_diagnosis]) / 365.0

    # add start point
    death_diagnosis = np.append(0, death_diagnosis)

    # for each time of death reduce survival by 1
    points = dict()
    for i, time in enumerate(sorted(death_diagnosis)):
        if time == 0:
            points[time] = len(df2)
        else:
            points[time] = points[sorted(death_diagnosis)[i - 1]] - 1

    # transform to percentage
    points = {time: (float(count) / len(df2)) * 100 for time, count in points.items()}

    # make df
    d = pd.DataFrame([points.keys(), points.values()], index=['years', 'survival']).T
    d['igvh_mutation_status'] = "mCLL" if status == 1 else "uCLL"
    df3 = pd.concat([df3, d])

# plot
fig = sns.lmplot(x='years', y='survival', data=df3, hue="igvh_mutation_status", fit_reg=False, scatter_kws={"s": 50})
fig.ax.set_xlabel("time after diagnosis (years)")
fig.ax.set_ylabel("% survival")
fig.ax.set_xlim(0, 40)
fig.ax.set_ylim(0, 100)
fig.fig.set_title("Kapplar-Mayer curve of cohort dependent on IGVH mutation status")
fig.fig.savefig(os.path.join(plotsDir, 'cohort_kapplar-mayer.pdf'), bbox_inches='tight')


# other option
# df3 = df3.sort("years")
# grid = sns.FacetGrid(df3, hue="igvh_mutation_status")
# grid.map(plt.plot, "years", "survival", marker="o")

# save plot data
df.to_csv(os.path.join("data", "cohort_kapplar-mayer.csv"), index=False)
