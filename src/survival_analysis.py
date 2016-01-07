#!/usr/bin/env python

"""
Survival analysis on the cohort
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines import NelsonAalenFitter
import itertools
from lifelines.statistics import logrank_test
import re
from lifelines import AalenAdditiveFitter
import patsy


# Set settings
pd.set_option("date_dayfirst", True)
sns.set_style("whitegrid")
sns.set_context("paper")
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Helvetica Neue")
matplotlib.rc("text", usetex="false")


def survival_plot(clinical, fitter, fitter_name, feature, time):
    T = [i.days / 365 for i in clinical[time]]  # duration of patient following
    # events:
    # True for observed event (death);
    # else False (this includes death not observed; death by other causes)
    C = [True if i is not pd.np.nan else False for i in clinical["patient_death_date"]]

    fitter = fitter()

    # plot
    ax = plt.subplot(111)

    fitter.fit(T, event_observed=C, label="all patients")
    fitter.plot(ax=ax, show_censors=True)

    # for each type: subset, fit, plot
    # clinical = clinical.reset_index(drop=True)
    x = clinical[feature].unique()
    x = x[~np.array(map(pd.isnull, x))]

    # diagnosis stage remove secondary CLL patients
    if feature == "diagnosis_stage":
        x = filter(lambda i: False if i in ['MBL', 'SLL', '?'] else True, x)
    for value in x:
        s = clinical[clinical[feature] == value].index.tolist()
        fitter.fit([T[i] for i in s], event_observed=[C[i] for i in s], label=str(value))
        fitter.plot(ax=ax, show_censors=True)
    if fitter_name == "survival":
        ax.set_ylim(0, 1.05)

    # Test pairwise differences
    p_values = list()
    for a, b in itertools.combinations(x, 2):
        a_ = clinical[clinical[feature] == a].index.tolist()
        b_ = clinical[clinical[feature] == b].index.tolist()
        p = logrank_test(
            [T[i] for i in a_], [T[i] for i in b_],
            event_observed_A=[C[i] for i in a_],
            event_observed_B=[C[i] for i in b_]).p_value  # .print_summary()
        p_values.append(" - ".join([str(a), str(b)]) + ": %f" % p)

    # Add p-values as anchored text
    ax.add_artist(AnchoredText("\n".join(p_values), loc=8, frameon=False))

    ax.set_title("%s - %s since %s" % (feature, fitter_name, name))
    plt.savefig(os.path.join(plots_dir, "%s_%s_since_%s.svg" % (feature, fitter_name, name)), bbox_inches="tight")
    plt.close("all")


# Get clinical info
clinical = pd.read_csv(os.path.join("metadata", "clinical_annotation.csv"))
clinical.index = clinical["patient_id"]

plots_dir = os.path.join("results", "plots")

# Data cleanup
# remove left censors (patients without birth or diagnosis date)
clinical = clinical[~clinical["patient_birth_date"].isnull()]
clinical = clinical[~clinical["diagnosis_date"].isnull()]
# remove weird stuff from "diagnosis_stage" column
clinical.diagnosis_stage.replace(to_replace=['MBL', 'SLL', '?'], value=pd.np.nan, inplace=True)
# assume all cases are non mono-allelic methylation except when explicitely said so
clinical.ZAP70_monoallelic_methylation.fillna("N", inplace=True)


# Lifelines analysis
# Get duration of patient observation
# but first, let"s replace missing "patient_last_checkup_date" with last date of treatment, death or collection
cols = ["patient_death_date", "sample_collection_date", "diagnosis_date"] + ["treatment_%i_date" % i for i in range(1, 5)]
subs = clinical[cols].apply(pd.to_datetime).apply(lambda x: max(x.dropna()), axis=1).groupby(clinical.index).max()
no_checkup = clinical[clinical["patient_last_checkup_date"].isnull()]

for i in no_checkup.index:
    if subs[i] is not pd.NaT:
        clinical.loc[i, "patient_last_checkup_date"] = subs[i]

clinical["duration_following"] = pd.to_datetime(clinical["patient_last_checkup_date"]) - pd.to_datetime(clinical["diagnosis_date"])
clinical["duration_life"] = pd.to_datetime(clinical["patient_last_checkup_date"]) - pd.to_datetime(clinical["patient_birth_date"])


# For time since birth and time since diagnosis
times = {"birth": "duration_life", "diagnosis": "duration_following"}

features = [
    "patient_gender",
    "diagnosis_disease",
    "diagnosis_stage",
    "igvh_mutation_status",
    "CD38_positive",
    "ZAP70_positive",
    "ZAP70_monoallelic_methylation",
]

# For each time since (birth, diagnosis)
for name, time in times.items():
    # For each clinical feature
    for feature in features:
        print(name, feature, time)
        survival_plot(clinical, KaplanMeierFitter, "survival", feature, time)
        survival_plot(clinical, NelsonAalenFitter, "hazard", feature, time)

# Regression using Aalen's additive model
combinations = [" + ".join(j) for i in range(1, len(features) + 1) for j in itertools.combinations(features, i)]

models = zip(
    [re.sub(r" \+ ", "-", c) for c in combinations],
    combinations
)

for model_name, model in models:
    # Create matrix
    X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")

    # Add the survival data
    # and since patsy removes rows with nan in variables used for matrix,
    # we have to remove those as well from the survival data (fortunately original index is kept in X yay!)
    X['T'] = [i.days / 365 for i in clinical["duration_following"].ix[X.index]]
    X['E'] = [True if i is not pd.np.nan else False for i in clinical["patient_death_date"].ix[X.index]]

    aaf = AalenAdditiveFitter(coef_penalizer=1.0, fit_intercept=True)
    aaf.fit(X, 'T', event_col='E')
