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
from lifelines.utils import k_fold_cross_validation
from joblib import Parallel, delayed


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
    T = [i.days / float(30) for i in clinical[time]]  # duration of patient following
    # events:
    # True for observed event (death);
    # else False (this includes death not observed; death by other causes)
    C = [True if i is not pd.NaT else False for i in clinical["patient_death_date"]]

    fitter = fitter()

    # plot
    ax = plt.subplot(111)

    fitter.fit(T, event_observed=C, label="all patients")
    fitter.plot(ax=ax, show_censors=True)

    # for each type: subset, fit, plot
    # Filter feature types which are nan
    x = clinical[feature].unique()
    x = x[~np.array(map(pd.isnull, x))]

    # diagnosis stage remove secondary CLL patients
    if feature == "diagnosis_stage":
        x = filter(lambda i: False if i in ['MBL', 'SLL', '?'] else True, x)

    # for each class plot curve
    clinical2 = clinical.reset_index(drop=True)
    for value in x:
        # get patients from class
        s = clinical2[clinical2[feature] == value].index.tolist()
        fitter.fit([T[i] for i in s], event_observed=[C[i] for i in s], label=str(value))
        fitter.plot(ax=ax, show_censors=True)
    if fitter_name == "survival":
        ax.set_ylim(0, 1.05)

    # Test pairwise differences
    p_values = list()
    for a, b in itertools.combinations(x, 2):
        a_ = clinical2[clinical2[feature] == a].index.tolist()
        b_ = clinical2[clinical2[feature] == b].index.tolist()
        p = logrank_test(
            [T[i] for i in a_], [T[i] for i in b_],
            event_observed_A=[C[i] for i in a_],
            event_observed_B=[C[i] for i in b_]).p_value  # .print_summary()
        p_values.append(" - ".join([str(a), str(b)]) + ": %f" % p)

    # Add p-values as anchored text
    # ax.add_artist(AnchoredText("\n".join(p_values), loc=8, frameon=False))

    ax.set_title("%s - %s since %s" % (feature, fitter_name, name))
    plt.savefig(os.path.join(plots_dir, "survival", "%s_%s_since_%s.svg" % (feature, fitter_name, name)), bbox_inches="tight")
    plt.close("all")


def test_model(model):
    X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")
    X['T'] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
    X['E'] = [True if i is not pd.np.nan else False for i in clinical["patient_death_date"].ix[X.index]]

    r = list()
    for penalty in range(1, 11):
        aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
        r.append((penalty, k_fold_cross_validation(aaf, X, 'T', event_col='E', k=5)))
    return r


#
plots_dir = os.path.join("results", "plots")
# Get clinical info
clinical = pd.read_csv(
    os.path.join("metadata", "clinical_annotation.csv"),
    parse_dates=[
        "patient_birth_date", "patient_death_date", "sample_collection_date",
        "diagnosis_date", "patient_last_checkup_date"] + ["treatment_%i_date" % i for i in range(1, 5)],
    dayfirst=True)

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
subs = clinical[cols].apply(lambda x: max(x.dropna()), axis=1).groupby(clinical.index).max()
no_checkup = clinical[clinical["patient_last_checkup_date"].isnull()]

for i in no_checkup.index:
    if subs[i] is not pd.NaT:
        clinical.loc[i, "patient_last_checkup_date"] = subs[i]

clinical["duration_following"] = clinical["patient_last_checkup_date"] - clinical["diagnosis_date"]
clinical["duration_life"] = clinical["patient_last_checkup_date"] - clinical["patient_birth_date"]


# For time since birth and time since diagnosis
times = {
    # "birth": "duration_life",
    "diagnosis": "duration_following"}

features = [
    "patient_gender",
    "diagnosis_disease",
    "diagnosis_stage",
    "igvh_mutation_status",
    "CD38_positive",
    "ZAP70_positive",
    "ZAP70_monoallelic_methylation",
]
muts = ['del13', 'del11', 'tri12']
muts += ['SF3B1', 'ATM', 'NOTCH1', 'BIRC3', 'BCL2', 'TP53', 'MYD88', 'CHD2', 'NFKIE']
features += muts

# For each time since (birth, diagnosis)
for name, time in times.items():
    # For each clinical feature
    for feature in features:
        print(name, feature)
        if len(feature) < 7:
            # add column to clinical with bool info about said mutation
            clinical[feature] = clinical.mutations.apply(lambda x: True if feature in str(x) else False)
        # else:
        #     clinical[feature] = clinical[feature].replace(1.0, False).replace(2.0, True)

        survival_plot(clinical, KaplanMeierFitter, "survival", feature, time)
        survival_plot(clinical, NelsonAalenFitter, "hazard", feature, time)


# Regression using Aalen's additive model
[features.pop(features.index(i)) for i in ['SF3B1', 'ATM', 'NOTCH1', 'BIRC3', 'BCL2', 'MYD88', 'CHD2', 'NFKIE']]
combinations = [" + ".join(j) for i in range(1, len(features) + 1) for j in itertools.combinations(features, i)]

# Test all models
performances = Parallel(n_jobs=-1)(delayed(test_model)(model) for model in combinations)

validated = dict(zip(
    [re.sub(r" \+ ", "-", c) for c in combinations],
    sum(performances, []))
)

# model concordance distribution
fig, axis = plt.subplots(3, sharex=False, sharey=False, figsize=(15, 8))
sns.distplot([np.mean(x[1]) for x in validated.values()], ax=axis[0])
# number of elements in model
axis[1].scatter([len(x.split("-")) for x in validated.keys()], [np.mean(x[1]) for x in validated.values()])
# penalty per model
axis[2].scatter([x[0] for x in validated.values()], [np.mean(x[1]) for x in validated.values()])
axis[0].set_xlabel("Mean model concordance")
axis[0].set_ylabel("Density")
axis[1].set_xlabel("Number of model terms")
axis[1].set_ylabel("Mean concordance")
axis[2].set_xlabel("Co-variate penalty")
axis[2].set_ylabel("Mean concordance")
fig.savefig(os.path.join(plots_dir, "survival", "model_selction_performance.svg"), bbox_inches="tight")


# sort by performance
best_models = sorted(validated.items(), key=lambda x: np.mean(x[1][1]))

# Plot 5 best models
fig, axis = plt.subplots(2, 5, sharex=False, sharey=False, figsize=(40, 15))
for m in range(1, 6):
    # get best model
    model = re.sub("-", " + ", best_models[-m][0])
    penalty = best_models[-m][1][0]

    # regress with best model
    X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")
    X['T'] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
    X['E'] = [True if i is not pd.np.nan else False for i in clinical["patient_death_date"].ix[X.index]]

    aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
    aaf.fit(X, 'T', event_col='E')
    aaf.predict_survival_function(X).plot(ax=axis[0][m - 1], legend=False)
    aaf.predict_cumulative_hazard(X).plot(ax=axis[1][m - 1], legend=False)
    axis[0][m - 1].set_title(best_models[-m][0])
    axis[1][2].set_xlabel("Time since diagnosis")
    axis[0][0].set_ylabel("Survival")
    axis[1][0].set_ylabel("Hazard")
    # [i[m - 1].legend_.remove() for i in axis]
fig.savefig(os.path.join(plots_dir, "survival", "best_5models_predictions_all_patients.svg"), bbox_inches="tight")

# get best model
model = re.sub("-", " + ", best_models[-1][0])
penalty = best_models[-1][1][0]

# regress with best model
X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")
X['T'] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
X['E'] = [True if i is not pd.np.nan else False for i in clinical["patient_death_date"].ix[X.index]]

aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
aaf.fit(X, 'T', event_col='E')

# predict survival and hazard for all patients
fig, axis = plt.subplots(2, sharex=True, sharey=False)
aaf.predict_survival_function(X).plot(ax=axis[0], legend=False)
aaf.predict_cumulative_hazard(X).plot(ax=axis[1], legend=False)
axis[0].set_title("Best model - prediction for all patients")
axis[1].set_xlabel("Time since diagnosis")
axis[0].set_ylabel("Survival")
axis[1].set_ylabel("Hazard")
fig.savefig(os.path.join(plots_dir, "survival", "best_model_predictions_all_patients.svg"), bbox_inches="tight")


# Predict survival probability and hazard of each sample
# at the time of sample collection
# based on best model

survival = aaf.predict_survival_function(X)
hazard = aaf.predict_cumulative_hazard(X)

clinical["duration_collection"] = clinical["sample_collection_date"] - clinical["diagnosis_date"]

s = list()
h = list()
for i in clinical.index:
    if i in survival.columns:
        t = clinical.ix[i]["duration_collection"]
        if t is not pd.NaT:
            s.append(np.interp(t.days, survival.index, survival[i]))
            h.append(np.interp(t.days, hazard.index, hazard[i]))
            continue
    s.append(np.nan)
    h.append(np.nan)

clinical['predicted_survival'] = s
clinical['predicted_hazard'] = h

# export
clinical[["patient_id", "sample_id", "predicted_survival", "predicted_hazard"]].to_csv(
    os.path.join("data", "survival_hazard_predictions.csv")
)
