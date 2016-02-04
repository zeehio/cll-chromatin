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
from lifelines import AalenAdditiveFitter, CoxPHFitter
import patsy
from lifelines.utils import k_fold_cross_validation
from joblib import Parallel, delayed


# Set settings
pd.set_option("date_dayfirst", True)
sns.set_style("white")
sns.set_context("paper")
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc("font", family="sans-serif")
matplotlib.rc("font", serif="Helvetica Neue")
matplotlib.rc("text", usetex="false")


def survival_plot(clinical, fitter, fitter_name, feature, time, axis=None):
    """
    Plot survival/hazard of all patients regardless of trait and dependent of trait.
    """
    # all time (duration of patient following) regardless of trait
    T = [i.days / float(30) for i in clinical[time]]
    # events:
    # True for observed event (death);
    # else False (this includes death not observed; death by other causes)
    C = [True if i is not pd.NaT else False for i in clinical["patient_death_date"]]

    # drop index (to have syncronised numbers between lists T, C and the index of clinical)
    # this is because we drop a few rows from clinical before during data cleanup
    clinical2 = clinical.reset_index(drop=True)

    # plot
    if axis is None:
        fig, axis = plt.subplots(1)
        save = True
    else:
        save = False

    # Plot survival of all patients regardless of trait
    fitter.fit(T, event_observed=C, label="all patients")
    fitter.plot(ax=axis, show_censors=True)

    # For each type: subset, fit, plot
    # Filter patients which feature is nan
    x = clinical2[feature].unique()
    x = x[~np.array(map(pd.isnull, x))]

    # for each class plot curve
    for value in x:
        # get patients from class
        s = clinical2[clinical2[feature] == value].index.tolist()
        fitter.fit([T[i] for i in s], event_observed=[C[i] for i in s], label=str(value))
        fitter.plot(ax=axis, show_censors=True)
    if fitter_name == "survival":
        axis.set_ylim(0, 1.05)

    # Test pairwise differences
    p_values = list()
    # test each against all
    for a in x:
        a_ = clinical2[clinical2[feature] == a].index.tolist()
        b_ = clinical2.index.tolist()
        p = logrank_test(
            [T[i] for i in a_], [T[i] for i in b_],
            event_observed_A=[C[i] for i in a_],
            event_observed_B=[C[i] for i in b_]).p_value  # .print_summary()
        p_values.append(" vs ".join([str(a), "all"]) + ": %f" % p)
    # test each pairwise combination
    for a, b in itertools.combinations(x, 2):
        a_ = clinical2[clinical2[feature] == a].index.tolist()
        b_ = clinical2[clinical2[feature] == b].index.tolist()
        p = logrank_test(
            [T[i] for i in a_], [T[i] for i in b_],
            event_observed_A=[C[i] for i in a_],
            event_observed_B=[C[i] for i in b_]).p_value  # .print_summary()
        p_values.append(" vs ".join([str(a), str(b)]) + ": %f" % p)

    # Add p-values as anchored text
    try:  # problem with matplotlib < 1.4
        axis.add_artist(AnchoredText("\n".join(p_values), loc=8, frameon=False))
    except:
        pass

    axis.set_title("%s" % feature)
    axis.set_xlabel("time (months)")
    axis.set_ylabel(fitter_name)
    sns.despine()
    if save:
        fig.savefig(os.path.join(plots_dir, "%s_%s_since_%s.svg" % (feature, fitter_name, time)), bbox_inches="tight")


def test_model(model, clinical2):
    X = patsy.dmatrix(model + " -1", clinical2, return_type="dataframe").reset_index(drop=True)
    if X.empty:
        return [None for _ in range(1, 11)]
    X["T"] = [i.days / float(30) for i in clinical2["duration_following"].ix[X.index]]
    X["E"] = [True if i is not pd.NaT else False for i in clinical2["patient_death_date"].ix[X.index]]

    r = list()
    for penalty in range(1, 11):
        aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
        r.append((penalty, k_fold_cross_validation(aaf, X, "T", event_col="E", k=10)))  # append = penalty, score
    return r


# plots output dir
plots_dir = os.path.join("results", "plots", "survival")
# Get clinical info
clinical = pd.read_csv(
    os.path.join("metadata", "annotation.csv"),
    parse_dates=[
        "patient_birth_date", "patient_death_date", "sample_collection_date",
        "diagnosis_date", "patient_last_checkup_date", "treatment_date"] + ["treatment_%i_date" % i for i in range(1, 5)],
    dayfirst=True)

# traits to investigate (used at several points)
traits = [
    "patient_gender",
    "diagnosis_disease",
    "diagnosis_stage_binet",
    "diagnosis_stage_rai",
    "ighv_mutation_status",
    "CD38_positive",
    "ZAP70_positive",
    "ZAP70_monoallelic_methylation",
]
chrom_abrs = ["del13", "del11", "tri12", "del17"]
mutations = ["SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "TP53", "MYD88", "CHD2", "NFKIE"]

# annotate with chromosomal aberrations
for abr in chrom_abrs:
    # add column to clinical with bool info about said mutation
    clinical[abr] = clinical["mutations"].apply(lambda x: 1 if abr in str(x) else 0)

# annotate with mutations
for mutation in mutations:
    # add column to clinical with bool info about said mutation
    clinical[mutation] = clinical["mutations"].apply(lambda x: 1 if mutation in str(x) else pd.np.nan)

# Get one record per patient
cols = [
    "patient_birth_date", "patient_death_date", "sample_collection_date",
    "diagnosis_date", "patient_last_checkup_date", "treatment_date"] + ["treatment_%i_date" % i for i in range(1, 5)]
clinical = clinical[cols + traits + chrom_abrs + mutations].drop_duplicates()


# Data cleanup
# remove left censors (patients without birth or diagnosis date)
clinical = clinical[~clinical["patient_birth_date"].isnull()]
clinical = clinical[~clinical["diagnosis_date"].isnull()]


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

features = traits + chrom_abrs + mutations

# Survival of each clinical feature
fig, axis = plt.subplots(3, 7, figsize=(50, 20))
axis = axis.flatten()
time = "duration_following"
for i, feature in enumerate(features):
    print(feature)
    survival_plot(clinical, KaplanMeierFitter(), "survival", feature, time, axis=axis[i])
fig.savefig(os.path.join(plots_dir, "all_traits.survival.svg"), bbox_inches="tight")


# Hazard of each clinical feature
fig, axis = plt.subplots(3, 7, figsize=(50, 20))
axis = axis.flatten()
time = "duration_following"
for i, feature in enumerate(features):
    print(feature)
    survival_plot(clinical, NelsonAalenFitter(nelson_aalen_smoothing=False), "hazard", feature, time, axis=axis[i])
fig.savefig(os.path.join(plots_dir, "all_traits.hazard.svg"), bbox_inches="tight")


# Regression using Aalen's additive model
features = traits + chrom_abrs + mutations
[features.pop(features.index(i)) for i in [
    "ZAP70_monoallelic_methylation", "diagnosis_stage_rai",
    "TP53", "SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "MYD88", "CHD2", "NFKIE"]]
combinations = [" + ".join(j) for i in range(1, len(features) + 1) for j in itertools.combinations(features, i)]

# Test all models
clinical2 = clinical.reset_index(drop=True)
performances = Parallel(n_jobs=-1, verbose=5)(delayed(test_model)(model, clinical2) for model in combinations)

validated = dict(zip(
    [re.sub(r" \+ ", "-", c) for c in combinations],
    sum(performances, []))
)
# filter out None (from models witough necessary info)
validated = {model: scores for model, scores in validated.items() if scores is not None}

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
fig.savefig(os.path.join(plots_dir, "model_selection_performance.svg"), bbox_inches="tight")

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
    X["T"] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
    X["E"] = [True if i is not pd.NaT else False for i in clinical["patient_death_date"].ix[X.index]]

    aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
    aaf.fit(X, "T", event_col="E")
    aaf.predict_survival_function(X).plot(ax=axis[0][m - 1], legend=False)
    aaf.predict_cumulative_hazard(X).plot(ax=axis[1][m - 1], legend=False)
    axis[0][m - 1].set_title(best_models[-m][0])
    axis[1][2].set_xlabel("Time since diagnosis")
    axis[0][0].set_ylabel("Survival")
    axis[1][0].set_ylabel("Hazard")
    # [i[m - 1].legend_.remove() for i in axis]
fig.savefig(os.path.join(plots_dir, "best_5models_predictions_all_patients.svg"), bbox_inches="tight")

# get best model
model = re.sub("-", " + ", best_models[-1][0])
# model = "patient_gender + ighv_mutation_status + ZAP70_positive + del11"
penalty = best_models[-1][1][0]
# penalty = 4

# regress with best model
X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")
X["T"] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
X["E"] = [True if i is not pd.NaT else False for i in clinical["patient_death_date"].ix[X.index]]

aaf = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
aaf.fit(X, "T", event_col="E")

# predict survival and hazard for all patients
fig, axis = plt.subplots(2, sharex=True, sharey=False)
aaf.predict_survival_function(X).plot(ax=axis[0], legend=False)
aaf.predict_cumulative_hazard(X).plot(ax=axis[1], legend=False)
axis[0].set_title("Best model - prediction for all patients")
axis[1].set_xlabel("Time since diagnosis")
axis[0].set_ylabel("Survival")
axis[1].set_ylabel("Hazard")
fig.savefig(os.path.join(plots_dir, "best_model_predictions_all_patients.svg"), bbox_inches="tight")


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
            s.append(np.interp(t.days / float(30), survival.index, survival[i]))
            h.append(np.interp(t.days / float(30), hazard.index, hazard[i]))
            continue
    s.append(np.nan)
    h.append(np.nan)

clinical["predicted_survival"] = s
clinical["predicted_hazard"] = h

# investigate
fig, axis = plt.subplots(2)
axis[0].scatter(
    clinical["duration_collection"].astype("timedelta64[M]"),
    clinical["predicted_survival"])
axis[1].scatter(
    clinical["duration_collection"].astype("timedelta64[M]"),
    clinical["predicted_hazard"])
fig.savefig(os.path.join(plots_dir, "predicted_at_collection_time.svg"), bbox_inches="tight")

# export
clinical[["patient_id", "sample_id", "predicted_survival", "predicted_hazard", "duration_collection"]].to_csv(
    os.path.join("data", "survival_hazard_predictions.csv"),
    index=False
)

# Regression on probability of survival, time-to-first-treatment, time-to-death
# first, get earlier date of treatment for each patient
cols = ["treatment_%i_date" % i for i in range(1, 5)]
clinical["first_treament_date"] = clinical[cols].apply(lambda x: min(x.dropna()) if not x.dropna().empty else np.nan, axis=1)
clinical["time_to_treatment"] = clinical["first_treament_date"] - clinical["diagnosis_date"]

# first, get earlier date of treatment for each patient
clinical["time_to_death"] = clinical["patient_death_date"] - clinical["sample_collection_date"]

# plot distribution
sns.distplot(pd.Series([i.days / 30. if i is not pd.NaT else pd.np.nan for i in clinical.time_to_treatment]).dropna(), bins=50)
sns.distplot(pd.Series([i.days / 30. if i is not pd.NaT else pd.np.nan for i in clinical.time_to_death]).dropna(), bins=50)

clinical[["patient_id", "sample_id", "predicted_survival", "predicted_hazard", "duration_collection", "time_to_treatment"]].to_csv(
    os.path.join("data", "survival_hazard_predictions_time2treatment.csv"),
    index=False
)


#


# Regression using Cox's additive hazard model
features = traits + chrom_abrs + mutations
[features.pop(features.index(i)) for i in [
    "patient_gender", "diagnosis_disease", "diagnosis_stage_binet",
    "ZAP70_monoallelic_methylation", "diagnosis_stage_rai",
    "TP53",
    "SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "MYD88", "CHD2", "NFKIE"]]

# features = ["ighv_mutation_status", "CD38_positive", "ZAP70_positive", "del13"]
# regress with complex model
X = patsy.dmatrix(" + ".join(features) + " -1", clinical, return_type="dataframe")
X["T"] = [i.days / float(30) for i in clinical["duration_following"].ix[X.index]]
X["E"] = [1 if i is not pd.NaT else 0 for i in clinical["patient_death_date"].ix[X.index]]
for col in X.columns:
    X[col] = X[col]# + (1. / 10000000000)
cf = CoxPHFitter()
cf.fit(X, "T", event_col="E")

# check coeficients
cf.print_summary()

# predict survival and hazard for all patients
fig, axis = plt.subplots(2, sharex=True, sharey=False)
cf.predict_survival_function(X).plot(ax=axis[0], legend=False)
cf.predict_cumulative_hazard(X).plot(ax=axis[1], legend=False)
axis[0].set_title("Best model - prediction for all patients")
axis[1].set_xlabel("Time since diagnosis")
axis[0].set_ylabel("Survival")
axis[1].set_ylabel("Hazard")
fig.savefig(os.path.join(plots_dir, "best_model_predictions_all_patients.svg"), bbox_inches="tight")
