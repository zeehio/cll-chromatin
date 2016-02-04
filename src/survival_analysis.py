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


def life_duration(x, kind="diagnosis_date"):
    if x["patient_death_date"] is pd.NaT:
        return x["patient_last_checkup_date"] - x[kind]
    else:
        return x["patient_death_date"] - x[kind]


def survival_plot(clinical, fitter, fitter_name, feature, axis=None):
    """
    Plot survival/hazard of all patients regardless of trait and dependent of trait.
    """
    # duration of life
    T = [i.days / 30. for i in clinical["duration"]]
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
        fig.savefig(os.path.join(plots_dir, "%s_%s.svg" % (feature, fitter_name)), bbox_inches="tight")


def test_model(model, fitter_name, df):
    X = patsy.dmatrix(model + " -1", df, return_type="dataframe").reset_index(drop=True)
    if X.empty:
        return [None for _ in range(1, 11)]
    X["T"] = [i.days / 30. for i in df["duration"].ix[X.index]]
    X["E"] = [True if i is not pd.NaT else False for i in df["patient_death_date"].ix[X.index]]

    scores = list()
    for penalty in range(1, 11):
        if fitter_name == "aalen":
            fitter = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
        else:
            fitter = CoxPHFitter(penalizer=penalty * 0.1)
        score = k_fold_cross_validation(fitter, X, "T", event_col="E", k=10)
        if score > 0:
            scores.append((penalty, score))  # append = penalty, score
    return scores


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
    "patient_id", "patient_birth_date", "patient_death_date", "sample_collection_date",
    "diagnosis_date", "patient_last_checkup_date", "treatment_date"] + ["treatment_%i_date" % i for i in range(1, 5)]
clinical = clinical[cols + traits + chrom_abrs + mutations].drop_duplicates()


# Data cleanup
# remove left censors (patients without birth or diagnosis date)
clinical = clinical[~clinical["patient_birth_date"].isnull()]
clinical = clinical[~clinical["diagnosis_date"].isnull()]


# Get duration of patient life or up to time of censorship
# but first, let"s replace missing "patient_last_checkup_date" with last date of treatment, death or collection
cols = ["patient_death_date", "sample_collection_date", "diagnosis_date"] + ["treatment_%i_date" % i for i in range(1, 5)]
subs = clinical[cols].apply(lambda x: max(x.dropna()), axis=1).groupby(clinical.index).max()
no_checkup = clinical[clinical["patient_last_checkup_date"].isnull()]

for i in no_checkup.index:
    if subs[i] is not pd.NaT:
        clinical.loc[i, "patient_last_checkup_date"] = subs[i]

clinical["duration"] = clinical.apply(life_duration, axis=1)


#


# Survival analysis
# For time since diagnosis
features = traits + chrom_abrs + mutations

# Survival of each clinical feature
fig, axis = plt.subplots(3, 7, figsize=(50, 20))
axis = axis.flatten()
time = "duration"
for i, feature in enumerate(features):
    survival_plot(clinical, KaplanMeierFitter(), "survival", feature, time, axis=axis[i])
fig.savefig(os.path.join(plots_dir, "all_traits.survival.svg"), bbox_inches="tight")


# Hazard of each clinical feature
fig, axis = plt.subplots(3, 7, figsize=(50, 20))
axis = axis.flatten()
time = "duration"
for i, feature in enumerate(features):
    survival_plot(clinical, NelsonAalenFitter(nelson_aalen_smoothing=False), "hazard", feature, time, axis=axis[i])
fig.savefig(os.path.join(plots_dir, "all_traits.hazard.svg"), bbox_inches="tight")


#


# Survival Regression
# exclude some traits
# due to being either rare or heavily assymetric between cohorts
features = traits + chrom_abrs + mutations
[features.pop(features.index(i)) for i in [
    "patient_gender",
    "diagnosis_disease",
    "diagnosis_stage_binet",
    "diagnosis_stage_rai",
    "ZAP70_monoallelic_methylation",
    "TP53", "del17", "SF3B1", "ATM", "NOTCH1", "BIRC3", "BCL2", "MYD88", "CHD2", "NFKIE"]]
combinations = [" + ".join(j) for i in range(1, len(features) + 1) for j in itertools.combinations(features, i)]

# Model selection
clinical2 = clinical.reset_index(drop=True)
fitters = ["aalen", "cox"]
for fitter_name in fitters:
    print(fitter_name)
    # Test all models
    performances = Parallel(n_jobs=-1, verbose=5)(delayed(test_model)(model, fitter_name, clinical2) for model in combinations)
    print("Done with cross-validation for model selection")

    # bit of post processing
    scores = pd.DataFrame(dict(zip(
        [re.sub(r" \+ ", "-", c) for c in combinations],
        sum(performances, []))
    )).T
    scores = scores.apply(lambda x: pd.Series(x[1]), axis=1)
    scores['penalty'] = scores[0]

    # mean  performance across folds
    scores['mean'] = scores[range(10)].apply(np.mean, axis=1)
    # model names
    scores = scores.reset_index()
    scores = scores.rename(columns={"index": "name"})
    # number of traits
    scores['traits'] = scores["name"].apply(lambda x: len(x.split("-")))

    # model concordance distribution
    fig, axis = plt.subplots(3, sharex=False, sharey=False, figsize=(15, 8))
    sns.distplot(scores["mean"], rug=True, ax=axis[0])
    # number of elements in model
    axis[1].scatter(scores["traits"], scores["mean"])
    # penalty per model
    axis[2].scatter(scores["penalty"], scores["mean"])
    axis[0].set_xlabel("Mean model concordance")
    axis[0].set_ylabel("Density")
    axis[1].set_xlabel("Number of model terms")
    axis[1].set_ylabel("Mean concordance")
    axis[2].set_xlabel("Co-variate penalty")
    axis[2].set_ylabel("Mean concordance")
    fig.savefig(os.path.join(plots_dir, "%s_model_selection_performance.svg" % fitter_name), bbox_inches="tight")

    # Sort models by performance
    best_models = scores.sort("mean")
    index_order = np.argsort(scores["mean"]).tolist()[-1]

    # Get 5 best models,
    # plot survival and hazard predictions
    fig, axis = plt.subplots(2, 5, sharex=False, sharey=False, figsize=(40, 15))
    fig2, axis2 = plt.subplots(5, sharex=False, sharey=True, figsize=(20, 15))
    # save
    predictions = pd.DataFrame()
    for m in range(1, 6):
        # get best model
        model = re.sub("-", " + ", best_models.ix[index_order[-m]]["name"])
        penalty = best_models.ix[index_order[-m]]["penalty"]

        # regress with best model
        X = patsy.dmatrix(model + " -1", clinical, return_type="dataframe")
        X["T"] = [i.days / 30. for i in clinical["duration"].ix[X.index]]
        X["E"] = [True if i is not pd.NaT else False for i in clinical["patient_death_date"].ix[X.index]]

        # Fit
        if fitter_name == "aalen":
            fitter = AalenAdditiveFitter(coef_penalizer=penalty, fit_intercept=True)
        elif fitter_name == "cox":
            fitter = CoxPHFitter(penalizer=penalty * 0.1)
        fitter.fit(X, "T", event_col="E")

        # Predict hazard, survival
        survival = fitter.predict_survival_function(X)
        survival["type"] = "survival"
        hazard = fitter.predict_cumulative_hazard(X)
        hazard["type"] = "hazard"

        # ...save
        for df in [survival, hazard]:
            df["model"] = model
            df["fitter"] = fitter_name
        predictions = predictions.append(survival)

        # ...and plot
        survival.plot(ax=axis[0][m - 1], legend=False)
        hazard.plot(ax=axis[1][m - 1], legend=False)
        axis[0][m - 1].set_title(model)
        axis[1][2].set_xlabel("Time since diagnosis")
        axis[0][0].set_ylabel("Survival")
        axis[1][0].set_ylabel("Hazard")
        # [i[m - 1].legend_.remove() for i in axis]

        # Hazard proportion of each trait
        if fitter_name == "cox":
            haz = fitter.hazards_.T.reset_index()
            haz = haz.sort('coef')
            sns.barplot("index", "coef", data=haz, ax=axis2[m])
            axis2[m].set_xlabel("Trait")
            axis2[m].set_ylabel("Hazard")

    #  save figures
    fig.savefig(os.path.join(plots_dir, "%s_model_best_5models_predictions_all_patients.svg" % fitter_name), bbox_inches="tight")
    fig2.savefig(os.path.join(plots_dir, "%s_model_hazard_per_trait.svg" % fitter_name), bbox_inches="tight")

    # save predictions
    predictions.columns = clinical.ix[predictions.columns]['patient_id']
    predictions = predictions.reset_index()
    predictions = predictions.rename(columns={"index": "time"})
    predictions.to_csv(os.path.join("data", "survival_hazard_predictions.csv"), index=False)
