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


def time_to_treatment(x, kind="diagnosis_date"):
    if x["date_first_treatment"] is pd.NaT:
        return x["patient_last_checkup_date"] - x[kind]
    else:
        return x["date_first_treatment"] - x[kind]


def survival_plot(clinical, fitter, fitter_name, feature, time="life_duration", event="patient_death_date", axis=None):
    """
    Plot survival/hazard of all patients regardless of trait and dependent of trait.
    """
    # duration of life
    T = [i.days / 30. for i in clinical[time]]
    # events:
    # True for observed event (death);
    # else False (this includes death not observed; death by other causes)
    C = [True if i is not pd.NaT else False for i in clinical[event]]

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
        axis.set_xlabel("time (months)")
    except:
        axis.set_xlabel("time (months)\n%s" % "\n".join(p_values))

    axis.set_title("%s" % feature)
    axis.set_ylabel(fitter_name)
    sns.despine()
    if save:
        fig.savefig(os.path.join(plots_dir, "%s_%s.svg" % (feature, fitter_name)), bbox_inches="tight")


# plots output dir
plots_dir = os.path.join("results", "plots", "survival")
# Get clinical info
clinical = pd.read_csv(
    os.path.join("metadata", "annotation.csv"),
    parse_dates=[
        "patient_birth_date", "patient_death_date", "date_first_treatment",
        "diagnosis_date", "patient_last_checkup_date"],
    dayfirst=True)

# traits to investigate (used at several points)
traits = [
    "patient_gender",
    "diagnosis_disease",
    "diagnosis_stage_binet",
    "ighv_mutation_status",
]

# Get one record per patient
cols = [
    "patient_id", "patient_birth_date", "patient_death_date", "date_first_treatment",
    "diagnosis_date", "patient_last_checkup_date"]
clinical = clinical[cols + traits].drop_duplicates()


# Data cleanup
# remove left censors (patients without birth or diagnosis date)
clinical = clinical[~clinical["patient_birth_date"].isnull()]
clinical = clinical[~clinical["diagnosis_date"].isnull()]


# Get duration of patient life or up to time of censorship
# but first, let"s replace missing "patient_last_checkup_date" with last date of treatment, death or collection
cols = ["patient_death_date", "diagnosis_date"]
subs = clinical[cols].apply(lambda x: max(x.dropna()), axis=1).groupby(clinical.index).max()
no_checkup = clinical[clinical["patient_last_checkup_date"].isnull()]

for i in no_checkup.index:
    if subs[i] is not pd.NaT:
        clinical.loc[i, "patient_last_checkup_date"] = subs[i]

# Get times
clinical["life_duration"] = clinical.apply(life_duration, axis=1)
clinical["time_to_treatment"] = clinical.apply(time_to_treatment, axis=1)


# Survival analysis
# For time since diagnosis
features = traits

for plot in [
    ("life_duration", "patient_death_date"),
        ("time_to_treatment", "date_first_treatment")]:
    time, event = plot
    # Survival of each clinical feature
    fig, axis = plt.subplots(1, 4, figsize=(20, 4))
    axis = axis.flatten()
    for i, feature in enumerate(features):
        survival_plot(clinical, KaplanMeierFitter(), "survival", feature, time=time, event=event, axis=axis[i])
    fig.savefig(os.path.join(plots_dir, "all_traits.%s.survival.svg" % time), bbox_inches="tight")

    # Hazard of each clinical feature
    fig, axis = plt.subplots(1, 4, figsize=(20, 4))
    axis = axis.flatten()
    for i, feature in enumerate(features):
        survival_plot(clinical, NelsonAalenFitter(nelson_aalen_smoothing=False), "hazard", feature, time=time, event=event, axis=axis[i])
    fig.savefig(os.path.join(plots_dir, "all_traits.%s.hazard.svg" % time), bbox_inches="tight")


# Take discovered sample groups
# read in
df = pd.read_csv(os.path.join(
    "data_submission", "cll_peaks.%s_significant.classification.random_forest.loocv.pca.sample_assignments.csv" % "IGHV"))
df["patient_id"] = df['sample'].apply(lambda x: x.split("_")[2]).astype(np.int64)

# get patients which all of their samples agree in the cluster
from collections import Counter
counts = pd.DataFrame(Counter(df[["patient_id", "cluster"]].drop_duplicates()['patient_id']).items())
patients = counts[counts[1] == 1][0]

# subset clinical with those patients
df2 = df[df["patient_id"].isin(patients)]
# merge with rest of clinical info
clinical_subset = pd.merge(
    clinical, df2[["patient_id", "cluster"]].drop_duplicates(), on=["patient_id"]).drop_duplicates()

feature = "cluster"

for plot in [
    ("life_duration", "patient_death_date"),
        ("time_to_treatment", "date_first_treatment")]:
    time, event = plot
    # Survival per ighv cluster
    fig, axis = plt.subplots(1, figsize=(8, 6))
    survival_plot(clinical_subset, KaplanMeierFitter(), "survival", feature, time=time, event=event, axis=axis)
    fig.savefig(os.path.join(plots_dir, "ighv_cluster.%s.survival.svg" % time), bbox_inches="tight")

    # Hazard per ighv cluster
    fig, axis = plt.subplots(1, figsize=(8, 6))
    survival_plot(clinical_subset, NelsonAalenFitter(nelson_aalen_smoothing=False), "hazard", feature, time=time, event=event, axis=axis)
    fig.savefig(os.path.join(plots_dir, "ighv_cluster.%s.hazard.svg" % time), bbox_inches="tight")
