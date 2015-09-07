
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

sns.set_style("whitegrid")
sns.set_context("paper")


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


df = pd.read_csv("mergePeaks.raw.txt", sep="\t")
df = df.reset_index(drop=True)
df = np.array(df)

df2 = pd.read_csv("mergePeaks.raw.qnorm.txt", sep="\t")
df2 = df2.reset_index(drop=True)
df2 = df2.drop(['Transcript_id', 'GeneSymbol'], axis=1)
df2 = np.array(df2)

df_qnorm = normalize_quantiles_r(df)

# Compare raw counts vs qnormalized data
fig, axis = plt.subplots(3)
[sns.distplot(np.log2(1 + df[:, i]), ax=axis[0], hist=False) for i in range(df.shape[1])]
[sns.distplot(np.log2(1 + df2[:, i]), ax=axis[1], hist=False) for i in range(df2.shape[1])]
[sns.distplot(np.log2(1 + df_qnorm[:, i]), ax=axis[2], hist=False) for i in range(df_qnorm.shape[1])]
axis[0].set_title("raw")
axis[1].set_title("their norm")
axis[2].set_title("our norm")

fig.savefig(os.path.join("cell_systems.coverage_vs_coverage_qnorm.pdf"))

# four particular site
fig, axis = plt.subplots(2)
[axis[0].scatter(np.log2(1 + df[100, i]), 1) for i in range(df.shape[1])]
[axis[1].scatter(np.log2(1 + df2[100, i]), 1) for i in range(df2.shape[1])]

[axis[0].scatter(np.log2(1 + df[120, i]), 2) for i in range(df.shape[1])]
[axis[1].scatter(np.log2(1 + df2[120, i]), 2) for i in range(df2.shape[1])]

[axis[0].scatter(np.log2(1 + df[140, i]), 3) for i in range(df.shape[1])]
[axis[1].scatter(np.log2(1 + df2[140, i]), 3) for i in range(df2.shape[1])]

[axis[0].scatter(np.log2(1 + df[160, i]), 4) for i in range(df.shape[1])]
[axis[1].scatter(np.log2(1 + df2[160, i]), 4) for i in range(df2.shape[1])]

fig.savefig(os.path.join("cell_systems.coverage_vs_coverage_qnorm.single_sites.pdf"), bbox_inches="tight")
