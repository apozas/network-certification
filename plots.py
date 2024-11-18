# Code for
#
# Guarantees on the structure of experimental quantum networks
# npj Quantum Inf. 10, 117 (2024)
# arXiv:2403.02376
#
# Authors: Alejandro Pozas-Kerstjens
#
# Requires: matplotlib, seaborn for plotting
#           numpy     for array operations
#           pandas    for dataframes
#           symengine for symbolic operations
#           itertools, pickle, os
# 
# Last modified: Sep, 2023

import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import pandas as pd
import pickle

from itertools import product
from utils import import_ineq, eval_ineq, prob_noin, prob_2in

ghz_results   = "inequality_evaluations_ghz.pkl"
ghz_samples   = "inequality_samples_ghz_2input.pkl"
ghz_ineq_path = "GHZInequalities"
trident_results   = "inequality_evaluations_trident.pkl"
trident_ineq_path = "TridentInequalities"

# Code from https://towardsdatascience.com/better-heatmaps-and-correlation-matrix-plots-in-python-41445d0f2bec
n_colors = 512
palette = sns.diverging_palette(20, 220, n=2*n_colors)
palette = palette[:len(palette)//2]
palette = list(reversed(sns.dark_palette((179/255, 49/255, 48/255),
                                         n_colors=n_colors)))
color_min, color_max = [-1, 0] # Range that will be mapped to the palette

def value_to_color(val, max_size):
    if (val == 0):
        return (234 / 255, 234 / 255, 242 / 255)
    elif val == 1:
        return (246 / 255, 210 / 255, 189 / 255)
    elif abs(val) < 1e-6:
        return (0.5, 1., 1.)
    elif val > 1e-6:
        return (0.24715576253545807, 0.49918708160096675, 0.5765599057376697)
    else:
        ind = int(val / max_size * (len(palette) - 1)) - 1
        return palette[ind]

def value_to_size(val, max_size):
    if abs(val) < 1e-6:
        return 0
    elif val > -1e-6:
        return 1
    else:
        return abs(val) / max_size

def heatmap(x, y, z, title=None, smallest_value=None, figsize=(15,15)):
    _, ax = plt.subplots(figsize=figsize)
    
    # Mapping from column names to integer coordinates
    x_labels = [v for v in sorted(x.unique())]
    y_labels = [v for v in sorted(y.unique())]
    x_to_num = {p[1]: p[0] for p in enumerate(x_labels)} 
    y_to_num = {p[1]: p[0] for p in enumerate(reversed(y_labels))} 
    
    
    plot_grid = plt.GridSpec(1, 35, hspace=0.2, wspace=0.4) # Setup a 1x15 grid
    # Use the leftmost 14 columns of the grid for the main plot
    ax = plt.subplot(plot_grid[:,:-1])

    size_scale = 100
    if not smallest_value:
        smallest_value = z[z < 0].min()
        max_size = abs(smallest_value)
    else:
        max_size = abs(smallest_value)
    ax.scatter(
        x=x.map(x_to_num),
        y=y.map(y_to_num),
        s=z.apply(value_to_size, args=[max_size]) * size_scale,
        c=z.apply(value_to_color, args=[max_size]),
        marker="s"
    )
    
    # Show column labels on the axes
    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment="right")
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels)
    ax.grid(False, "major");
    ax.grid(True, "minor");
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5]) 
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    ax.set_xlabel("Data measurement basis", fontsize=20)
    ax.set_ylabel("Witness measurement basis", fontsize=20)
    if title is not None:
        ax.set_title(title, fontsize=16)
    
    # Add color legend on the right side of the plot
    # Use the rightmost column of the plot
    ax = plt.subplot(plot_grid[:,-1], frameon=False)

    col_x = [0]*len(palette) # Fixed x coordinate for the bars
    bar_y = np.linspace(color_min, color_max, n_colors)

    bar_height = bar_y[1] - bar_y[0]
    ax.barh(
        y=bar_y,
        width=[5]*len(palette), # Make bars 5 units wide
        left=col_x, # Make bars start at 0
        height=bar_height,
        color=palette,
        linewidth=0
    )
    ax.set_xlim(1, 2) # Bars go from 0 to 5, so crop somewhere in the middle
    ax.grid(False) # Hide grid
    ax.set_facecolor("white") # Make background white
    ax.set_xticks([]) # Remove horizontal ticks
    ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3),
                  labels=[np.round(smallest_value, 4),
                          np.round(smallest_value/2, 4),
                          "0.0000"],
                  fontsize=20) # Show vertical ticks for min, middle and max
    ax.yaxis.tick_right() # Show vertical ticks on the right

###############################################################################
# GHZ: Two-input distributions (Figure 3)
###############################################################################
sns.set(style="ticks", font_scale=1.8)
plt.rcParams.update({"font.size": 18,
                     "text.usetex": True,
                     "font.family": "serif"})
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")

colors = sns.color_palette()

with open(ghz_results, "rb") as file:
    data = pickle.load(file)
    
with open(ghz_samples, "rb") as file:
    samples = pickle.load(file)

ineq = [key for key in data.keys() if key.startswith("twoin")][0]
twoin_data = data[ineq]
percentages = np.asarray(list(twoin_data.keys()))
means = np.array([twoin_data[perc][0] for perc in percentages]).astype(float)
stds  = np.array([twoin_data[perc][1] for perc in percentages])
plt.plot(percentages*290464, means,
        color=colors[0]
        )
plt.fill_between(percentages*290464,
                means-5*stds,
                means+5*stds, alpha=0.3,
                color=colors[0]
                )
for perc in percentages:
    plt.scatter([perc*290464]*100, samples[ineq][perc],
                color=colors[1], alpha=0.3, s=10)

plt.scatter(percentages*290464, means,
            color=colors[0])
plt.hlines(0, 0.01*290464, 290464,
        linestyles="dashed", linewidth=0.5, color="black")
plt.xscale("log")
plt.xlabel(r"Number of events")
plt.ylabel(r"$\mathcal{W}_{1}$");
plt.xticks([3e3, 1e4, 1e5], ["$3\!\!\\times\!\!10^3$", "$10^4$", "$10^5$"])
plt.savefig(f"Plots/GHZ_{ineq}.pdf", bbox_inches="tight")
plt.clf();

###############################################################################
# GHZ: Single-input distributions (Figures 4, S3, S4, S5)
###############################################################################
sns.set(font_scale=1.)
plt.rcParams.update({"font.size": 12})
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.rc("font", family="serif")

# # Figures 4a, S3
allprobs = {"".join(measurements): prob_noin(1, "ghz", measurements)
            for measurements in product(["X", "Z"], repeat=6)}
short = ["ZZXXZX", "ZZXXZZ", "ZZXZXX", "ZZXZXZ", "ZZXZZX", "ZZXZZZ", "ZZZXXX",
         "ZZZXXZ", "ZZZXZX", "ZZZXZZ", "ZZZZXX", "ZZZZXZ", "ZZZZZX", "ZZZZZZ"]
inequalities = [path for path in os.listdir(ghz_ineq_path)
                if path[:6] != "twoin_"]
ineq_bases = [ineq[:6] for ineq in inequalities]

minvals = []
for fig_bases in [allprobs.keys(), short]:
    theo_dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value"])
    for ineq_bases in fig_bases:
        ineq_path = [path for path in inequalities
                    if path.startswith(ineq_bases)]
        for eval_bases in fig_bases:
            if len(ineq_path) == 0:
                value = 0
            else:
                value = eval_ineq(
                    import_ineq(f"{ghz_ineq_path}/{ineq_path[0]}"),
                                allprobs[eval_bases])
            theo_dataframe = pd.concat([theo_dataframe,
                                        pd.DataFrame({"Ineq": ineq_bases,
                                                      "Data": eval_bases,
                                                      "Value": value},
                                                     index=[0])],
                                       ignore_index=True)
    minvals.append(theo_dataframe.Value.min())

    if len(fig_bases) == 64:
        extra = "all"
        figsize = (15, 15)
    else:
        extra = "short"
        figsize = (4, 4)

    heatmap(
        x=theo_dataframe["Data"],
        y=theo_dataframe["Ineq"],
        z=theo_dataframe["Value"],
        figsize=figsize
    )

    plt.savefig(f"Plots/theoretical_ghz_{extra}.pdf", bbox_inches="tight")
    plt.clf();


# Figures 4b, S4
with open(ghz_results, "rb") as file:
    data = pickle.load(file)

for pos, fig_bases in enumerate([allprobs.keys(), short]):
    dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value"])
    for ineq_bases in fig_bases:
        ineq_path = [path for path in inequalities
                    if path.startswith(ineq_bases)]
        for eval_bases in fig_bases:
            if len(ineq_path) == 0:
                value = 0
            else:
                value = float(data[ineq_bases][eval_bases][1.][0])
            dataframe = pd.concat([dataframe,
                                   pd.DataFrame({"Ineq": ineq_bases,
                                                 "Data": eval_bases,
                                                 "Value": value},
                                                index=[0])],
                                  ignore_index=True)

    extra   = "all" if len(fig_bases) == 64 else "short"
    figsize = (15, 15) if len(fig_bases) == 64 else (4, 4)
    heatmap(
        x=dataframe["Data"],
        y=dataframe["Ineq"],
        z=dataframe["Value"],
        smallest_value=minvals[pos],
        figsize=figsize
    )
    plt.savefig(f"Plots/experimental_ghz_{extra}.pdf", bbox_inches="tight")
    plt.clf();

# Figure S5
sns.set(style="ticks", font_scale=3.0)
plt.rcParams.update({"font.size": 30})
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.rc("font", family="serif")

ineqs = ["ZZXXXZ", "ZZXXZX", "ZZXXZZ", "ZZXZXX", "ZZXZXZ", "ZZXZZX", "ZZXZZZ"]
datas = ["ZZZXZX", "ZZZXZZ", "ZZZZXX", "ZZZZXZ", "ZZZZZX", "ZZZZZZ"]
counts = [4156, 4329, 4245, 4404, 4281, 4220]

markers = "ox*dsph"
colors = sns.color_palette()

for bases, c in zip(datas, counts):
    plt.figure()
    for idx, ineq in enumerate(ineqs):
        basis_data = data[ineq][bases]
        percentages = np.asarray(list(basis_data.keys()))
        means = np.array([basis_data[perc][0] for perc in percentages]
                        ).astype(float)
        stds = np.array([basis_data[perc][1] for perc in percentages])
        plt.plot(percentages*c, means,
                 color=colors[idx],
                 marker=markers[idx])
        plt.scatter(percentages*c, means,
                    color=colors[idx], marker=markers[idx])
        plt.fill_between(percentages*c,
                         means-5*stds,
                         means+5*stds, alpha=0.3, color=colors[idx])
        plt.hlines(0, 0.01*c, c,
                   linestyles="dashed", linewidth=0.5, color="black")
    plt.xlabel("Number of events")
    plt.title(f"{bases} data")
    plt.savefig(f"Plots/GHZ_{bases}.pdf", bbox_inches="tight")
    plt.clf();

###############################################################################
# Trident: Single-input distributions (Figure S7)
###############################################################################
sns.set(font_scale=1.)
plt.rcParams.update({"font.size": 12})
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
plt.rc("font", family="serif")

# Theoretical plot (figure S7a)
ineqs_trident_l13 = ["YYYYYY", "YYZXZX", "YYXYYY", "YYZZYY", "ZXZZYY", "YYXZYY",
                     "YYXZZX", "YYZXYY", "ZXZYYY", "YYXXYY", "ZXZZZX", "ZXYZYY",
                     "YYZYZX", "ZXXXZX", "YYYZZX", "ZXYYZX", "YYZYYY", "YYYZYY",
                     "ZXZXYY", "YYZZZX", "YYYXYY", "ZXXZYY"]
datas_trident_l13 = ["YYXXYY", "YYXYYY", "YYXZYY", "YYYXYY", "YYYYYY", "YYYZYY",
                     "YYZXYY", "YYZYYY", "YYZZYY", "ZYYZYY"]
allprobs_l13 = {measurements: prob_noin(1, "trident", measurements)
                for measurements in datas_trident_l13}

inequalities = [path for path in os.listdir(trident_ineq_path)
                if path[:6] != "twoin_"]
theo_dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value"])
for ineq_bases, eval_bases in product(ineqs_trident_l13, datas_trident_l13):
    ineq_path = [path for path in inequalities if path.startswith(ineq_bases)]
    theo_dataframe = pd.concat(
        [theo_dataframe,
         pd.DataFrame({"Ineq": ineq_bases,
                       "Data": eval_bases,
                       "Value": eval_ineq(import_ineq(
                           f"{trident_ineq_path}/{ineq_path[0]}"),
                           allprobs_l13[eval_bases])},
                       index=[0])],
        ignore_index=True)

heatmap(
    x=theo_dataframe["Data"],
    y=theo_dataframe["Ineq"],
    z=theo_dataframe["Value"],
    figsize=(2.5,5)
)
plt.savefig("Plots/theoretical_trident.pdf", bbox_inches="tight")

# Experimental plot (figure S7b)
with open(trident_results, "rb") as file:
    data = pickle.load(file)
dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value"])
for ii, jj in product(ineqs_trident_l13, datas_trident_l13):
    dataframe = pd.concat([dataframe,
                           pd.DataFrame(
                               {"Ineq": [ii],
                                "Data": [jj],
                                "Value": [float(data[ii][jj][1.][0])]})],
                          ignore_index=True)

heatmap(
    x=dataframe["Data"],
    y=dataframe["Ineq"],
    z=dataframe["Value"],
    smallest_value=theo_dataframe.Value.min(),
    figsize=(2.5,5)
)
plt.savefig("Plots/experimental_trident.pdf", bbox_inches="tight")

###############################################################################
# Trident: Two-input distributions (Table I)
###############################################################################
# Theoretical distributions
allprobs = {measurements: prob_2in(1, "trident", measurements)
            for measurements in ["XY", "XZ", "YX", "YZ", "ZX", "ZY"]}

theo_dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value"])
inequalities = [path for path in os.listdir(trident_ineq_path)
                if path.startswith("twoin")]
for ineq_bases in ["XY", "XZ", "YZ"]:
    ineq_path = [path for path in inequalities if ineq_bases in path][0]
    for bases in allprobs.keys():
        value = eval_ineq(import_ineq(f"{trident_ineq_path}/{ineq_path}"),
                          allprobs[bases])
        theo_dataframe = pd.concat([theo_dataframe,
                                    pd.DataFrame({"Ineq":  ineq_bases,
                                                  "Data":  bases,
                                                  "Value": value},
                                                 index=[0])],
                                   ignore_index=True)

print("Trident two-input distributions: theory")
print(theo_dataframe)

# Experimental distributions (Table I)
available_bases = list(data.keys())
exp_dataframe = pd.DataFrame(columns=["Ineq", "Data", "Value", "Std"])
for ineq in data.keys():
    if ineq.startswith("twoin"):
        for bases, val in data[ineq].items():
            exp_dataframe = pd.concat([exp_dataframe,
                                       pd.DataFrame({"Ineq": ineq.split('_')[1],
                                                     "Data": bases,
                                                     "Value": val[0],
                                                     "Std":   val[1]},
                                                    index=[0])],
                                      ignore_index=True)
print("Trident two-input distributions: experiment")
print(exp_dataframe)
