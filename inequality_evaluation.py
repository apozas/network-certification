# Code for
#
# Guarantees on the network structure of experimental quantum networks
#
# arXiv:2310.xxxxx
#
# Authors: Alejandro Pozas-Kerstjens
#
# Requires: numpy     for array operations
#           pandas    for dataframes
#           symengine for symbolic operations
#           tqdm      for progress bars
#           csv, itertools, pickle, os
#
# Last modified: Sep, 2023

import numpy as np
import pandas as pd
import pickle
import os

from itertools import product
from symengine import symbols
from tqdm import tqdm
from utils import expand_ineq, expand_ineq_2in, error_expr, \
                  get_means_stds, get_distr_from_data, read_ineq

trident_data_path = "TRIDENT_DATA_PATH"
trident_ineq_path = "TRIDENT_INEQ_PATH"
ghz_data_path     = "GHZ_DATA_PATH"
ghz_ineq_path     = "GHZ_INEQ_PATH"

###############################################################################
# Empirical distributions
###############################################################################
percentages = np.logspace(-2, 0, 30)

ghz     = [ghz_data_path, ghz_ineq_path, percentages, 100]
trident = [trident_data_path, trident_ineq_path, [1.0], 1]

for data_path, ineq_path, percs, runs in [ghz, trident]:
    exp_probs = {}
    prob_samples = {}
    n_counts = 0
    for data in tqdm(os.listdir(data_path),
                     desc=f"Building distributions for {ineq_path[:-12]}"):
        data_name = data[:6]
        exp_probs[data_name] = {}
        prob_samples[data_name] = {}
        for perc in percs:
            allprobs = []
            for _ in range(runs):
                prob, counts = get_distr_from_data(f"{data_path}/{data}",
                                                   ineq_path[:-12],
                                                   perc)
                allprobs.append(prob)
            prob_mean = np.mean(allprobs, axis=0)
            prob_std  = np.std(allprobs,  axis=0)
            exp_probs[data_name] = {**exp_probs[data_name],
                                    **{perc: [prob_mean, prob_std]}}
            prob_samples[data_name] = {**prob_samples[data_name],
                                       **{perc: allprobs}}
        n_counts += counts

    print(f"Total counts for {ineq_path[:-12]}: {n_counts}")

    with open(f"prob_estimations_{ineq_path[:-12].lower()}.pkl", "wb") as file:
        pickle.dump(exp_probs, file)

    if ineq_path.startswith("GHZ"):
        with open(f"prob_samples_ghz.pkl", "wb") as file:
            pickle.dump(prob_samples, file)

###############################################################################
# Inequalities for GHZ data
###############################################################################
names = ["".join(letters) for letters in product(["T","R"], repeat=6)]

with open(f"prob_estimations_ghz.pkl", "rb") as file:
    exp_probs = pickle.load(file)
available_bases = list(exp_probs.keys())

exp_ineqs = {}

# Single-input inequalities
inequalities = [path for path in os.listdir(ghz_ineq_path)
                if not path.startswith("twoin")]

for inequality in tqdm(inequalities, desc="GHZ, single input"):
    ineq_name = inequality[:6]
    ineq  = expand_ineq(read_ineq(f"{ghz_ineq_path}/{inequality}"))
    error = error_expr(ineq)
    exp_ineqs[ineq_name] = {}
    for data in exp_probs.keys():
        ineq_evaluations = {}
        for perc in exp_probs[data].keys():
            prob, std = exp_probs[data][perc]
            prob_subs = {symbols("p_{" + name + "}"): pr
                         for name, pr in zip(names, prob)}
            exp_ineq = ineq.subs(prob_subs)
            error_subs = {symbols("\Delta_{" + name + "}"): pr
                          for name, pr in zip(names, std)}
            exp_error = np.sqrt(float(error.subs({**prob_subs, **error_subs})))
            ineq_evaluations = {**ineq_evaluations,
                                **{perc: [exp_ineq, exp_error]}}
        exp_ineqs[ineq_name][data] = ineq_evaluations

with open(f"inequality_evaluations_ghz.pkl", "wb") as file:
    pickle.dump(exp_ineqs, file)

# Two-input inequalities
with open(f"prob_samples_ghz.pkl", "rb") as file:
    prob_samples = pickle.load(file)
ineq_samples = {}

inequality_names = [file
                    for file in os.listdir(ghz_ineq_path)
                    if file.startswith("twoin")]

for inequality in inequality_names:
    # Average values
    ineq_name = inequality.split(".")[0]
    ineq = read_ineq(f"{ghz_ineq_path}/{inequality}", num_inputs=2)
    twoin_ineq  = expand_ineq_2in(ineq, available_bases, "ZX")
    twoin_error = error_expr(twoin_ineq)

    exp_2in_ineq = {}
    for perc in exp_probs[list(exp_probs.keys())[0]].keys():
        prob_subs = {}
        error_subs = {}
        for bases in exp_probs.keys():
            prob_subs = {**prob_subs,
                        **{symbols("p_{" + name + "}(" + bases + ")"): pr
                            for name, pr in zip(names,
                                                exp_probs[bases][perc][0])}}
            error_subs = {**error_subs,
                        **{symbols("\Delta_{" + name + "}(" + bases + ")"): e
                            for name, e in zip(names,
                                                exp_probs[bases][perc][1])}}
        exp_ineq  = twoin_ineq.subs(prob_subs)
        exp_error = np.sqrt(float(twoin_error.subs({**prob_subs,
                                                    **error_subs})))
        exp_2in_ineq[perc] = [exp_ineq, exp_error]

        exp_ineqs[f"{ineq_name}"] = exp_2in_ineq

    # Individual datapoints
    samples = {}
    for perc in tqdm(prob_samples[list(prob_samples.keys())[0]].keys(),
                     desc=f"Evaluating GHZ {ineq_name} in samples of data"):
        evals = []
        for ii in range(100):
            prob_subs = {}
            for bases in prob_samples.keys():
                prob_subs = {**prob_subs,
                             **{symbols("p_{" + name + "}(" + bases + ")"): pr
                                for name, pr in zip(names,
                                                prob_samples[bases][perc][ii])}}
            evals.append(twoin_ineq.subs(prob_subs))

        samples[perc] = evals
    ineq_samples[ineq_name] = samples

with open(f"inequality_evaluations_ghz.pkl", "wb") as file:
    pickle.dump(exp_ineqs, file)

with open(f"inequality_samples_ghz_2input.pkl", "wb") as file:
    pickle.dump(ineq_samples, file)

###############################################################################
# Inequalities for Trident data
###############################################################################
with open(f"prob_estimations_trident.pkl", "rb") as file:
    exp_probs = pickle.load(file)
available_bases = list(exp_probs.keys())

# Single-input inequalities
exp_ineqs = {}
inequalities = [path for path in os.listdir(trident_ineq_path)
                if not path.startswith("twoin")]

for inequality in tqdm(inequalities, desc="Trident, single input"):
    ineq_name = inequality[:6]
    ineq  = expand_ineq(read_ineq(f"{trident_ineq_path}/{inequality}"))
    error = error_expr(ineq)
    exp_ineqs[ineq_name] = {}
    for data in exp_probs.keys():
        ineq_evaluations = {}
        for perc in exp_probs[data].keys():
            prob, std = exp_probs[data][perc]
            prob_subs = {symbols("p_{" + name + "}"): pr
                         for name, pr in zip(names, prob)}
            exp_ineq = ineq.subs(prob_subs)
            ineq_evaluations = {**ineq_evaluations, **{perc: [exp_ineq, 0.]}}
        exp_ineqs[ineq_name][data] = ineq_evaluations

with open(f"inequality_evaluations_trident.pkl", "wb") as file:
    pickle.dump(exp_ineqs, file)

# Two-input inequalities
inequality_fnames = [file
                     for file in os.listdir(trident_ineq_path)
                     if file.startswith("twoin")]
for inequality in inequality_fnames:
    ineq = read_ineq(f"{trident_ineq_path}/{inequality}", num_inputs=2)
    err = error_expr(ineq)
    ineq_bases = inequality.split("_")[1]
    exp_ineqs[f"twoin_{ineq_bases}"] = {}
    for eval_bases in ["XY", "YX", "XZ", "ZX", "YZ", "ZY"]:
        means, stds = get_means_stds(ineq,
                                     available_bases,
                                     eval_bases,
                                     exp_probs)
        ineq_eval   = ineq.subs(means)
        err_eval    = np.sqrt(float(err.subs({**means, **stds})))
        exp_ineqs[f"twoin_{ineq_bases}"][eval_bases] = [ineq_eval, err_eval]

with open(f"inequality_evaluations_trident.pkl", "wb") as file:
    pickle.dump(exp_ineqs, file)
