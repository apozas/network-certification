# Code for
#
# Guarantees on the structure of experimental quantum networks
# npj Quantum Inf. 10, 117 (2024)
# arXiv:2403.02376
#
# Authors: Alejandro Pozas-Kerstjens
#
# Requires: inflation for setting up and solving the problems
#           numpy for array operations
#           qutip for quantum operations
#           sympy for symbolic operations
#           tqdm  for progress bars
#           itertools, numbers
#
# Last modified: Sep, 2023

import numpy as np
import qutip as qt

from inflation import InflationProblem, InflationSDP, max_within_feasible
from itertools import product
from numbers import Real
from sympy import Symbol
from tqdm import tqdm
from utils import export_inequality, prob_noin, rho, rho_ghz_list

meas = [[0.5*(qt.qeye(2)+qt.sigmax()),0.5*(qt.qeye(2)-qt.sigmax())],
        [0.5*(qt.qeye(2)+qt.sigmaz()),0.5*(qt.qeye(2)-qt.sigmaz())]]

dag = {"h1": ["A", "B", "C"],
       "h2": ["B", "C", "D", "E"],
       "h3": ["D", "E", "F"]}

vis = Symbol("v")

ineq_path = "GHZInequalities"
measurements_list = ["".join(m)
                     for m in product(["X", "Z"], repeat=6)]

###############################################################################
# Single-input inequalities
###############################################################################
InfProb = InflationProblem(dag=dag,
                           outcomes_per_party=[2, 2, 2, 2, 2, 2],
                           inflation_level_per_source=2,
                           verbose=0
                           )

InfSDP = InflationSDP(InfProb)

Local1Len2 = InfSDP.build_columns("local1", max_monomial_length=2)
info = "Local1Len2"

InfSDP.generate_relaxation(Local1Len2)
try:
    with open(f"ghzvisibilities_INF2{info}.txt", "r") as file:
        visibilities = file.read()
except FileNotFoundError:
    visibilities = ""
for measurements in tqdm(measurements_list, desc=f"Inequalities for {info}"):
    if measurements in visibilities:
        pass
    else:
        InfSDP.set_distribution(prob_noin(vis, "ghz", measurements))
        if any([not isinstance(val, Real)
                for val in InfSDP.known_moments.values()]):
            vcrit = max_within_feasible(InfSDP, InfSDP.known_moments, "dual")
            if abs(vcrit - 1) > 1e-3:
                InfSDP.reset("all")
                v = np.ceil(vcrit * 1000) / 1000
                InfSDP.set_distribution(prob_noin(v, "ghz", measurements))
                InfSDP.solve(feas_as_optim=True)
                if InfSDP.solution_object["status"] == "feasible":
                    cert = InfSDP.certificate_as_probs()
                    export_inequality(cert,
                                f"{ineq_path}/{measurements}_INF2{info}.csv")
                    visibilities += \
                        f"Critical visibility for {measurements} is {vcrit}\n"
                else:
                    print(f"Problem for {measurements}, failed with status "
                          + InfSDP.solution_object["status"])
            else:
                visibilities += f"Critical visibility for {measurements} is 1\n"
        else:
            visibilities += f"{measurements} produces a constant distribution\n"

    with open(f"ghzvisibilities_INF2{info}.txt", "w") as fileexport:
        fileexport.write(visibilities)

###############################################################################
# Binary-input inequality
###############################################################################
def prob_2in(vis):
    prob_array = np.zeros((2,2,2,2,2,2,2,2,2,2,2,2))
    if isinstance(vis, Real):
        state = rho("ghz", vis)
        for a,b,c,d,e,f,x,y,z,t,u,v in np.ndindex((2,2,2,2,2,2,2,2,2,2,2,2)):
            prob_array[a,b,c,d,e,f,x,y,z,t,u,v] \
                = qt.expect(state, qt.tensor(meas[x][a],
                                             meas[y][b],
                                             meas[z][c],
                                             meas[t][d],
                                             meas[u][e],
                                             meas[v][f]))
    else:
        states = rho_ghz_list()
        prob_array = np.asarray(prob_array, dtype=object)
        for a,b,c,d,e,f,x,y,z,t,u,v in np.ndindex((2,2,2,2,2,2,2,2,2,2,2,2)):
            operator = qt.tensor(meas[x][a], meas[y][b],
                                 meas[z][c], meas[t][d],
                                 meas[u][e], meas[v][f])
            prob_array[a,b,c,d,e,f,x,y,z,t,u,v] = (
                vis**3 * qt.expect(states[0], operator)
                + vis**2 * qt.expect(states[1], operator)
                + vis * qt.expect(states[2], operator)
                + qt.expect(states[3], operator))
    return prob_array

InfProb = InflationProblem(dag=dag,
                       outcomes_per_party=[2, 2, 2, 2, 2, 2],
                       settings_per_party=[2, 2, 2, 2, 2, 2],
                       inflation_level_per_source=2,
                       verbose=0
                       )

InfSDP = InflationSDP(InfProb)
info = "NPA1"
InfSDP.generate_relaxation(info)
InfSDP.set_distribution(prob_2in(vis))
if any([not isinstance(val, Real)
        for val in InfSDP.known_moments.values()]):
    vcrit = max_within_feasible(InfSDP, InfSDP.known_moments, "dual")
    InfSDP.reset("all")
    InfSDP.set_distribution(prob_2in(min(vcrit, 1.)))
    InfSDP.solve(feas_as_optim=True)
    if InfSDP.solution_object["status"] == "feasible":
        cert = InfSDP.certificate_as_probs()
        export_inequality(cert, f"{ineq_path}/twoin_INF2{info}.csv")
    print(f"Critical visibility for INF2{info} is {vcrit}")
