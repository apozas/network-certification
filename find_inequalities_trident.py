# Code for
#
# Guarantees on the structure of experimental quantum networks
#
# arXiv:2403.02376
#
# Authors: Alejandro Pozas-Kerstjens
#
# Requires: inflation for setting up and solving the problems
#           numpy for array operations
#           qutip for quantum operations
#           sympy for symbolic operations
#           tqdm  for progress bars
#           numbers, os
# 
# Last modified: Mar, 2024

import os
import numpy as np
import qutip as qt

from inflation import InflationProblem, InflationSDP, max_within_feasible
from numbers import Real
from sympy import Symbol
from tqdm import tqdm
from utils import export_inequality, prob_noin, rho, rho_trident_list

dag = {"h1": ["A", "B", "C", "D"],
       "h4": ["C", "D", "E", "F"]}

vis = Symbol("v")

data_path = "DATA_PATH"
ineq_path = "TridentInequalities"
measurements_list = [file[:6] for file in os.listdir(data_path)]

###############################################################################
# Single-input inequalities
###############################################################################
InfProb = InflationProblem(dag=dag,
                           outcomes_per_party=[2, 2, 2, 2, 2, 2],
                           inflation_level_per_source=2,
                           verbose=0
                           )

InfSDP = InflationSDP(InfProb)

Local1Len3 = InfSDP.build_columns('local1', max_monomial_length=3)
info = "Local1Len3"

InfSDP.generate_relaxation(Local1Len3)
try:
    with open(f"tridentvisibilities_INF2{info}.txt", "r") as file:
        visibilities = file.read()
except FileNotFoundError:
    visibilities = ""
for measurements in tqdm(measurements_list, desc=f"Inequalities for {info}"):
    if measurements in visibilities:
        pass
    else:
        InfSDP.set_distribution(prob_noin(vis, "trident", measurements))
        if any([not isinstance(val, Real)
                for val in InfSDP.known_moments.values()]):
            vcrit = max_within_feasible(InfSDP, InfSDP.known_moments, "dual")
            if abs(vcrit - 1) > 1e-3:
                InfSDP.reset("all")
                v = np.ceil(vcrit * 1000) / 1000
                InfSDP.set_distribution(prob_noin(v, "trident", measurements))
                InfSDP.solve(feas_as_optim=True)
                if InfSDP.solution_object["status"] == "feasible":
                    cert = InfSDP.certificate_as_probs()
                    export_inequality(cert,
                                f'{ineq_path}/{measurements}_INF2{info}.csv')
                    visibilities += \
                        f"Critical visibility for {measurements} is {vcrit}\n"
                else:
                    print(f"Problem for {measurements} failed with status "
                          + InfSDP.solution_object["status"])
            else:
                visibilities += f"Critical visibility for {measurements} is 1\n"
        else:
            visibilities += f"{measurements} produces a constant distribution\n"

    with open(f"tridentvisibilities_INF2{info}.txt", "w") as fileexport:
        fileexport.write(visibilities)

###############################################################################
# Binary-input inequalities
###############################################################################
meas = [[0.5*(qt.qeye(2)+qt.sigmax()), 0.5*(qt.qeye(2)-qt.sigmax())],
        [0.5*(qt.qeye(2)+qt.sigmay()), 0.5*(qt.qeye(2)-qt.sigmay())],
        [0.5*(qt.qeye(2)+qt.sigmaz()), 0.5*(qt.qeye(2)-qt.sigmaz())]]

def prob_2in(vis, bases):
    prob_array = np.zeros((2,2,2,2,2,2,2,2,2,2,2,2))
    msmnts = [meas[bases[0]], meas[bases[1]]]
    if isinstance(vis, Real):
        state = rho("trident", vis)
        for a,b,c,d,e,f,x,y,z,t,u,v in np.ndindex((2,2,2,2,2,2,2,2,2,2,2,2)):
            prob_array[a,b,c,d,e,f,x,y,z,t,u,v] \
                = qt.expect(state, qt.tensor(msmnts[x][a],
                                             msmnts[y][b],
                                             msmnts[z][c],
                                             msmnts[t][d],
                                             msmnts[u][e],
                                             msmnts[v][f]))
    else:
        states = rho_trident_list()
        prob_array = np.asarray(prob_array, dtype=object)
        for a,b,c,d,e,f,x,y,z,t,u,v in np.ndindex((2,2,2,2,2,2,2,2,2,2,2,2)):
            operator = qt.tensor(msmnts[x][a], msmnts[y][b],
                                 msmnts[z][c], msmnts[t][d],
                                 msmnts[u][e], msmnts[v][f])
            prob_array[a,b,c,d,e,f,x,y,z,t,u,v] = (
                vis**4 * qt.expect(states[0], operator)
                + vis**3 * qt.expect(states[1], operator)
                + vis**2 * qt.expect(states[2], operator)
                + vis * qt.expect(states[3], operator)
                + qt.expect(states[4], operator))
    return prob_array

InfProb = InflationProblem(dag=dag,
                           outcomes_per_party=[2, 2, 2, 2, 2, 2],
                           settings_per_party=[2, 2, 2, 2, 2, 2],
                           inflation_level_per_source=2,
                           verbose=0
                           )

InfSDP = InflationSDP(InfProb)

Local1Len2 = InfSDP.build_columns('local1', max_monomial_length=2)
info = "Local1Len2"

InfSDP.generate_relaxation(Local1Len2)
for measurements in tqdm([[0, 1], [0, 2], [1, 2]],
                         desc="Getting 2-input inequalities for INF2" + info):
    if measurements == [0, 1]:
        input_names = 'XY'
    elif measurements == [0, 2]:
        input_names = 'XZ'
    elif measurements == [1, 2]:
        input_names = 'YZ'
    else:
        print("An error with the parsing of measurements has occurred")
        break
    InfSDP.reset("all")
    InfSDP.set_distribution(prob_2in(vis, measurements))
    if any([not isinstance(val, Real)
            for val in InfSDP.known_moments.values()]):
        vcrit = max_within_feasible(InfSDP, InfSDP.known_moments, "dual")
        if abs(vcrit - 1) > 1e-3:
            InfSDP.reset("all")
            v = np.ceil(vcrit * 1000) / 1000
            InfSDP.set_distribution(prob_2in(v, measurements))
            InfSDP.solve(feas_as_optim=True)
            if InfSDP.solution_object["status"] == "feasible":
                print(f"Critical visibility for {info} and {input_names} is ",
                      vcrit)
                cert = InfSDP.certificate_as_probs()
                export_inequality(cert,
                            f"{ineq_path}/twoin_{input_names}_INF2{info}.csv")
            else:
                print("Problem for " + input_names + " failed with status "
                      + InfSDP.solution_object["status"])

        else:
            print(f"Critical visibility for {input_names} is 1")
    else:
        print("The measurements " + input_names +
              " give a trivial distribution at this level")
