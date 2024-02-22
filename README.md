## Code to accompany *[Guarantees on the network structure of experimental quantum networks](https://www.arxiv.org/abs/2310.xxxxx)*
#### Andrés Ulibarrena, Jonathan W. Webb, Alexander Pickston, Joseph Ho, Christopher L. Morrison, Peter Barrow, Francesco Graffitti, Alessandro Fedrizzi, and Alejandro Pozas-Kerstjens

This is a repository containing the computational appendix of the article "*Guarantees on the network structure of experimental quantum networks*. Andrés Ulibarrena, Jonathan W. Webb, Alexander Pickston, Joseph Ho, Christopher L. Morrison, Peter Barrow, Francesco Graffitti, Alessandro Fedrizzi, and Alejandro Pozas-Kerstjens [arXiv:2403.xxxxx](https://www.arxiv.org/abs/2403.xxxxx)." It provides the codes for obtaining the results depicted in the figures in the manuscript, and the Bell-like inequalities found.

The code is written in Python.

Libraries required:

- [inflation](https://www.github.com/ecboghiu/inflation) (and its requirements) for setting up and solving the compatibility problems.
- [matplotlib](https://matplotlib.org) and [seaborn](https://seaborn.pydata.org) for plots.
- [numpy](https://www.numpy.org) for math operations.
- [pandas](https://pandas.pydata.org) for operations with dataframes in the reading and writing of `.csv` files.
- [qutip](https://qutip.org) for operations with quantum states and measurements.
- [sympy](https://www.sympy.org) and [symengine](https://pypi.org/project/symengine/) for symbolic computations. Sympy is needed for the inflation codes, while symengine is faster when it comes to symbolic manipulation.
- [tqdm](https://tqdm.github.io/) for progress bars.
- [csv](https://docs.python.org/3/library/csv.html), [itertools](https://docs.python.org/3/library/itertools.html), [numbers](https://docs.python.org/3/library/numbers.html), [os](https://docs.python.org/3/library/os.html), [pickle](https://docs.python.org/3/library/pickle.html).

Files and folders:

  - [find_inequalities_ghz.py](https://github.com/apozas/yyyy/blob/main/find_inequalities_ghz.py): Code for obtaining the inequalities for the network corresponding to the experiment of [arXiv:2311.14158](https://arxiv.org/abs/2311.14158).

  - [find_inequalities_trident.py](https://github.com/apozas/yyyy/blob/main/find_inequalities_trident.py): Code for obtaining the inequalities for the network corresponding to the experiment of [npj Quantum Inf. **9**, 82 (2023)](https://doi.org/10.1038/s41534-023-00750-4).

  - [evaluate_inequalities.py](https://github.com/apozas/yyyy/blob/main/evaluate_inequalities.py): Code for obtaining the inequalities for the network corresponding to the experiment of [arXiv:2311.14158](https://arxiv.org/abs/2311.14158).

  - [plots.py](https://github.com/apozas/yyyy/blob/main/plots.py): Code for creating the figures in the manuscript.

  - [utils.py](https://github.com/apozas/yyyy/blob/main/utils.py): Helper functions for importing, exporting and evaluating the inequalities.

  - [GHZInequalities](https://github.com/apozas/yyyy/blob/main/GHZInequalities): Folder with the inequalities obtained for the network corresponding to the experiment of [arXiv:2311.14158](https://arxiv.org/abs/2311.14158).

  - [TridentInequalities](https://github.com/apozas/yyyy/blob/main/TridentInequalities): Folder with the inequalities obtained for the network corresponding to the experiment of [npj Quantum Inf. **9**, 82 (2023)](https://doi.org/10.1038/s41534-023-00750-4).

#### How to read the inequalities
Each inequality is stored in a separate `.csv` file, with 13 columns each. The structure is the following:

  - Column 1: coefficient of the term.

  - Columns 2-7: outputs of the measurements performed by the parties.

  - Columns 8-13: inputs of the measurements performed by the parties.

The value -1 in the columns represents the fact that the party does not play a role in the corresponding term. For example, the term ![](https://latex.codecogs.com/svg.latex?p_{BD}(0,0|1,0)) is represented by the values `-1, 1, 0, -1, 0, -1, -1, 1, -1, 0, -1, -1` in columns 2-13.

Rows are arranged by blocks. Each block begins in a row with a nonempty first column, and finishes in the row before the next row with a nonempty first column. Each row represents a probability, and all probabilities in the block are multiplied. For example, the block

| -0.5 | -1 | 0  | 0  | 0  | -1 | -1 | -1 | 0  | 1  | 0  | -1 | -1 |
|------|----|----|----|----|----|----|----|----|----|----|----|----|
|      | 0  | 0  | -1 | -1 | -1 | -1 | 0  | 0  | -1 | -1 | -1 | -1 |
|      | -1 | -1 | -1 | -1 | -1 | 0  | -1 | -1 | -1 | -1 | -1 | 1  |
|      | -1 | -1 | -1 | -1 | -1 | -1 | -1 | -1 | -1 | -1 | -1 | -1 |

corresponds to the term ![](https://latex.codecogs.com/svg.latex?0.5p_{BCD}(0,0,0|0,1,0)p_{AB}(0,0|0,0)p_{F}(0|1)). The inequality is built by summing all the blocks, and it witnesses incompatibility if the evaluation on a given distribution is negative. The functions `import_ineq`, `read_ineq` and `eval_ineq` in [utils.py](https://github.com/apozas/yyyy/blob/main/utils.py) can be used for manipulating the inequalities and give them a human-readable form.

#### Citing
If you would like to cite this work, please use the following format:

A. Ulibarrena, J. W. Webb, A. Pickston, J. Ho, C. L. Morrison, P. Barrow, F. Graffitti, A. Fedrizzi, and A. Pozas-Kerstjens, _Guarantees on the network structure of experimental quantum networks_, arXiv:2403.xxxxx

```
@misc{ulibarrena2024guarantees,
  title = {Guarantees on the network structure of experimental quantum networks},
  author = {Ulibarrena, Andrés and Webb, Jonathan W. and Pickston, Alexander and Ho, Joseph and Fedrizzi, Alessandro and Pozas-Kerstjens, Alejandro},
  archivePrefix = {arXiv},
  eprint = {2403.xxxxx},
  year = {2024}
}
```
