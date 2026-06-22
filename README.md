# QCombo: Automated Commutator Calculation for Quantum Many-Body Operators

QCombo is a Python library for automated computation of commutators for normal-ordered many-body operators based on the generalized Wick theorem. Designed specifically for nuclear physics and condensed matter applications, it provides efficient and accurate tools for handling complex quantum many-body calculations.

## Table of Contents

- Overview
- Features
- Installation
- Directory Structure
- Quick Start
- Command Line Interface
- Python API Usage
- Index Conventions
- Examples
- Applications
- Testing
- Citation
- License
- Support

## Overview

### Background and Motivation

In quantum mechanics, operator commutators form the foundation of theoretical calculations. In nuclear physics, *ab initio* methods (such as NCSM, CC, IMSRG, etc.) are essential tools for studying nuclear many-body systems. These methods start from nuclear forces constructed from chiral effective field theory and provide approximate but rigorous solutions to quantum many-body problems.

As the number of many-body terms increases, the number of terms in commutators grows exponentially, making manual calculations extremely tedious and error-prone. QCombo addresses this challenge by providing an automated solution for computing commutators of normal-ordered many-body operators.

## Features

- **Based on Generalized Wick Theorem**: Implements automated computation of commutators for normal-ordered many-body operators
- **Many-Body Coupling Support**: Handles commutators of arbitrary-order many-body operators
- **Automated Output Generation**:
  - Generates LaTeX files for commutator expressions
  - Produces input files for the AMC software package to obtain expressions in J-scheme
- **Flexible Configuration**: Allows specification of output to specific many-body ranks
- **Symbolic Computation**: Built on SymPy for precise algebraic manipulations
- **Command Line Interface**: Easy-to-use terminal commands for quick calculations
- **Interactive Mode**: Step-by-step guided interface for exploratory calculations

## Installation

### Install from PyPI (Recommended)

```
pip install qcombo
```

### Install from Source

```
git clone https://github.com/chenlh73/qcombo.git
cd qcombo
pip install -e .
```
Note that for the first execution, running "pip install -e ."n the terminal may be required to activate the qcombo library in the current directory before the examples can be run.

### System Requirements

- Python >= 3.12
- SymPy >= 1.13.3

You can install Sympy from PyPI
```
pip install sympy
```

as you run the jupter notebook, you may need to install ipython libraty
```
pip install ipython
```

## Directory Structure

The directory structure of the qcombo package is organized to clearly separate core code, packaging configurations, and example resources. Below is a breakdown of key components:  

```
qcombo/                        # Root directory of the Python package.
├── qcombo/                    # Core computation modules.
│   ├── __init__.py            # Marks the qcombo directory as a Python package.
│   ├── __main__.py            # Entry point for command-line interface (enables qcombo to run via terminal).
│   ├── canonical.py           # Implements core algorithms for canonical commutator computations.
│   ├── output.py              # Handles generation of LaTeX/AMC-formatted output files.
│   ├── simplify.py            # Provides symbolic simplification utilities using SymPy.
│   ├── tools.py               # Contains helper functions for operator manipulation.
│   └── wickcalculate.py       # Applies the generalized Wick theorem for commutator calculation.
├── examples/                  # Example notebooks and scripts demonstrating various use cases.
│   ├── easyCombo_Guide.ipynb  # Comprehensive guide for using the easyCombo function.
│   ├── MR_IMSRG2.ipynb        # Multi-reference IMSRG(2) flow equation example.
│   ├── SR_IMSRG3.ipynb        # Single-reference IMSRG(3) flow equation example.
│   ├── Brillouin.ipynb        # Demonstration of Brillouin's theorem.
│   ├── Product_1B1B.ipynb     # Product of two one-body operators.
│   ├── easyCombo.py           # Python script example for easyCombo usage.
│   ├── easyCombo_110.py       # Python script example for [1B,1B] commutator coupled to 0-body.
│   ├── qcombo_to_amc.ipynb    # Workflow for converting QCombo results to AMC input.
│   └── results/               # Output files generated from example scripts.
├── test/                      # Test notebooks for validating correctness.
│   ├── FunctionTest.ipynb     # Comprehensive function tests.
│   ├── Jacobi_Identity.ipynb  # Verification of Jacobi identity for commutators.
│   ├── CommutatorAntisymmetry.ipynb  # Tests for antisymmetry properties.
│   ├── Indices_test.ipynb     # Tests for index manipulation functions.
│   └── parallelTest.ipynb     # Tests for parallel computing functionality.
├── results/                   # Example output files from typical calculations.
│   ├── commutator_*.tex       # LaTeX output files for various commutator calculations.
│   └── commutator_*.amc       # AMC input files for various commutator calculations.
├── LICENSE                    # Open-source license file (MIT License).
├── README.md                  # Main documentation file (this document).
├── setup.cfg                  # Configuration file for setuptools (specifies build/distribution rules).
└── setup.py                   # Script for building and distributing the Python package.
```

## Quick Start

### Command Line Interface (Recommended)

After installation, you can use QCombo directly from the command line:

```
# Compute the [2B, 2B] commutator with all possible contractions
qcombo 2 2

# Compute only 0-body and 1-body contractions
qcombo 2 2 -c 0,1

# Compute with range of contractions (0-2 body)
qcombo 1 3 -c 0-2

# Specify custom output filenames
qcombo 2 2 -c 0,1,2 --latex-output my_result.tex --amc-output my_result.txt

# Use specific indices (Python list syntax)
qcombo "[['a','b'], ['c','d']]" "[['e','f','g'], ['h','i','j']]" -c 0,1,2

# Use single-reference Wick's theorem mode
qcombo 2 2 -w SR

# Enable parallel computing
qcombo 2 2 --parallel

# Start interactive mode
qcombo
# or explicitly:
qcombo -i

# Show help message
qcombo --help

# Show version information
qcombo --version
```

### Python API

You can also use QCombo within Python scripts:

```
import qcombo

# Compute the [2B, 2B] commutator coupled to all possible many-body ranks
qcombo.easyCombo(2, 2)

# Output only terms coupled to 0-body and 1-body ranks
qcombo.easyCombo(2, 2, [0, 1])
```

## Command Line Interface

### Basic Usage

QCombo provides a powerful command-line interface that allows you to perform calculations directly from your terminal:

```
qcombo LEFT RIGHT [OPTIONS]
```

Where:
- `LEFT`: Number of particles in left operator (integer) or index list (Python list string)
- `RIGHT`: Number of particles in right operator (integer) or index list (Python list string)

### Command Line Options

| Option | Short | Description | Example |
|--------|-------|-------------|---------|
| `--contraction` | `-c` | Contraction body numbers: single integer(0), range(0-2), list(0,1,2), or 'all' (calculate all) | `qcombo 2 2 -c 0,1,2` |
| `--latex-output` | `-lo` | LaTeX output filename | `qcombo 2 2 -lo result.tex` |
| `--amc-output` | `-ao` | AMC program input file output name | `qcombo 2 2 -ao result.txt` |
| `--wick-mode` | `-w` | Wick's theorem mode: 'SR' (single-reference) or 'MR' (multi-reference). Default: MR | `qcombo 2 2 -w SR` |
| `--no-process` |`-ns` | Disable the process bar of calculation | `qcombo 2 2 --no-process` |
| `--parallel` | `-p` | Enable parallel computing using ProcessPoolExecutor | `qcombo 2 2 --parallel` |
| `--version` | `-v` | Show version information | `qcombo --version` |
| `--quiet` | `-q` | Quiet mode, reduce output (also disables process bar) | `qcombo 2 2 -q` |
| `--interactive` | `-i` | Start interactive mode | `qcombo -i` |

### Interactive Mode

For exploratory calculations, QCombo offers an interactive mode that guides you through the process step-by-step:

```
qcombo
```

Or explicitly:
```
qcombo -i
```

In interactive mode, you will be prompted for:
1. Left operator (integer or index list)
2. Right operator (integer or index list)
3. Contraction body numbers
4. Wick's theorem mode ('SR' or 'MR')
5. Show process bar (y/n)
6. Enable parallel computing (y/n)
7. Output filenames (optional)

### Examples

**Basic calculations:**
```
# 2-body with 2-body operator commutator
qcombo 2 2

# 1-body with 3-body operator commutator
qcombo 1 3

# 2-body with 2-body, only 0-body and 1-body terms
qcombo 2 2 -c 0,1

# Enable parallel computing for large calculations
qcombo 2 3 --parallel
```

**With custom indices:**
```
# Using specific index lists
qcombo "[['a','b'], ['c','d']]" "[['e','f'], ['g','h']]"

# Using range contractions
qcombo "[['i','j'], ['k','l']]" "[['m','n','o'], ['p','q','r']]" -c 0-2
```

**With output files:**
```
# Save results to specific files
qcombo 2 2 --latex-output commutator.tex --amc-output commutator.txt

# Using short options
qcombo 2 2 -c 0,1,2 -lo result.tex -ao result.txt
```

## Python API Usage

### Function Parameters

```
qcombo.easyCombo(left, right, contraction=None, latexOutput=None, amcOutput=None, **kwargs)
```

- `left`: Many-body rank of the first operator (integer) or index list
- `right`: Many-body rank of the second operator (integer) or index list
- `contraction`: Optional parameter specifying output ranks
  - Example: `[0, 1]` outputs only zero-body and one-body terms
  - `None` (default) outputs all possible many-body terms
- `latexOutput`: Custom LaTeX output filename
- `amcOutput`: Custom AMC input filename
- `**kwargs`: Optional keyword arguments
  - `wick_mode`: `'SR'` (single-reference) or `'MR'` (multi-reference). Default: `'MR'`
  - `show_process`: Whether to display progress bar. Default: `True`
  - `parallel`: Enable parallel computing using ProcessPoolExecutor. Default: `False`

### Basic Examples

```
import qcombo

# Compute [2B, 2B] commutator
result = qcombo.easyCombo(2, 2)

# Compute with specific contractions
result = qcombo.easyCombo(2, 2, contraction=[0, 1])

# Compute with custom output files
result = qcombo.easyCombo(2, 2, 
                          contraction=[0, 1, 2],
                          latexOutput="my_commutator.tex",
                          amcOutput="my_commutator.txt")

# Compute with single-reference mode and disable parallel computing
result = qcombo.easyCombo(2, 2, contraction=[0, 1],
                          wick_mode='SR', show_process=False, parallel=False)

# Compute with specific indices
left_indices = [['a', 'b'], ['c', 'd']]
right_indices = [['e', 'f', 'g'], ['h', 'i', 'j']]
result = qcombo.easyCombo(left_indices, right_indices, contraction=[0, 1])
```

### Output Description

After execution, QCombo generates:

1. **Console Output**: Displays the computed commutator expressions in formatted text
2. **LaTeX Files**: Generated expressions in LaTeX format (ready for inclusion in papers)
3. **AMC Input Files**: Input files for the AMC software package to obtain expressions in J-scheme

The output filenames follow the pattern:
- LaTeX: `commutator_[left]B[right]B_to_[contractions]B.tex`
- AMC: `commutator_[left]B[right]B_to_[contractions]B.amc`

For example: `commutator_2B2B_to_0_1_2B.tex`, `commutator_2B2B_to_0_1_2B.amc`

### Output Format

The returned result is a dictionary where keys are strings in the format `'{filter_body}B_lambda{lambda_body}B'` and values are SymPy expressions:

```
{
    '0B_lambda1B': sympy_expression_1,
    '0B_lambda2B': sympy_expression_2,
    # ...
}
```

### Example Output

For example, compute the commutator of two normal-ordered one-body operators contracting to zero-body:

```python
qcombo.easyCombo(1, 1, 0)
```

In the output file `commutator_1B1B_to_0B.tex`, the equation is expected to be:

$$R^{}_{} = \sum_{ab} (n^{}_{a}-n^{}_{b})  G^{a}_{b}  H^{b}_{a}$$

In the output file `commutator_1B1B_to_0B.amc`, the equation is expected to be:

```
declare G{mode= (1,1),latex ="G" } 
declare H{mode= (1,1),latex ="H" } 
declare R0{mode= (0,0),latex ="R" } 
declare n {  mode=2, diagonal=true, latex="n"} 
 
# commutator [1B,1B]-0B 
# lambda_1B
R0 = 1*sum_ab((n_a-n_b)*G_ab*H_ba);
```


## Index Conventions

QCombo uses a specific index notation based on the generalized Wick theorem for normal-ordered many-body operators. Understanding these conventions is essential for correctly interpreting input and output expressions.

### Upper and Lower Indices

In second quantization, a many-body operator is expressed in terms of creation and annihilation operators. QCombo represents the matrix elements of such operators using **upper indices** (creation/particle indices) and **lower indices** (annihilation/hole indices):

```
Operator:  O = Σ  O^{a₁a₂...}_{b₁b₂...}  a†_{a₁} a†_{a₂} ... a_{b₂} a_{b₁}
```

In QCombo's SymPy representation:
```python
O[(a1, a2, ...), (b1, b2, ...)]   # upper tuple, lower tuple
```

**Example**: A two-body operator `G` with upper indices `a, b` and lower indices `c, d`:
```python
G[(a, b), (c, d)]   # represents G^{ab}_{cd}
```

### Normal Ordering and the A Operator

QCombo expresses all results in terms of **normal-ordered products** with respect to a reference state (Fermi vacuum). The normal-ordering symbol is represented by the tensor `A`:

```
A[(a, b), (c, d)]   # normal-ordered product: {a†_a a†_b a_d a_c}
```

The indices of `A` are the **external (free) indices** of the commutator result and are NOT summed over. They define the many-body rank of the output term:
- `A[(), ()]` → 0-body (scalar) term
- `A[(a,), (b,)]` → 1-body term
- `A[(a, b), (c, d)]` → 2-body term

### Special Tensors

The generalized Wick theorem introduces the following contraction tensors:

| Symbol | Name | Meaning |
|--------|------|---------|
| `λ` (lambda) | Irreducible density matrix | `λ^{a₁...}_{b₁...}`,one-body is particle density matrix `λ^i_j = ⟨a†_i a_j⟩`|
| `ξ` (xi) | Hole density | `ξ^i_j = λ^i_j - δ^i_j`, only has one-body|
| `δ` (delta) | Kronecker delta | `δ^a_b = 1` if `a = b`, else `0` |
| `n` | Occupation number | `n_a = ⟨a†_a a_a⟩`, diagonal part of the one-body density |

**In MR (multi-reference) mode**, `λ` can be multi-body (e.g., `λ^{ab}_{cd}` for 2-body density).

**In SR (single-reference) mode**, only the one-body `λ^a_b` is non-zero, 

**In Natural Orbital Basis** the one-body lambda simplifies to:
```
λ^a_b = n_a · δ^a_b
```

### Dummy Indices vs. Free Indices

- **Free indices**: Indices that appear exactly once in a term. These are the indices of the `A` operator and define the output many-body rank.
- **Dummy indices**: Indices that appear exactly twice (once upper, once lower) and are summed over. They can be freely renamed without changing the value of the expression.

**Example**:
```
Σ_p  G^{ap}_{cd} · H^{eb}_{pf}   # p is a dummy index (summed)
                                   # a, b, c, d, e, f are free indices
```

### Antisymmetry

Matrix elements in QCombo are **fully antisymmetric** under exchange of upper indices or lower indices:

```
G^{ab}_{cd} = -G^{ba}_{cd} = -G^{ab}_{dc} = G^{ba}_{dc}
```

This property is automatically exploited during simplification to combine equivalent terms.

### Input Format

When specifying operators by integer rank, QCombo automatically assigns indices from a predefined alphabetical list `['a', 'b', 'c', ...]`:

```python
# left = 2 means a 2-body operator with indices [['a','b'], ['c','d']]
# right = 2 means a 2-body operator with indices [['e','f'], ['g','h']]
qcombo.easyCombo(2, 2)
```

You can also specify custom indices explicitly:
```python
left  = [['i', 'j'], ['k', 'l']]   # upper indices [i,j], lower [k,l]
right = [['p', 'q'], ['r', 's']]
qcombo.easyCombo(left, right)
```

## Examples

QCombo includes a rich collection of example notebooks in the `examples/` directory:

### Getting Started

| Notebook | Description |
|----------|-------------|
| [easyCombo_Guide.ipynb](examples/easyCombo_Guide.ipynb) | Comprehensive guide covering all features of the `easyCombo` function |
| [Product_1B1B.ipynb](examples/Product_1B1B.ipynb) | Simple example: product of two one-body operators |

### Physics Applications

| Notebook | Description |
|----------|-------------|
| [MR_IMSRG2.ipynb](examples/MR_IMSRG2.ipynb) | Multi-reference IMSRG(2) flow equations  |
| [SR_IMSRG3.ipynb](examples/SR_IMSRG3.ipynb) | Single-reference IMSRG(3) flow equations |
| [Brillouin.ipynb](examples/Brillouin.ipynb) | Verification of Brillouin generator |

### Advanced Usage

| Notebook | Description |
|----------|-------------|
| [qcombo_to_amc.ipynb](examples/qcombo_to_amc.ipynb) | Complete workflow from QCombo calculation to AMC input file generation |
| [easyCombo_110.ipynb](examples/easyCombo_110.ipynb) | [1B, 1B] commutator coupled to 0-body term |

### Python Script Examples

For users who prefer Python scripts over Jupyter notebooks:

- [easyCombo.py](examples/easyCombo.py) - Basic usage example
- [easyCombo_110.py](examples/easyCombo_110.py) - [1B, 1B] commutator example

## Applications

### Nuclear Physics Applications

In nuclear *ab initio* methods, QCombo can be used for:

1. **IMSRG Flow Equations**: Computing commutators between generators and Hamiltonians
2. **Similarity Renormalization Group (SRG)**: Evaluating flow equations
3. **Configuration Interaction**: Deriving effective interactions

### Education and Research

QCombo is also valuable for quantum mechanics education and research:
- Verifying theoretical derivations
- Exploring contributions from higher-order many-body terms
- Automating tedious symbolic calculations
- Teaching quantum many-body theory concepts

## Testing

QCombo includes test notebooks in the `test/` directory to validate correctness:

| Notebook | Description |
|----------|-------------|
| [FunctionTest.ipynb](test/FunctionTest.ipynb) | Comprehensive tests for all core functions |
| [Jacobi_Identity.ipynb](test/Jacobi_Identity.ipynb) | Verification that commutators satisfy the Jacobi identity |
| [CommutatorAntisymmetry.ipynb](test/CommutatorAntisymmetry.ipynb) | Tests for antisymmetry properties of commutators |
| [Indices_test.ipynb](test/Indices_test.ipynb) | Tests for index manipulation and simplification functions |
| [parallelTest.ipynb](test/parallelTest.ipynb) | Validation of parallel computing results |

### Running Tests

Open the test notebooks in Jupyter:

```
jupyter notebook test/
```

Run all cells in each notebook to verify correctness. The notebooks include:
- Expected output comparisons
- Algebraic identity verifications
- Parallel vs. serial result consistency checks

## Citation

If you use QCombo in your research, please cite:

```
@misc{chen2026qcombo,
      title={Qcombo: A Python Package for Automated Commutator Calculations of Quantum Many-Body Operators}, 
      author={L. H. Chen and Y. Li and H. Hergert and J. M. Yao},
      year={2026},
      eprint={2603.24399},
      archivePrefix={arXiv},
      primaryClass={nucl-th},
      url={https://arxiv.org/abs/2603.24399}, 
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

- **Issues**: Report bugs or request features on GitHub Issues: https://github.com/chenlh73/qcombo/issues
- **Documentation**: Visit the GitHub repository for more information: https://github.com/chenlh73/qcombo
- **Questions**: Contact the development team via email

---

*QCombo - Quantum Commutator of Many-Body operator *