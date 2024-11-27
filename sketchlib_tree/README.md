# Sketch Tree

Sketch Tree is a workflow that utilizes Sketchlib to estimate the core genome and generate a RapidNJ tree. It is designed to process input from assemblies and produce plots

A common workflow model would be
```
input_fastqs
| ASSEMBLER
| collect
| SKETCHLIB_TREE
```
## Key Features
- Automatically estimates the core genome from input assemblies.
- Computes an ANI matrix using Sketchlib.
- Generates and visualizes a RapidNJ tree from the matrix.
## Requirements
This workflow does not require any specific parameters to run.

---

## Included Scripts
The bin directory contains the following scripts:

1. plot_tree.py
    - Description: Plots a tree from an input tree file and saves it as a PNG.
    - Requirements:
        - Python library: ETE3.

2. ani_tree_tools.py
    - Description: Provides tools for working with ANI matrices and RapidNJ.
    - Requirements:
        - Python dependencies: numpy, pandas.

3. tree_builder.py
    - Description: Builds trees using RapidNJ, relying on functions from ani_tree_tools.py.
    - Requirements:
        - Executable: RapidNJ.
        - Python dependencies: dendropy

## Installation and Dependencies
To use this tool effectively, ensure the following are installed and accessible in your environment:

1. Python Libraries:
    - ETE3
    - dendropy
    - numpy
    - pandas

2. External Tools:
    - RapidNJ: https://github.com/somme89/rapidNJ