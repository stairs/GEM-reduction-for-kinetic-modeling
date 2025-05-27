# Genome Scale Metabolic Model Contextualisation and Reduction

## Overview
This repository provides MATLAB-based workflows to:
1. **Extract** an HT-29 cell line–specific genome-scale metabolic model from the generic Human-GEM using DepMap RNAseq data.
2. **Reduce** the extracted GEM by pruning (and optionally compressing) the network while preserving key metabolic functions.

## Installation
1. Clone this repository
2. Ensure MATLAB R2020b or newer is installed.
3. Install RAVEN Toolbox v2.7.5 or newer: [GitHub](https://github.com/SysBioChalmers/RAVEN).
4. Download and install CellNetAnalyzer: [Official Site](https://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html).
5. Clone Human-GEM repo (we used v1.18): [GitHub](https://github.com/SysBioChalmers/Human-GEM).
6. Create a `.env` file in the project root containing:
   ```ini
   HUMAN_GEM_PATH=/path/to/Human-GEM
   ```

## Scripts Description
- **`extract_gem.m`**: Extracts cell line–specific GEM from Human-GEM.
- **`reduce_gem.m`**: Prunes/compresses the extracted GEM while preserving critical functions.

## Data
- **`ht29-depmap-rnaseq-tpm.csv`**: TPM matrix (genes × samples) for HT-29 from DepMap v23Q4.
- **`medium.csv`**: Exchange reaction IDs and lower bounds used for medium definition in reduction.
