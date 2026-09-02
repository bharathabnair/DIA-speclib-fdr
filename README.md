# DIA spectral-library FDR analysis

Analysis scripts used to study how spectral libraries affect false discovery
rate (FDR) estimation in data-independent acquisition (DIA) proteomics.

This repository contains research code rather than an installable software
package. Input data and generated results are not included.

## Files

- `estimate_fdp.R`: prepares DIA-NN reports for FDRBench, estimates false
  discovery proportions (FDP), and provides plotting functions.
- `batch_run_fdp.R`: finds DIA-NN result files, runs the
  FDP analysis in batches, and creates summary plots.
- `properties_entrapped.ipynb`: compares sequence and physicochemical
  properties of target and entrapped peptide databases.
- `calculate_inflation.R`: compares discoveries at selected FDR and FDP
  thresholds and creates inflation plots.

The numerical prefixes reflect the order of these scripts in the larger study;
not every step of that study is included here.

## Requirements

- R with `tidyverse`, `readr`, `ggplot2`, `ggpubr`, `DBI`, `here`, `fs`,
  `purrr`, and `arrow`
- Python with `numpy`, `pandas`, `matplotlib`, and optionally `scipy`
- Java 11 or newer and [FDRBench](https://github.com/Noble-Lab/FDRBench)
- DIA-NN report files and the corresponding target/entrapment databases

Before running the FDP scripts, set the path to the FDRBench JAR:

```bash
export FDRBENCH_JAR=/path/to/fdrbench-1.1.1.jar
```

## Expected data layout

The scripts expect DIA-NN outputs below `output/diann/` and reference peptide
databases below `data/reference-fasta/`. See the path construction in each
script for the expected filenames. These files must be supplied by the user.

## Usage

Run the batch analysis from the repository root:

```bash
Rscript batch_run_fdp.R \
  "$PWD" \
  "$PWD/data/reference-fasta" \
  "$PWD/estimate_fdp.R" \
  peptide
```

Use `protein` instead of `peptide` for protein-level analysis.

Calculate discovery inflation after the FDP result files have been generated:

```bash
Rscript calculate_inflation.R "$PWD/output/diann" "$PWD/results"
```

Open `properties_entrapped.ipynb` from the repository root so its relative
data paths resolve correctly.

## Acknowledgements and third-party code

Parts of the analysis code in this repository were adapted from
[FDRBench](https://github.com/Noble-Lab/FDRBench), developed by the Noble Lab
and contributors. FDRBench is distributed under the
[Apache License 2.0](https://github.com/Noble-Lab/FDRBench/blob/main/LICENSE).
The adapted code was modified for the DIA-NN analyses presented in this
repository.

## License

Licensed under the Apache License 2.0. See `LICENSE`.
