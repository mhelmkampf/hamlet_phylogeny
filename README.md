# README – Code and data package

This repository contains data and code used for the manuscript titled **"Radiation with reproductive isolation in the near-absence of phylogenetic signal"** by Helmkampf, Coulmance, et al. published in *Science Advances* (2025). Larger datasets, particularly VCF files, are available in the accompanying **Dryad repository**: [DOI 10.5061/dryad.gxd2547z7](https://doi.org/10.5061/dryad.gxd2547z7). The raw sequencing data can be found at the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home) under the project accession numbers PRJEB75494, PRJEB74501, and PRJEB76552.

## Repository Structure

The repository is organized as follows:

```plaintext
.
├── code       # Scripts used for workflows and analyses
├── data       # Intermediate data required for the analyses
├── figures    # Figures and associated plotting scripts
├── metadata   # Metadata, including sample information
└── results    # Final output files generated from analyses (other than plots)
```

### Code Language and File Extensions

Scripts in this repository use the following programming languages, as indicated by the file extensions:
- `.sh`  : Shell script
- `.py`  : Python
- `.R`   : R

## Usage Notes

This repository serves to document the **parameters and procedures** used in the analyses presented in the paper. However, re-running the code may not be straightforward for users on different systems, due to the following reasons:

- The code was developed over several years by different authors and executed on various systems.
- Code blocks were often run as **sequential Slurm scripts** with manual checks between steps, and were not designed to work as continuously automated workflows.

### Potential Issues

1. **Paths**: Most file paths are system-specific and do not match the file structure of this repository.
2. **Sample Names**: Several samples have been renamed to reflect updated species/population assignments. Older names may persist in some scripts or data files. See `metadata/relabel.txt` for details on these changes.
3. **Dataset and File Names**: Some dataset and file names have changed over the course of the project, and occasionally older names may still appear in some scripts (see table below).
4. **Prerequisite Files**: Providing prerequisite files may have been overlooked in individual cases. Please contact the corresponding author to request any missing files.

### Dataset and File Name Changes

Here are some common changes to dataset names throughout the analysis pipeline:

| Final name                | Working name               |
|------------------------------|-----------------------------|
| `phylo-snp_5kb`               | `phylo2e_m2k5`              |
| `phylo-all_mtg`               | `phylo2_mtg`                |
| `phylo-snp_phased`            | `phylo2e_phased`            |
| `phyps-snp_LDfilt`     | `phyps2e_m2n1l5`     |
| `phyps-snp`            | `phyps2e_m2`         |

## Acknowledgments

If you encounter any issues or require additional information about the dataset or code, please contact the corresponding author (martin.helmkampf@leibniz-zmt.de).
