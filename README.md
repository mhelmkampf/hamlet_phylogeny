# README

This repository contains the data and code used for the manuscript titled **"Radiation with reproductive isolation in the near-absence of phylogenetic signal"** by **Helmkampf, Coulmance, et al.**, published in *Science Advances* (2025). Larger datasets, particularly VCF files, are available in the accompanying **Dryad repository**: [DOI 10.5061/dryad.gxd2547z7](https://doi.org/10.5061/dryad.gxd2547z7).

## Repository Structure

The repository is organized as follows:

```plaintext
.
├── code       # Scripts used for workflows and analyses
├── data       # Intermediate data required for the analyses
├── figures    # Figures and associated plotting scripts
├── metadata   # Metadata, including sample information
└── results    # Final output files generated from analyses
```

### Code Language and File Extensions

Files in this repository use the following programming languages, as indicated by the file extensions:
- `.sh`  : Shell script
- `.py`  : Python
- `.R`   : R

## Usage Notes

This repository serves to document the **parameters** and **procedures** used in the analyses presented in the paper. However, **re-running the code** may not be straightforward for users on different systems, due to the following reasons:

- The code was developed over several years by different authors and executed on various systems.
- Code blocks were often run as **sequential Slurm scripts** with **manual checks** between steps, and may not function as fully automated workflows.

### Potential Issues

1. **System-Specific Paths**: Most file paths are system-specific and do not match the file structure of this repository.
2. **Inconsistent Dataset and File Names**: The dataset and file names have changed multiple times, and older names may still appear in some scripts (see list of common new vs. old names below).
3. **Sample Renaming**: Several samples have been renamed to reflect updated species/population assignments. Some older names may persist in scripts. See the file `metadata/relabel.txt` for details on these changes.
4. **Missing Prerequisite Files**: Some prerequisite files may be missing. Please contact the corresponding author at martin.helmkampf@leibniz-zmt.de if you need any missing files.

## Common File Name Changes

Here are some common changes to dataset names throughout the analysis pipeline:

| Old File Name                | New File Name               |
|------------------------------|-----------------------------|
| `phylo-snp_5kb`               | `phylo2e_m2k5`              |
| `phylo-all_mtg`               | `phylo2_mtg`                |
| `phylo-snp_phased`            | `phylo2e_phased`            |
| `phyps-snp_LDfilt`     | `phyps2e_m2n1l5`     |
| `phyps-snp`            | `phyps2e_m2`         |

## Acknowledgments

If you encounter any issues or require additional information about the dataset or code, please contact the corresponding author.
