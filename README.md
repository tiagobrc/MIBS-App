# MIBS-App

A Shiny app for exploring Multiwell Incidence-based sequencing (MIBS) BCR
datasets. It reads the `sulk1` dataset shipped with the
[MIBS](https://github.com/victoraLab/MIBS) package.

The code dates from 2020-21 (the git history starts in April 2020) and was
recovered in November 2026.

## What it does

**Levenshtein tab**

- Filters: minimum reads per sequence, minimum number of clone-wells a clone
  spans, collapse in-mouse clones, V gene.
- CDR3 amino-acid Levenshtein distance from a chosen pattern (default
  `CAREGFAYW`), plotted as counts or ratio and split by any suitable column.
- Plain or histogram table (Excel download), boxplots, and a Dunn
  (Kruskal-Wallis) multiple-comparison test.
- Plot download as SVG.

**Alignment tab**

- A random sample of sequences, optionally filtered by V and J gene, is
  translated and aligned with Clustal Omega and shown in an interactive
  viewer (`msaR`).
- Nucleotide pattern search with an amino-acid lookahead, and a clonality
  prediction that reports the sequences the pattern missed.

## Requirements

- R >= 4.3 (developed against 4.5).
- **MIBS >= 0.0.0.9004.** Older versions fail on R >= 4.3 and do not import
  `stringdist`; the app stops at start-up with install instructions if the
  version is too old.

  ```r
  # the package lives on this branch until it is merged to master
  remotes::install_github("victoraLab/MIBS@package-conversion")
  ```

- CRAN: `shiny`, `shinythemes`, `bslib`, `plotly`, `DT`, `ggplot2`,
  `stringdist`, `FSA`, `tidyr`, `dplyr`, `htmlwidgets`, `msaR`, `bios2mds`,
  `openxlsx`, `ggpubr`, `ape`, `scales`.
- Bioconductor: `Biostrings`, `msa`.

## Run

```r
shiny::runApp("path/to/MIBS-App")
```

## Known issues

- The database selector offers only `sulk1`; the upload box is commented out.

## License

MIT, see [LICENSE](LICENSE).
