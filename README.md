# Code for Benchmarking Differential Splicing Tools at the Event Level

<font size="2"> _public repository containing code for reproducing the published results_ </font>

------------------------------------------------------------------------

## Table of Contents

> -   [Publication](#publication)
> -   [Simulation](#simulation)
> -   [Run](#run)
> -   [Figures](#figures)


## Publication

[**A comprehensive benchmarking of differential splicing tools for RNA-seq analysis at the event level**](https://doi.org/10.1093/bib/bbad121) by Jiang et al.

## Simulation

-   **code**: download [SRR493366](https://www.ncbi.nlm.nih.gov/sra/?term=SRR493366), map to the human genome, and generate simulated fastq files.
-   **files**: files required when simulating fastq files.
-   **profiles**: output directories of simulated fastq files.
-   \***/single_cell**: code for single-cell data simulation.

## Run

-   `benchmark.R` for execution of conventional workflows.
-   `benchmark_j.R` for execution of novel junction workflows.
-   `benchmark_ss.R` for execution of novel splice site workflows.
-   `benchmark_f.R` for execution of transcript filtering workflows.
-   \***/single_cell**: code for single-cell data workflows.

## Figures

-   **tidy**: code for tidying results.
-   **plot**: code for generating figures.
-   **rds**: RDS files of each result.
-   \***/single_cell**: code and RDS files for single-cell tool benchmarking.
