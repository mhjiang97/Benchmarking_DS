# Code for benchmarking differential splicing tools at the event level

------------------------------------------------------------------------

## Table of Contents

> -   [simulation](#simulation)
> -   [run](#run)
> -   [figures](#figures)

## simulation

-   **code**: download [SRR493366](https://www.ncbi.nlm.nih.gov/sra/?term=SRR493366), mapping to the genome, and generate simulated fastq files.
-   **files**: files required when simulating fastq files.
-   **profiles**: the output directory of simulated fastq files.
-   \***/single_cell**: code for single-cell data simulation.

## run

-   `benchmark.R` for execution of conventional workflows.
-   `benchmark_j.R` for execution of novel junction workflows.
-   `benchmark_ss.R` for execution of novel splice site workflows.
-   `benchmark_f.R` for execution of transcript filtering workflows.
-   \***/single_cell**: code for single-cell data workflows.

## figures

-   **tidy**: code for result tidy.
-   **plot**: code for figure generation.
-   **rds**: RDS files of each result.
-   \***/single_cell**: code and RDS files for single-cell tool benchmarking.
