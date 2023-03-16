# Code for benchmarking differential splicing tools at the event level

<font size="2"> _Public repository containing code for reproducing the research_ </font>

------------------------------------------------------------------------

## Table of Contents

> -   [publication](#publication)
> -   [simulation](#simulation)
> -   [run](#run)
> -   [figures](#figures)


## pubication

[**A comprehensive benchmarking of differential splicing tools for RNA-seq analysis at the event level**](https://doi.org/10.1093/bib/bbad121) by Jiang et al.

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
