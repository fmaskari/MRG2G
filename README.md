MR-G2G

This R code implements MR-G2G, a two-stage Mendelian randomization framework for estimating the causal effects of CpG methylation sites on a continuous outcome, using SNP genotypes as instruments. It exploits the gene and region structure of the data: SNPs are grouped into linkage-disequilibrium (LD) blocks and CpG sites into regions within genes, and both groupings are built into the instrument-selection penalty.
It supports simulation of correlated genotype, methylation, confounder, and outcome data, cross-validated tuning of the Stage 1 penalty, Stage 2 causal effect estimation with several estimators, and a regional Wald test for which CpG regions affect the outcome. The current release covers individual-level data with a continuous outcome, for a single gene or multiple genes.

![MR-G2G overview](./Fig1.jpg)
