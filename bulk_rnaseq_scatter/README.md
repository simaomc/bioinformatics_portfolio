# Bulk RNA-seq Scatter Plot Analysis

## Objective
This project aims to compare gene expression levels between different stages of spermatogenesis (mitosis, meiosis, and post-meiosis) using bulk RNA-seq data.

## Background
Understanding gene expression dynamics during spermatogenesis is essential to identify stage-specific transcriptional programs, especially post-meiotic gene activation.

## Methods
- Normalized RNA-seq count data
- Log-scale transformation
- Scatter plot visualization
- Highlighting selected gene sets
- Median expression comparison

## Data
The dataset contains normalized counts from:
- 4 mitotic replicates
- 4 meiotic replicates
- 3 post-meiotic replicates

## Tools
- Python (pandas, matplotlib)
- R (optional)

## Output
- Scatter plots comparing expression levels
- Highlighted genes of interest
- Median expression lines

## How to run
1. Load normalized counts table
2. Run the script in `/scripts`
3. Generate plots in `/results`

## Notes
This analysis is part of a broader study on post-meiotic gene expression in Drosophila spermatogenesis.
