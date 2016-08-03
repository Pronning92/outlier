# Outlier Analysis
Outlier analysis module to identify aberrantly highly expressed genes.
Takes in a Gene(rows) X Sample(columns) matrix and identifies genes that have outlier expression in the cohort

### Usage
	Rscript --vanilla outlier.R -l gene_list -m gene_matrix -o output_directory 
	Rscript --vanilla outlier_scaled.R -l gene_list -m gene_matrix -o output_directory 
