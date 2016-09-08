# Outlier Analysis
Outlier analysis module to identify aberrantly highly expressed genes.
Takes in a Gene(rows) X Sample(columns) matrix and identifies genes that have outlier expression in the cohort

### Usage
Rscript --vanilla outlier.R -l gene_list -m gene_matrix -o output_directory 
Rscript --vanilla outlier_scaled.R -l gene_list -m gene_matrix -o output_directory 

# using example files
Rscript --vanilla outlier.R -l gene_list_example.txt -m gene_matrix_example.tsv -o out
# with log transformation
Rscript --vanilla outlier.R -l gene_list_example.txt -m gene_matrix_example.tsv -o out -t 1 
# with outlier score cut-off: default 1.5
Rscript --vanilla outlier.R -l gene_list_example.txt -m gene_matrix_example.tsv -o out -t 1 -c 1.5
