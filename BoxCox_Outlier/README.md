# BoxCox Outlier Analysis #

This outlier analysis pipeline is specially suitable for calculating outlier Z-scores for skewed distributions. The pipeline consists of
two main stateps:

1. Perform a BoxCox transformation on the data vector after calculating an optimized parameter-lambda.
2. Calculate Z-scores for the BoxCox transformed data set.

### Prerequisites ###

```
install.packages("EnvStats")
```

### Example ###

```
Rscript Simple_BoxCox_Outlier_Analysis.R BoxCox_example1.txt
```
