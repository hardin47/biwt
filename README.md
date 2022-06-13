# biwt: Compute the Biweight Mean Vector and Covariance & Correlation Matrices


### Overview
The `biwt` package includes a base function to compute multivariate location, scale, and correlation estimates based on Tukey's biweight M-estimator.  Using the base function, the computations can be applied to a large number of observations to create either a matrix of biweight distances or biweight correlations. 



In its current state, the main functions in this package are:

- `biwt.cor`,
- `biwt.cor.matrix`,
- `biwt.dist.matrix`, and
- `biwt.est`.

### Installation

Running the following lines of code in `R` will install this package from Github:

```{r}
library(devtools)
devtools::install_github(repo = "hardin47/biwt")
```  

### Instructions
See `documentation.pdf` for information on how to use this package. A portion of the example given in the documentation is reproduced below for convenience.

```{r}
set.seed(4747)
samp.data <- MASS::mvrnorm(30, mu=c(0,0,0),
                           Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1), 
                          ncol=3))

r<-0.2 # breakdown

biwt.cor(samp.data[,1:2], r=.2)$biwt.cor
biwt.cor(samp.data[,c(1,3)], r=.2)$biwt.cor
biwt.cor(samp.data[,c(2,3)], r=.2)$biwt.cor


biwt.cor(samp.data, r=.2)

biwt.cor.matrix(samp.data, r=.2)

biwt.dist.matrix(samp.data, r=.2, absval=TRUE)
biwt.dist.matrix(samp.data, r=.2, absval=FALSE)
```

### License
See `DESCRIPTION` for information.

### Author
Johanna Hardin

### References
* Hardin, J., Mitani, A., Hicks, L., VanKoten, B.; A Robust Measure of Correlation Between Two Genes on a Microarray, *BMC Bioinformatics, 8*:220; 2007. https://doi.org/10.1186/1471-2105-8-220
