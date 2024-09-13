## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 5,
  dpi = 75
)
library(biwt)
library(colorspace)
library(gplots)
library(dendextend)
library(tidyverse)

## -----------------------------------------------------------------------------
set.seed(4747)
# The `contamination` function generates a dataset with `p` percent
# contaminated data.
# `n`: number of observations
# `p`: percent of contaminated observations
# `pc`: correlation of the population from which the data are generated
# `u1x`: lower bound of uniform dist on x-axis (contaminated data)
# `u2x`: upper bound of uniform dist on x-axis (contaminated data)
# `u1y`: lower bound of uniform dist on y-axis (contaminated data)
# `u2y`: upper bound of uniform dist on y-axis (contaminated data)
contamination <- function(n = 100, p = 0.1, pc = 0.75, 
                          u1x = -3, u2x = -1, u1y = 15, u2y = 25){
    temp2 <- rbind(data.frame(MASS::mvrnorm(n*(1-p), mu = c(0,0), 
                          Sigma = matrix(c(1,pc,pc,1), ncol= 2))), 
               cbind(X1 = runif(n*p, u1x, u2x), X2 = runif(n*p, u1y, u2y)))
    return(data = temp2)
}


# The `correlation_func` function simulates contaminated data
# and returns both the biweight and the Pearson correlation
# `r`: the breakdown (% of allowed contamination)
correlation_func <- function(r = 0.2, n = 100, p = 0.1, pc = 0.75, 
                             u1x = -3, u2x = -1, u1y = 15, u2y = 25){
  
  samp_data <- contamination(n = n, p = p, pc = pc, 
                             u1x = u1x, u2x = u2x, u1y = u1y, u2y = u2y)

  return(data.frame(bwcorrelation = biwt_cor(samp_data, r)$biwt_corr, 
                    pcorrelation = cor(samp_data)[1,2]))

}

# First, run the correlation function on clean (uncontaminated) data
correlation_func(p = 0) |> 
  gt::gt(caption = "Biweight and Pearson correlations from a sample of 100 observations in an uncontaminated bivariate normal distribution with population correlation of 0.75.") |> 
  gt::cols_label(bwcorrelation = "biweight correlation", 
             pcorrelation = "Pearson correlation") |> 
  gt::tab_options(table.width = gt::pct(85))

## ----fig.cap = "Plot of uncontaminated bivariate normal sample of size 100 from a population with true correlation of 0.75."----
contamination(n = 100, p = 0, pc = 0.75, u1x = -3, u2x = -1, u1y = 15, u2y = 25) |>   
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Uncontaminated dataset") + 
  theme_minimal()

## -----------------------------------------------------------------------------
# using the default values
correlation_func(r = 0.2, n = 100, p = 0.1, pc = 0.75,
                 u1x = -3, u2x = -1, u1y = 15, u2y = 25)|> 
  gt::gt(caption = "Biweight and Pearson correlations from a sample of 100 observations in a 10% contaminated bivariate normal distribution with population correlation of 0.75.") |> 
  gt::cols_label(bwcorrelation = "biweight correlation", 
             pcorrelation = "Pearson correlation") |> 
  gt::tab_options(table.width = gt::pct(85))

## ----fig.cap = "Plot of 10% contaminated bivariate normal sample of size 100 from a population with true correlation of 0.75."----

# using the default values
contamination(n = 100, p = 0.1, pc = 0.75, 
              u1x = -3, u2x = -1, u1y = 15, u2y = 25) |> 
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Dataset with 10% contamination") + 
  theme_minimal()

## ----fig.cap = "One simulated dataset with n = 100, 20% compressed contamination, and a population correlation of 0.75."----
data1 <- contamination(n = 100, p = 0.2, pc = 0.75, 
                      u1x = -3, u2x = -1, u1y = 15, u2y = 25)

data1 |>
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Dataset with 20% compressed contamination") + 
  theme_minimal()

## ----fig.cap = "Histogram of correlations for 500 random samples.  Each of the 500 datasets has n = 100 and 20% compressed contamination. The pink vertical line is at the population correlation, 0.75."----
correlation1 <- 1:500 |> 
  map_dfr(~correlation_func(r = 0.3, n = 100, p = 0.2, pc = 0.75, 
                            u1x = -3, u2x = -1, u1y = 15, u2y = 25))

correlation_long1 <- correlation1 |>
  pivot_longer(cols = c ("bwcorrelation", "pcorrelation"), 
               names_to = "statistic", values_to = "value")

correlation_long1 |>
  mutate(statistic = case_when(
    statistic == "bwcorrelation" ~ "biweight correlation",
    statistic == "pcorrelation" ~ "Pearson correlation"
  )) |> 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  facet_wrap(~statistic) + 
  theme_minimal() +
  labs(x = "correlation value", y = "")

## ----fig.cap = "Scatterplot describing the biweight and Pearson correlations calculated for each of 500 contaminated datasets. The pink vertical and horizontal lines are at the population correlation, 0.75."----
correlation1 |>
  ggplot(aes()) +
  geom_point(aes(x = bwcorrelation, y = pcorrelation)) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  geom_hline(yintercept = 0.75, color = "pink") + 
  theme_minimal() +
  labs(x = "biweight correlation",
       y = "Pearson correlation",
       title = "biweight and Pearson correlation for each \nof 500 simulated datasets")

## ----fig.cap = "One simulated dataset with n = 100, 20% diffuse contamination, and a population correlation of 0.75."----
data2 <- contamination(n = 100, p = 0.2, pc = 0.75, 
                       u1x = -10, u2x = 10, u1y = -10, u2y = 10)

data2 |>
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Dataset with 20% diffuse contamination") + 
  theme_minimal()

## ----fig.cap = "Histogram of correlations for 500 random samples.  Each of the 500 datasets has n = 100 and 20% diffuse contamination. The pink vertical line is at the population correlation, 0.75."----
correlation2 <- 1:500 |> 
  map_dfr(~correlation_func(r = 0.3, n = 100, p = 0.2, pc = 0.75, 
                            u1x = -10, u2x = 10, u1y = -10, u2y = 10))

correlation_long2 <- correlation2 |>
  pivot_longer(cols = c("bwcorrelation", "pcorrelation"), 
               names_to = "statistic", values_to = "val")

correlation_long2 |>
  mutate(statistic = case_when(
    statistic == "bwcorrelation" ~ "biweight correlation",
    statistic == "pcorrelation" ~ "Pearson correlation"
  )) |> 
  ggplot(aes(x = val)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  facet_wrap(~statistic) + 
  theme_minimal() +
  labs(x = "correlation value", y = "")

## ----fig.cap = "Scatterplot describing the biweight and Pearson correlations calculated for each of 500 20% diffuse contaminated datasets. The pink vertical and horizontal lines are at the population correlation, 0.75."----

correlation2 |>
  ggplot(aes()) +
  geom_point(aes(x = bwcorrelation, y = pcorrelation)) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  geom_hline(yintercept = 0.75, color = "pink") + 
  theme_minimal() +
  labs(x = "biweight correlation",
       y = "Pearson correlation",
       title = "biweight and Pearson correlation for each \nof 500 simulated datasets")

## ----fig.cap = "One simulated dataset with n = 100, 40% compressed contamination, and a population correlation of 0.75."----
data3 <- contamination(n = 100, p = .4, pc = 0.75, 
                       u1x = -3, u2x = -1, u1y = 15, u2y = 25)

data3 |>
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Dataset with 40% compressed contamination") + 
  theme_minimal()

## ----fig.cap = "Histogram of correlations for 500 random samples.  Each of the 500 datasets has n = 100 and 40% compressed contamination. The pink vertical line is at the population correlation, 0.75."----
correlation3 <- 1:500 |> 
  map_dfr(~correlation_func(r = 0.3, n = 100, p = .4, pc = 0.75, 
                            u1x = -3, u2x = -1, u1y = 15, u2y = 25))

correlation_long3 <- correlation3 |>
  pivot_longer(cols = c("bwcorrelation", "pcorrelation"), 
               names_to = "statistic", values_to = "value")

correlation_long3 |>
  mutate(statistic = case_when(
    statistic == "bwcorrelation" ~ "biweight correlation",
    statistic == "pcorrelation" ~ "Pearson correlation"
  )) |> 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  facet_wrap(~statistic) + 
  theme_minimal() +
  labs(x = "correlation value", y = "")

## ----fig.cap = "Scatterplot describing the biweight and Pearson correlations calculated for each of 500 40% compressed contaminated datasets. The pink vertical and horizontal lines are at the population correlation, 0.75."----

correlation3 |>
  ggplot(aes()) +
  geom_point(aes(x = bwcorrelation, y = pcorrelation)) +
  geom_vline(xintercept = 0.75, color = "pink") + 
  geom_hline(yintercept = 0.75, color = "pink") + 
  theme_minimal() +
  labs(x = "biweight correlation",
       y = "Pearson correlation",
       title = "biweight and Pearson correlation for each \nof 500 simulated datasets")

## ----fig.cap = "One simulated dataset with n = 100, 20% compressed contamination, and a population correlation of 0.1."----
data4 <- contamination(n = 100, p = 0.2, pc = 0.1, 
                       u1x = -3, u2x = -1, u1y = 15, u2y = 25)

data4 |>
  ggplot() +
  geom_point(aes(x = X1, y = X2)) + 
  labs(x = "first coordinant",
       y = "second coordinant",
       title = "Dataset with 20% compressed contamination") + 
  theme_minimal()

## ----fig.cap = "Histogram of correlations for 500 random samples.  Each of the 500 datasets has n = 100 and 20% compressed contamination. The pink vertical line is at the population correlation, 0.1."----

correlation4 <- 1:500 |> 
  map_dfr(~correlation_func(r = 0.3, n = 100, p = 0.2, pc = 0.1, 
                            u1x = -3, u2x = -1, u1y = 15, u2y = 25))

correlation_long4 <- correlation4 |>
  pivot_longer(cols = c ("bwcorrelation", "pcorrelation"), names_to = "statistic", values_to = "value")

correlation_long4 |>
  mutate(statistic = case_when(
    statistic == "bwcorrelation" ~ "biweight correlation",
    statistic == "pcorrelation" ~ "Pearson correlation"
  )) |> 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0.1, color = "pink") + 
  facet_wrap(~statistic) + 
  theme_minimal() +
  labs(x = "correlation value", y = "")

## ----fig.cap = "Scatterplot describing the biweight and Pearson correlations calculated for each of 500 contaminated datasets. The pink vertical and horizontal lines are at the population correlation, 0.1."----

correlation4 |>
  ggplot(aes()) +
  geom_point(aes(x = bwcorrelation, y = pcorrelation)) +
  geom_vline(xintercept = 0.1, color = "pink") + 
  geom_hline(yintercept = 0.1, color = "pink") + 
  theme_minimal() +
  labs(x = "biweight correlation",
       y = "Pearson correlation",
       title = "biweight and Pearson correlation for each \nof 500 simulated datasets")

## -----------------------------------------------------------------------------
# combining the two datasets and creating 
# a single data frame with the relevant information
votes_repub <- cluster::votes.repub

newer_president <- read_csv("data/1976-2020-president.csv") |> 
  mutate(percent = candidatevotes / totalvotes * 100) |> 
  filter(party_detailed == "REPUBLICAN" & !writein) |> 
  filter(state != "DISTRICT OF COLUMBIA") |> 
  dplyr::select(year, state, percent) |> 
  pivot_wider(id_cols = state, names_from = year, values_from = percent)

all_votes_repub <- cbind(votes_repub, newer_president) |> 
  dplyr::select(-state, -X1976) 
years <- as.numeric(gsub("X", "", colnames(all_votes_repub)))

# For states that did not exist or vote during certain years
# (i.e. Alaska and Hawaii), the NA data is ignored.
dend_NA <- votes_repub |> 
  is.na() |>
  dist() |> 
  hclust() |> 
  as.dendrogram() |> 
  ladderize()

color_gradient <- colorspace::diverge_hcl

## ----fig.width = 9------------------------------------------------------------
dist.repub <- all_votes_repub |> 
  biwt.cor(r = 0.2, absval = FALSE, output = "distance") |>
  as.dist() |>
  usedist::dist_setNames(rownames(votes_repub)) |>
  hclust(method = "complete") |> 
  as.dendrogram() |>
  plot(ylab = "biweigtht correlation distance")

## ----fig.width = 9------------------------------------------------------------
dist.repub <- all_votes_repub |> 
  dist()|>
  hclust(method = "complete") |> 
  as.dendrogram() |>
  plot(ylab = "Euclidean distance")

## ----fig.width = 9------------------------------------------------------------
pcor_dist <- function(data){
  1 - abs(cor(t(data), use = "pairwise.complete.obs"))
}

dist.repub <- all_votes_repub |> 
  pcor_dist() |>
  as.dist() |>
  usedist::dist_setNames(rownames(votes_repub)) |>
  hclust(method = "complete") |> 
  as.dendrogram() |>
  plot(ylab = "Pearson correlation distance")

## ----fig.height = 8, fig.width = 8--------------------------------------------
par(mar = c(1,1,1,1))
bdend <- all_votes_repub |> 
  biwt.cor(r = 0.2, absval = FALSE, output = "distance") |>
  as.dist() |>
  usedist::dist_setNames(rownames(all_votes_repub)) |> 
  hclust(method = "complete") |> 
  as.dendrogram() |> 
  rotate(labels(dend_NA)) |>
  color_branches(k=4)

gplots::heatmap.2(as.matrix(all_votes_repub),
          main = "Votes for\n Republican Presidential Candidate\n (clustered using complete)",
          srtCol = 60,
          dendrogram = "row",
          Rowv = bdend,
          Colv = "NA", # this to make sure the columns are not ordered
          trace="none",          
          key.xlab = "% Votes for Republican\n Presidential Candidate",
          labCol = years,
          denscol = "grey",
          density.info = "density",
          col = color_gradient
         )

## ----fig.height = 8, fig.width = 8--------------------------------------------
par(mar = c(1,1,1,1))
edend <- votes_repub |> 
  dist() |>
  hclust(method = "complete") |> 
  as.dendrogram() |>
  rotate(labels(dend_NA)) |>
  color_branches(k=4)

gplots::heatmap.2(as.matrix(all_votes_repub),
          main = "Votes for\n Republican Presidential Candidate\n (clustered using complete)",
          srtCol = 60,
          dendrogram = "row",
          Rowv = edend,
          Colv = "NA", # this to make sure the columns are not ordered
          trace="none",          
          key.xlab = "% Votes for Republican\n Presidential Candidate",
          labCol = years,
          denscol = "grey",
          density.info = "density",
          col = color_gradient
         )

## ----fig.height = 8, fig.width = 8--------------------------------------------
par(mar = c(1,1,1,1))
pdend <- all_votes_repub |> 
  pcor_dist() |>
  as.dist() |>
  usedist::dist_setNames(rownames(all_votes_repub)) |>
  hclust(method = "complete") |> 
  as.dendrogram() |>
  rotate(labels(dend_NA)) |>
  color_branches(k=4)

gplots::heatmap.2(as.matrix(all_votes_repub),
          main = "Votes for\n Republican Presidential Candidate\n (clustered using complete)",
          srtCol = 60,
          dendrogram = "row",
          Rowv = pdend,
          Colv = "NA", # this to make sure the columns are not ordered
          trace="none",          
          key.xlab = "% Votes for Republican\n Presidential Candidate",
          labCol = years,
          denscol = "grey",
          density.info = "density",
          col = color_gradient
         )

