This repo contains the data and scripts to reproduce the results and figures from the manuscript https://doi.org/10.1101/586362.

You can get the figures and a report in HTML format by running:

**WARNING:** this takes several days to run on my computer (4 cores @2.80GHz and 32Go memory)

```r
# install.packages("remotes")
# remotes::install_github("richfitz/remake")
remake::make()
```
