## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
# Load the main package
library(bonsaiforest2)
library(brms)


## ----ex-data------------------------------------------------------------------
# Load the shrink_data package dataset
shrink_data <- bonsaiforest2::shrink_data

print(head(shrink_data))

