## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## -----------------------------------------------------------------------------
# Load data and prepare offset variable
library(bonsaiforest2)
shrink_data <- bonsaiforest2::shrink_data

# Add log follow-up duration as the offset (accounts for varying exposure time)
shrink_data$log_fup_duration <- log(shrink_data$fup_duration)

