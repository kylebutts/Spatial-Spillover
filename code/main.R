#' # Main script to run the complete project
#' 
#' `render_file` will run the R file and log the results of the script 
#' in the `logbook`. The logbook can be viewed with 
#' `quarto preview logbook` from the terminal.
library(here)
source(here("logbook/render_file.R"))

# Examples ---------------------------------------------------------------------
0 && render_file(here("code/rings-example/rings_example.R"))

# TVA --------------------------------------------------------------------------
1 && render_file(here("code/TVA/analysis.R"))

# CHC --------------------------------------------------------------------------
1 && render_file(here("code/chc/analysis.R"))

# Opportunity Zones ------------------------------------------------------------
0 && render_file(here("code/OZ/replication.R"))
0 && render_file(here("code/OZ/spillovers.R"))

# Simulations ------------------------------------------------------------------
0 && render_file(here("code/simulations/misspecified_exposure_mapping.R"))



