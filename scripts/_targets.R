# library(tarchetypes) # Load other packages as needed.
library(targets)
devtools::load_all(".")
# library(flowFun)

# Set target options:
tar_option_set(
  packages = c("flowCore", "PeacoQC", "Biobase", "ggplot2", "flowFun",
               "flowDensity", "flowCut"), # Packages that your targets need for their tasks.
  format = "qs"
  # Set other options as needed.
)

# Replace the target list below with your own:
file <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw/Ab_PHA_Ctrl AWB4.fcs"
comp <- system.file("extdata", "compensation_matrix.csv", package = "flowFun")
comp <- read.csv(comp, check.names = FALSE)

flowPreprocessing(file,
                  compensation = comp,
                  ld_channel = "BUV496-A",
                  transformation_type = "logicle")
