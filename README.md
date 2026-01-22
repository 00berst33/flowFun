# flowFun
R package containing functions for analysis of flow cytometry data. This package can be installed from the Rstudio console as follows: 
```{r, eval = FALSE}
install.packages("BiocManager")
install.packages("devtools")

devtools::install_github("00berst33/flowFun", dependencies = TRUE, build_vignettes = TRUE)
```

User testing:
- Load in user-gated FlowJo workspace, go right to DA and DE testing
  - is R aware of any transformation applied in FlowJo? can it back-transform?
- Try manually editing a few gates after default preprocessing
    - redraw boundary gates to polygon gates, polygon gates to rectangle, etc.
- Try adding additional steps to default gating scheme
- Try exporting gating scheme created in R to FlowJo
- Try to break pipeline; make note of errors


Priority:
- add option to create, use, and edit default .csv for openCyto gating
- normalization
- reduce size of sample data
- make output of doDEAnalysis() more user friendly, particularly the `data` component
- reincorporate weights
- plotGroupUMAPs() and plotLabeled2DScatter() each have a case that still needs
  to be accounted for
