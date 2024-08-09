# flowFun
R package containing functions for analysis of flow cytometry data. This package can be installed from the Rstuio console as follows: 
```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("00berst33/flowFun")
```

Priority:
- change plotBeforeAfter() to draw gates on flowDensity plots instead of 
highlighting cells
- improve automated gating of preprocessing, and allow for more gating schemes
- allow for a single sample's gates to be redrawn manually after failed automated attempt
- relational database for data management?
- make output of doDEAnalysis() more user friendly, particularly the `data` component
- reincorporate weights
- plotGroupUMAPs() and plotLabeled2DScatter() each have a case that still needs
  to be accounted for
- time permitting, make all plotting functions unspecific to FlowSOM
