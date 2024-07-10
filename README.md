# flowFun
R package containing functions for analysis of flow cytometry data.

Priority:
- change column names of data table to include both markers and channels (specifying
markers instead of channel names to cluster on does not work with data table
currently)
- make output of doDEAnalysis() more user friendly, particularly the `data` component
- reincorporate weights
- make incorporating FMO and isotype controls more straightforward
- plotGroupUMAPs() and plotLabeled2DScatter() each have a case that still needs
  to be accounted for
- generate data for examples in functions used after initial clustering
- links in documentation need to be added or cleaned up
- time permitting, make plotting functions unspecific to FlowSOM
- implement tidytable instead of data.frame for big computations
