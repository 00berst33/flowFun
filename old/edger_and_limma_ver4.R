# Perform differential marker expression analysis of flow cytometry data using 
# limma, and differential count analysis using edgeR. Before using this script, 
# you should have preprocessed all .fcs files, and created a FlowSOM object 
# whose metaclusters you would like to test.
# 
# This script requires the packages flowCore, FlowSOM, limma, and edgeR to run.

library(flowCore)
library(FlowSOM)
library(edgeR)
library(limma)

# Note that this workflow creates a directory structure as follows:
# 
# home ---------- Data ----- SetUp
#           +            +-- Preprocessed SetUp
#           +            +-- Raw
#           +            +-- Preprocessed Raw
#           +            +-- Aggregates
#           +            +-- Clustered Raw
#           +            +-- directories for any specific cell types of interest (e.g. CD4 T cells)
#           +
#           + --- Info
#           +
#           + --- RDS ------ Unedited
#           +            +-- Edited
#           +
#           + --- Preprocessing Results ----- QC plots
#           +
#           + --- Analysis Results ----- limma Analysis ----- Tables
#                                    +-- edgeR Analysis ----- Tables

# In this file, we are only concerned with the directories defined below. You 
# may change their names if you prefer, but do not change the names of the 
# variables. If the directories have already been created, and you would like 
# to rename them, you must do so in your file editor, then edit the 
# variables below accordingly.

# Contains all preprocessed .fcs files. 
dir_prepr = "Data/Preprocessed Raw/"

# Contains all aggregate .fcs files. Only relevant if your FlowSOM object was
# created by clustering on principal components.
dir_agg = "Data/Aggregates/"

# Contains edited FlowSOM objects.
dir_rds_edited = "RDS/Edited/"

# Will contain all plots and tables that are created from analysis.
dir_result = "Analysis Results/"

# Will contain any figures and tables relevant to edgeR analysis.
dir_edger = "Analysis Results/edgeR Analysis/"

# Will contain .csv files containing results of the edgeR analysis.
dir_edger_tables = "Analysis Results/edgeR Analysis/Tables/"

# Will contain all figures and tables relevant to limma analysis.
dir_limma = "Analysis Results/limma Analysis/"

# Will contain .csv files containing the scripts' results.
dir_limma_tables = "Analysis Results/limma Analysis/Tables/"

# Checks whether the directories that will contain results already exist, and 
# creates them if not.
for (dir in c(dir_result, dir_limma, dir_limma_tables, dir_edger, dir_edger_tables)) {
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
}

# Read in FlowSOM object you would like to use for analysis. 
fSOM = readRDS(paste0(dir_rds_edited, "fsom_cd8_pca_edited_A.rds"))

# Get the names of each metacluster.
meta_names = levels(fSOM$metaclustering)

# Define the names of any metaclusters you are uninterested in for this analysis.
# For example, you may have a metacluster that you identified as debris. If you
# do not remember the names of your metaclusters, simply call `meta_names`.
# ("MC" is not part of the metacluster name. So if you call `meta_names` and see
# a cluster called "MC pDC", the metacluster's actual name is just "pDC".)
meta_to_exclude = c("debris")

# Remove metaclusters not of interest, if necessary.
if (length(meta_to_exclude) > 0) {
  meta_names = meta_names[-(which(meta_names %in% meta_to_exclude))]
}

# If you have a .csv file specifying sample information, it can be read in here.
# It MUST have columns for file names and sample names, and any additional 
# columns may specify parameters of interest, like NAC vs. No NAC. Each column 
# should be named after what it represents. This file will be read in as a data 
# frame and will be used to construct the design matrix. 
sample_df = read.csv(file = "Info/MR242 Sample Information.csv",
                     header = TRUE)

# If you do not have a .csv file specifying sample information, you must instead
# manually specify information about each sample. This is done by creating a 
# data frame where each row is a sample, and each column is a parameter of 
# interest. So, constructing this data frame from scratch would look something
# like the following: 
# 
# sample_df = data.frame(sample = c("Ctrl AWB2", "MIBC MR95"),
#                        filename = c("Ab_PHA_Ctrl AWB2", "Ab_PHA_MIBC MR95.fcs"),
#                        sex = c("female", "male"),
#                        disease = c("Ctrl", "MIBC"),
#                        NAC = c("NA", "NAC"))
# 
# In this case, our columns will be sample, filename, sex, disease, and NAC, and 
# we will have two rows corresponding to each sample and its characteristics. 

# Note that the script will not run if any of your parameters' column names or 
# values are not R-friendly. So if there are any spaces, or invalid characters like 
# $, @, !, in your data frame, these should be removed or replaced. You may do
# this quickly by using the function "make.names" on any columns of concern.
# For example, if you have a column titled "Tumor", defining Tumor vs. No Tumor,
# this column may be of concern, as it is not applicable to your control samples,
# and will be considered as NA values by R. You would fix this by calling the 
# following:
# 
# sample_df$Tumor = make.names(sample_df$Tumor)
# 
# Note that all you have to provide is the name of your column, after each $ symbol.
sample_df$NAC = make.names(sample_df$NAC)

# Define the name of the column containing your sample names. This is 
# necessary to ensure that samples and their cells are assigned to the proper 
# groups. 
sample_col = "Sample.Name"

# Define the name of the column containing your file names.
sample_files = "File.Name"

# Move the column containing sample names so that it is the first column in 
# the data frame.
sample_df = sample_df[, c(sample_col, setdiff(names(sample_df), sample_col))]

# Define the columns of "sample_df" that are relevant for analysis. Please make 
# sure that the column names you type here are identical to how they appear in 
# the data frame you defined above. It is expected that at the very least, 
# sample names are included. So, if you were to specify this variable according 
# to the two-sample data frame defined in the example above, you would type:
# 
# vars_of_interest = c("sample", "sex", "disease", "NAC")
# 
# Note that we excluded the parameter "filename", as in terms of the analysis
# it is made redundant by our "sample" column, and is only included to 
# ensure that our data is labeled properly.
vars_of_interest = c("Sample.Name", "Sex", "Disease", "NAC")

# Define the comparisons you are interested in testing during the analysis, as
# a list of lists. Each nested list should include at least two levels of a 
# factor to be used for the comparison. If you would like to compare two groups
# within another group, you should specify the level(s) of interest for both of 
# these factors. For example, if you were interested in control vs MIBC in 
# females only, and male vs female for all samples, you would define the 
# variable as follows: 
# 
# comparisons = list(
#   female_mibc_vs_ctrl = list(Sex = "female", Disease = list("Ctrl", "MIBC")),
#   male_vs_female = list(Sex = list("male", "female"))
# )
# 
# You may name each comparison as you see fit, just ensure that it is R-friendly. 
# The factor names, however (which would be "Sex" and "Disease" in the above 
# example), should be identical to their column names in our "sample_df" data 
# frame. Furthermore, the factor levels (which would be "female", "male", "Ctrl",
# and "MIBC" in the above example), should be listed identical to how they 
# appear in the data frame (AFTER correcting for any non-R-friendly values). 
comparisons = list(
  male_vs_female = list(Sex = list("male", "female")),
  male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female")),
  ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

# Get the names of all preprocessed .fcs files. Although file names are included 
# in the data frame defined above, this variable is included in case any files 
# were discarded during preprocessing due to issues with the sample, like low 
# cell viability. For example, you may have a data frame which contains information
# on all 10 of your samples, however, during preprocessing you discarded 3 samples
# due to too few live cells. By checking the file names contained in this variable
# against the file names contained in the data frame, we are able to determine
# which samples were discarded, and extract only the relevant information from
# the data frame. 
# 
# It is also necessary to include this variable because it reads 
# the preprocessed files in in the same order as the function we used to create 
# our aggregate files. This is significant because the aggregate .fcs files do 
# not preserve information about file names; instead, each file is given a number 
# based on the order it was read in. Therefore, this variable also serves as a 
# way to determine which sample corresponds to which number in the aggregate 
# file, ensuring that our cells aren't attributed to the wrong samples.
prepr_files = list.files(path = dir_prepr)

# Remove samples that were excluded from the analysis from the data frame 
# `sample_df`, and reorder its rows such that they are the same as 
# the file order contained in the `prepr_files` variable.
matched_ind = match(prepr_files, sample_df[[sample_files]])
sample_df = sample_df[matched_ind, ]
rownames(sample_df) = seq(1, length(prepr_files))

# Remove any samples you are uninterested in by defining a vector of row indices.
# Even if you do not wish to remove any samples, please still set `rows_to_remove`
# equal to an empty vector.
rows_to_remove = c()

if (length(rows_to_remove) > 0) {
  sample_df = sample_df[-rows_to_remove, ]
}

#########
# edgeR #
#########

# NOTE !!! If you are running multiple edgeR analyses on the same preprocessed
# files (i.e. first running an analysis on all cell types, then running an analysis
# on a reclustering of CD8 T cells), and are performing the same tests between 
# groups, it is technically not necessary to redefine many of the variables 
# below (setting `normalize` to `TRUE` multiple times may feel redundant).
# However, it is highly recommended that you rerun all of the code below between 
# analyses to ensure that no logic errors occur. If you attempt to guess which 
# variables do not need to be defined again, it may result in you missing 
# variables which are crucial to redefine. By being thorough, you ensure that 
# you avoid errors like using normalization factors from your CD4 T cells for 
# your CD8 T cells.

# Set to TRUE if you would like to normalize your data. The way in which the data
# is normalized is decided by the `norm_factors` variable below.
normalize = TRUE

# By default, this is set to "TMM" to normalize using the "trimmed mean of M-values" 
# (TMM) method. However, if you wish, you may also set it to your own vector of 
# normalization factors, and those will be used instead. Note that the vector's 
# product must be 1.
norm_factors = "TMM"

# Set the minimum number of cells a metacluster should have in a specified number 
# of samples to be included in the analysis.
min_cells = 3

# Set the minimum number of samples a metacluster should have at least `min_cells`
# events in to be included in the analysis. By default, this is half the total
# number of samples.
min_samples = length(prepr_files)/2

if (TRUE) {
  # Create count matrix, where rows are cell types (defined by your FlowSOM object's
  # metaclusters) and columns are samples.
  count_mat = c()
  for (i in 1:length(prepr_files)) {
    name_ind = which(sample_df[[sample_files]] == prepr_files[i]) 
    
    if (length(name_ind) > 0) {
      ind = which(fSOM$data[, "File"] == i)
      meta_assignments = GetMetaclusters(fSOM)[ind]
      
      meta_counts = c()
      for (j in meta_names) {
        count = length(which(meta_assignments == j))
        meta_counts = c(meta_counts, count)
      }
      
      count_mat = cbind(count_mat, meta_counts)
      num_col = ncol(count_mat)
      
      if (num_col == 2) {
        name_ind = which(sample_df[[sample_files]] %in% prepr_files[c(1,2)])
        colnames(count_mat) = sample_df[name_ind, sample_col]
      } else if (num_col > 2) {
        colnames(count_mat)[num_col] = sample_df[name_ind, sample_col]
      }
    }
  }
  
  # Name rows of count matrix.
  rownames(count_mat) = meta_names
  
  # Filter count matrix based on user defined minimum number of cells and samples.
  ind = count_mat >= min_cells
  meta_to_keep = apply(ind, 1, function(i) {
    sum(i) >= min_samples
  })
  count_mat = count_mat[meta_to_keep, , drop = FALSE]


  # Create factors data frame, which will be used to create our design matrix.
  factors = data.frame(row.names = sample_df[, 1], check.names = FALSE)
  for (i in 1:length(sample_df[, 1])) {
    idx = grep(paste(sample_df[i, 1], collapse = "|"), rownames(factors))
    for (a in vars_of_interest) {
      if (!(a %in% attributes(factors)$names)) {
        factors[[a]] = NA
      }
      factors[[a]][idx] = sample_df[i, a]
    }
  }
  
  # Determine the factors that should be used to create a "group" factor that
  # combines individual factors.
  comp_factors = c()
  for (i in 1:length(comparisons)) {
    new_atts = attributes(comparisons[[i]])$names
    new_atts = new_atts[!(new_atts %in% comp_factors)]
    comp_factors = append(comp_factors, new_atts)
  }
  
  # Create a new "group" column that concatenates the factors for comparisons.
  factors$group = do.call(paste, c(factors[comp_factors], sep = "_"))
  
  # Convert the columns in the "factors" data frame to factors.
  factors[] = lapply(factors, as.factor)
  
  # factors$NAC = relevel(factors$NAC, "No.NAC")
  
  # Create design matrix. 
  design_matrix = model.matrix(~ 0 + factors$group)
  colnames(design_matrix) = gsub("factors$group", "", colnames(design_matrix), fixed = TRUE)
  
  # Define comparisons by groups.
  grp_comps = list()
  i = 0
  for (comp in comparisons) {
    i = i+1
    factors_idx1 = rep(TRUE, nrow(factors))
    factors_idx2 = rep(TRUE, nrow(factors))
    for (a in attributes(comp)$names) {
      # If an attribute in comp has only one element, this implies it is the same
      # for both factor levels to be compared.
      if (length(comp[[a]])==1) {comp[[a]] = list(comp[[a]], comp[[a]])}
      # Get the logical row indices of the factors dataframe that correspond to the
      # groups to be compared.
      factors_idx1_a = rep(FALSE, nrow(factors))
      for (grp_var in comp[[a]][[1]]) {
        factors_idx1_a = factors_idx1_a | (factors[[a]]==grp_var)
      }
      factors_idx1 = factors_idx1 & factors_idx1_a
      factors_idx2_a = rep(FALSE, nrow(factors))
      for (grp_var in comp[[a]][[2]]) {
        factors_idx2_a = factors_idx2_a | (factors[[a]]==grp_var)
      }
      factors_idx2 = factors_idx2 & factors_idx2_a
    }
    # Store the names of groups for each comparison in a list of vectors.
    ### Without prepending group factor names by "group":
    grp1 = unique(factors$group[factors_idx1])
    grp2 = unique(factors$group[factors_idx2])
    ### With prepending group factor names by "group":
    # grp1 = paste0("group", unique(factors$group[factors_idx1]))
    # grp2 = paste0("group", unique(factors$group[factors_idx2]))
    if (!is.null(attributes(comparisons))) {
      att_name = attributes(comparisons)$names[i]
      if (is.character(att_name) && att_name != "" && !is.null(att_name)) {
        grp_comps[[att_name]] = list(grp1, grp2)
      } else {
        grp_comps[[i]] = list(grp1, grp2)
      }
    } else {
      grp_comps[[i]] = list(grp1, grp2)
    }
  }
  
  # Define comparison names to be used in the contrasts matrix.
  comparison_names = list()
  for(i in 1:length(grp_comps)){
    comp_name = gsub("[^[:alnum:].]", "_", paste(unlist(lapply(lapply(
      grp_comps[[i]], gsub, pattern="^group", replacement=""), paste, collapse="_OR_")), collapse = "__vs__"))
    if (!is.null(attributes(grp_comps))) {
      comparison_names[[attributes(grp_comps)$names[[i]]]] = comp_name
    } else {
      comparison_names[[i]] = comp_name
    }
  }
  
  # Create contrasts matrix.
  Contrasts = c()
  for (i in 1:length(grp_comps)) {
    grp1 = grp_comps[[i]][[1]]
    if (length(grp1) > 1) {
      grp1 = paste0("(", paste(grp1, collapse = "+"), ")/", length(grp1))
    }
    grp2 = grp_comps[[i]][[2]]
    if (length(grp2) > 1) {
      grp2 = paste0("(", paste(grp2, collapse = "+"), ")/", length(grp2))
    }
    Contrasts = append(Contrasts, paste(grp1, grp2, sep = " - "))
  }
  Contrasts = limma::makeContrasts(contrasts = Contrasts, levels = design_matrix)
  colnames(Contrasts) = comparison_names


  # Calculate normalization factors if desired by the user.
  # If there is an error, double check that norm_factors has the correct value.
  if (!is.numeric(norm_factors)) {
    if (normalize & norm_factors == "TMM") {
      norm_factors = calcNormFactors(count_mat, method = "TMM")
    }
  }
  
  # Create DGEList object.
  if (normalize) {
    dge_list = DGEList(counts = count_mat,
                       samples = sample_df,
                       group = factors$group,
                       norm.factors = norm_factors,
                       remove.zeros = TRUE)
  } else {
    dge_list = DGEList(counts = count_mat,
                       samples = sample_df,
                       group = factors$group,
                       remove.zeros = TRUE)
  }
  
  # Estimate dispersions. The option `trend_method = "none"` is used, as the 
  # dispersion-mean relationship of flow data typically does not resemble that of 
  # RNAseq data.
  disp = estimateDisp(dge_list, design = design_matrix, trend.method = "none")
  
  # Fit GLM models.
  glm_fit = glmFit(disp, design = design_matrix)
}

# Perform likelihood ratio tests, and save results as a .csv file for each
# comparison in the directory "Analysis Results/edgeR Analysis/Tables/".
# If you will be running an edgeR analysis on multiple FlowSOM objects, or will
# be running the script multiple times for any reason, you may wish to create
# sub-directories within "Tables/" to stay organized. To do this, you may
# simply set the "sub_directory" parameter directly below to whatever you would
# like this sub-directory to be called. You do not have to include the forward 
# slash in the name. It is recommended you name the directory after the FlowSOM
# object you are testing.
#
# Note: in edgeR user guide, only one contrast is tested at a time. The same is 
# done here.
sub_directory = "week19_lymph_B_norm"
if (!dir.exists(paste0(dir_edger_tables, sub_directory, "/"))) {
  dir.create(paste0(dir_edger_tables, sub_directory, "/"))
}
for (i in 1:ncol(Contrasts)) {
  glm_lrt = glmLRT(glm_fit, contrast = Contrasts[, i])
  top = topTags(glm_lrt, adjust.method = "BH", sort.by = "none")
  
  write.csv(top, file = paste0(dir_edger_tables, sub_directory, "/", 
                               names(comparisons)[i], "_", sub_directory, 
                               "_edger.csv"))
}

#########
# limma #
#########

# The code below uses the package limma to perform differential marker expression
# analysis on a FlowSOM object.
#
# Within each identified cell type in the FlowSOM object, for each marker 
# specified by the user, a linear model is fit, with each coefficient 
# corresponding to a group of interest (e.g. female MIBC). The models created 
# are linear mixed effects models, where sample ID is treated as a random effect, 
# and group is treated as a fixed effect. This type of model has been shown to 
# notably reduce the family-wise type I error rate in single cell data. For 
# each comparison of interest, the results are returned as a .csv file, where 
# tests are ranked by their adjusted p-values (corrected using the Benjamini & 
# Hochberg method) alongside the responsible channel and cell type.
# 
# Not all of the variables discussed below are necessary for every analysis.
# However, it is crucial that `markers_of_interest` is defined, and that before
# every run of the loop, `pval_dfs` and `df_full` are reassigned to empty data
# frames.

# Define the markers you are interested in checking for differential expression.
# Please make sure that you type the marker names the same as they appear in 
# the .fcs files.
# IMPORTANT: There should be NO cell markers in common with those you used for
# clustering. 
markers_of_interest = c("PHA-L")

# Define path to aggregate file used for this FlowSOM object, so that the transformed
# medians may be calculated.
agg = read.FCS("Data/Aggregates/aggregate_A.fcs")

# Define path to the transformation used for preprocessing the data, again so
# the transformed medians may be calculated.
trans = readRDS("RDS/logicle_transformation.rds")
inv = inverseLogicleTransform(trans)

####################################################
# VARIABLES FOR INCORPORATING FMO/ISOTYPE CONTROLS #
####################################################

# Do you have FMO or ISO files you would like to use in the analysis? If so,
# you must first preprocess these files as you would for any other regular
# data .fcs file. Then, for each marker for which these are applicable, define
# which files are relevant, in a list. Make sure the full filepath is included
# in the filename.
marker_list = list("PHA-L" = list.files("Data/Preprocessed FMO/",
                                        full.names = TRUE),
                   "IL10R" = list.files("Data/Preprocessed Iso/",
                                        full.names = TRUE))

# It is necessary to specify which FMO file each data file corresponds to. This
# information will be contained in a data frame, which is created in the same manner
# as the object `sample_df` above; you may either read in a .csv file, or create
# a data frame manually in R. You may also decide to add a column specifying 
# the appropriate FMO file for each filename. In that case, simply set fmo_info
# equal to `sample_df`.
fmo_info = sample_df

# Specify the name of the column containing your original file names.
filename_col = "File.Name"

# Specify the names of the columns containing your FMO/Iso file names. The order
# in which you list them should correspond to the order of the markers in
# `marker_list`.
fmo_col = c("FMO", "Iso")
names(fmo_col) = names(marker_list)

# Remove samples that were excluded from the analysis from the data frame, and 
# reorder its rows such that they are the same as the file order contained in 
# the `prepr_files` variable.
matched_ind = match(prepr_files, fmo_info[[filename_col]])
fmo_info = fmo_info[matched_ind, ]
rownames(fmo_info) = seq(1, length(prepr_files))

# Remove unwanted rows, if necessary.
if (nrow(fmo_info) > nrow(sample_df) && length(rows_to_remove) > 0) {
  fmo_info = fmo_info[-rows_to_remove, ]
}

##################################################################
# VARIABLES FOR USING FMO/ISOTYPE CONTROLS WITH RECLUSTERED DATA #
##################################################################

# If you are using FMO/Iso controls, but not with reclustered data, set the
# variables below equal to their empty values. 

# If the FlowSOM object you are testing is a reclustering of your initial
# FlowSOM object, please specify the names of the metaclusters that were used
# for this instance of reclustering. This is necessary to appropriately subset
# your FMO/Iso files.
subsetted_meta = c("CD8 T cells")

# You must specify the names of the directories containing your clustered files
# for each group of FMO/Isotype controls. The order in which you list them
# should correspond to the order in the `marker_list` and `fmo_col` objects. 
fmo_clustr_dir = c("Data/Clustered FMO Files/", "Data/Clustered Iso Files/")

# Specify a string that will be used to help name the subsetted FlowSOM objects.
# For example, if you are examining subsets of CD4 T cells, you may simply type
# "CD4 T cells", or "CD4".
str = "CD8 T cells"

##########################################
# FOR PREVIOUSLY CREATED FMO/ISO OBJECTS #
##########################################

# If you are using previously created FlowSOM objects, set this variable to
# `TRUE`. Otherwise set to `FALSE`.
prev_obj = FALSE

# If you would instead like to use previously created FlowSOM objects for your 
# FMO/Iso controls, you may specify their names here, and they will be 
# used for the analysis. Make sure that you list the names of the objects corresponding
# to the FlowSOM object defined at the start of the file - i.e. if you are analyzing
# monocytes, make sure that you are using the names of the FlowSOM objects composed
# only of monocytes, rather than the full clustering.
flow_names = c("FMO_fsom.rds", "Iso_fsom.rds")

# Specify aggregate filepaths as well.
agg_names = c(paste0(dir_agg, c("aggregate_FMO.fcs", "aggregate_Iso.fcs")))

# A helper to subset .fcs files.
filter_by_meta = function(ff, metaclusters) {
  filter = exprs(ff)[,"FlowSOM_meta"] %in% metaclusters
  filtered = ff[filter]
  return(filtered)
}

# Function for creating matrices of interest from `df_full` object.
get_sample_mat = function(sample_cluster_df, num_col, prepr_files, sample_df, 
                          sample_col, meta_names) {
  vec = sample_cluster_df[, num_col]
  
  mat = matrix(vec, ncol = length(vec)/length(prepr_files))
  
  df = as.data.frame(mat, row.names = sample_df[[sample_col]])
  
  colnames(df) = meta_names
  
  return(df)
}

# Test for differential marker expression between the specified groups of 
# interest for each metacluster.
if (TRUE) {
  # Set a seed for reproducibility.
  set.seed(42)
  
  # Create empty data frames for each comparison that will be made. These will
  # store the results of the analysis and be used to create the .csv files.
  pval_dfs = lapply(1:length(comparisons), function(i) {data.frame()})
  
  # Master data frame that will contain data about all metaclusters.
  # Maybe change to list of data frames, one for each marker?
  df_full = lapply(1:length(markers_of_interest), function(i) {data.frame()})
  
  for (m in markers_of_interest) {
    if (m %in% names(marker_list)) {
      print("Preparing controls...")
      
      ind = which(names(fmo_col) == m)
      
      col = fmo_col[ind]
      
      # change the order of the if statements
      if (prev_obj && exists("flow_names") && length(flow_names) > 0) {
        print("Using the FlowSOM objects defined in `flow_names` object...")
        agg_fmo = read.FCS(agg_names[ind])
        fmo_fsom = readRDS(paste0(dir_rds_edited, flow_names[ind]))
      } else if (exists("subsetted_meta") && length(subsetted_meta) > 0) {
        frames = lapply(list.files(fmo_clustr_dir[ind], full.names = TRUE), read.FCS)
        
        meta_name = readRDS(paste0(dir_rds_edited, "full_meta_names.rds"))
        meta_ind = which(meta_name %in% subsetted_meta)
        
        frames = lapply(frames, filter_by_meta, meta_ind)
        
        names(frames) = marker_list[[names(col)]]
        
        fs = as(frames, "flowSet")
        
        agg_fmo = AggregateFlowFrames(fs,
                                      cTotal = nrow(fSOM$data),
                                      writeOutput = TRUE,
                                      outputFile = paste0(dir_agg, "aggregate_", 
                                                          col, "_", str, ".fcs"))
        
        fmo_fsom = NewData(fSOM, agg_fmo)
        
        saveRDS(fmo_fsom, paste0(dir_rds_edited, col, "_", str, ".rds"))
      } else {
        # Create aggregate files from cells in FMO/Iso files.
        agg_fmo = AggregateFlowFrames(marker_list[[m]],
                                      cTotal = nrow(fSOM$data),
                                      writeOutput = TRUE,
                                      keepOrder = TRUE,
                                      outputFile = paste0(dir_agg, "aggregate_", col, ".fcs"))
        
        # Map above aggregate file to FlowSOM object of interest.
        fmo_fsom = NewData(fSOM, agg_fmo)
        
        saveRDS(fmo_fsom, paste0(dir_rds_edited, col, "_fsom.rds"))
        saveRDS(meta_names, paste0(dir_rds_edited, "full_meta_names.rds"))
        
        dir_save = paste0("Data/Clustered ", col, " Files/")
        
        if(!dir.exists(dir_save)) {
          dir.create(dir_save)
        }
        
        SaveClustersToFCS(fsom = fmo_fsom,
                          originalFiles = marker_list[[names(col)]], 
                          outputDir = dir_save)
      }
      # Convert the file name to numbers and store them as vectors. This is necessary
      # because FlowSOM objects represent each file as a number rather than a name.
      fmo_num = match(fmo_info[[col]], unique(fmo_info[[col]]))
    }
    for (k in 1:length(meta_names)) {
      # Get a subset of cells belonging to each sample.
      print(paste0("Analyzing ", meta_names[k], " ..."))
      
      if (length(rows_to_remove > 0)) {
        idx1 = GetMetaclusters(fSOM) == meta_names[k]
        idx2 = !(fSOM$data[, "File"] %in% rows_to_remove)
        inds = which(idx1 & idx2)
        fSOM_sub = FlowSOMSubset(fSOM, inds)
      } else {
        inds = which(GetMetaclusters(fSOM) == meta_names[k])
        fSOM_sub = FlowSOMSubset(fSOM, inds)
      }
      
      # Get metacluster-sample medians.
      sample_medians = data.frame(fSOM_sub$data,
                                  sample = fSOM_sub$data[, "File"],
                                  check.names = FALSE) %>%
        dplyr::group_by(sample, .drop = FALSE) %>%
        dplyr::summarise_all(stats::median) %>%
        dplyr::select(-sample) %>%
        data.frame(row.names = levels(sample),
                   check.names = FALSE)
      
      # Transform medians to linear scale.
      cols = which(colnames(sample_medians) %in% colnames(agg))
      
      trans_sample_medians = new("flowFrame", exprs = as.matrix(sample_medians[,cols]), 
                                 parameters = agg@parameters, 
                                 description = agg@description)
      
      trans_sample_medians = transform(trans_sample_medians, inv)
      trans_sample_medians = exprs(trans_sample_medians)
      
      # Subtract FMO or ISO if necessary
      if (m %in% names(marker_list)) {
        sub_inds = which(GetMetaclusters(fmo_fsom) == meta_names[k])
        fmo_fsom_sub = FlowSOMSubset(fmo_fsom, sub_inds)
        
        agg_fmo_sub = agg_fmo[sub_inds, ]
        agg_fmo_sub_trans = transform(agg_fmo_sub, inv)
        agg_fmo_sub_trans = exprs(agg_fmo_sub_trans)
        
        fmo_sample_medians = data.frame(agg_fmo_sub_trans,
                                        sample = agg_fmo_sub_trans[, "File"],
                                        check.names = FALSE) %>%
          dplyr::group_by(sample, .drop = FALSE) %>%
          dplyr::summarise_all(stats::median) %>%
          dplyr::select(-sample) %>%
          data.frame(row.names = levels(sample),
                     check.names = FALSE)

        chan = GetChannels(fSOM, names(marker_list[ind]))
        all_diff = c()
        for (i in 1:length(fmo_num)) {
          diff = fmo_sample_medians[fmo_num[i], chan]
          all_diff = c(all_diff, diff)
          trans_sample_medians[i, chan] = trans_sample_medians[i, chan] - diff
        }
      }
      
      
      # Checks for samples with too few cells. (figure out better way to do this)
      missing_samples = which(!(rownames(sample_df) %in% sample_medians$File))
      
      factors = data.frame(row.names = sample_df[, 1], check.names = FALSE)
      for (i in 1:length(sample_df[, 1])) {
        idx = grep(paste(sample_df[i, 1], collapse = "|"), rownames(factors))
        for (a in vars_of_interest) {
          if (!(a %in% attributes(factors)$names)) {
            factors[[a]] = NA
          }
          factors[[a]][idx] = sample_df[i, a]
        }
      }
      
      # Set data frame of factors.
      if (length(missing_samples) > 0) {
        factors = factors[-missing_samples, vars_of_interest]
      }
      
      # Create a new "group" column that concatenates the factors for comparisons.
      factors$group = do.call(paste, c(factors[comp_factors], sep = "_"))
      
      # Convert the columns in the "factors" data frame to factors.
      factors[] = lapply(factors, as.factor)
      
      # Create formula for design matrix.
      model_formula = "~ 0 + factors$group"
      model_formula = as.formula(model_formula)
      
      # Create design matrix.
      design_matrix = model.matrix(model_formula)
      colnames(design_matrix) = gsub("factors$group", "", colnames(design_matrix), fixed = TRUE)
      
      # Get data matrix made up of all cells that were randomly sampled.
      expr_matrix = t(rbind(sample_medians[, GetChannels(fSOM, m)]))
      expr_matrix = as.matrix(t(expr_matrix))
      rownames(expr_matrix) = m

      trans_matrix = t(rbind(trans_sample_medians[, GetChannels(fSOM, m)]))
      trans_matrix = as.matrix(t(trans_matrix))
      rownames(trans_matrix) = "Transformed Exp."
      
      if (length(missing_samples) > 0) {
        count_row = count_mat[k, -missing_samples]
      } else {
        count_row = count_mat[k, ]
      }
      
      # Make data frame to store info about this metacluster.
      df = data.frame(t(expr_matrix),
                      t(trans_matrix),
                      diff = all_diff,
                      counts = count_row,
                      cell_type = meta_names[k],
                      Group = factors$group,
                      check.names = FALSE)
      
      for (var in vars_of_interest[-1]) {
        df = data.frame(df, 
                        var = factors[[var]])
      }
      
      colnames(df)[7:(5+length(vars_of_interest))] = vars_of_interest[-1]
      
      # Append to greater data frame.
      df_ind = which(markers_of_interest == m)
      df_full[[df_ind]] = rbind(df_full[[df_ind]], df)
      
    }
    colnames(df_full[[df_ind]])[1] = markers_of_interest[df_ind]
  }
  # Create expression matrix for testing.
  expr_matrix = data.frame()
  for (i in 1:length(df_full)) {
    new_df = as.matrix(t(get_sample_mat(df_full[[i]], 1, prepr_files, sample_df, sample_col, meta_names)))
    expr_matrix = as.matrix(rbind(expr_matrix, new_df))
  }
  
  # Create linear models.
  # NOTE: weights questionable
  lm_model = lmFit(object = expr_matrix, 
                   design = design_matrix)
  
  # Perform statistical tests.
  contrasts_fit = contrasts.fit(lm_model, Contrasts)
  limma_ebayes = eBayes(contrasts_fit, trend = TRUE)
  
  # Create tables containing results of our statistical tests, and add them
  # to the data frame corresponding to the relevant comparison.
  for (i in 1:ncol(Contrasts)) {
    table = limma::topTable(limma_ebayes,
                            colnames(Contrasts)[i],
                            sort.by = "p",
                            number = 12,
                            p.value = 0.20)
    
    if (nrow(table) != 0) {
      table$marker = rep(1, nrow(table))
      
      for (j in 1:nrow(table)) {
        curr = rownames(table)[j]
        ix = ceiling(as.integer(curr)/length(meta_names))
        table$marker[j] = markers_of_interest[ix]
      }
      
      temp_ind = nrow(pval_dfs[[i]]) + 1
      pval_dfs[[i]] = rbind(pval_dfs[[i]], table)
    }
  }
}

# Write the data frames that were created during the analysis to .csv files,
# which are saved to the directory "Analysis Results/limma Analysis/Tables/".
# If you will be running a limma analysis on multiple FlowSOM objects, or will
# be running the script multiple times for any reason, you may wish to create
# sub-directories within "Tables/" to stay organized. To do this, you may
# simply set the "sub_directory" parameter directly below to whatever you would
# like this sub-directory to be called. You do not have to include the forward 
# slash in the name.
sub_directory = "human_cd8_fmo_and_iso_B_pseudobulk"
if (!dir.exists(paste0(dir_limma_tables, sub_directory, "/"))) {
  dir.create(paste0(dir_limma_tables, sub_directory, "/"))
}
pval_dfs = lapply(1:length(pval_dfs), function(i) {
  pval_df = pval_dfs[[i]]
  if (nrow(pval_df) != 0) {
    pval_df = pval_df[order(pval_df$adj.P.Val), ]
  }
  write.csv(pval_df, file = paste0(dir_limma_tables, sub_directory, "/", 
                                   names(comparisons)[i], "_", sub_directory, 
                                   "_limma.csv"))
  return(pval_df)
})

# Function to plot abundances by group for each metacluster. The following is
# an example of how it might be called:
# 
# abundance_boxplot(groups_of_interest = c("Female_BBN_N", 
#                                          "Female_BBN_Y", 
#                                           "Male_BBN_N", 
#                                           "Male_BBN_Y"), 
#                   facet_by = "Sex", 
#                   sample_cluster_df = df_full)
#
# The parameters are as follows:
# 
# groups_of_interest: A vector defining the groups whose data you would like
#   to use for plotting. The names of these groups must match those found in
#   the "Group" column of the `df_full` object created during the limma analysis.
#   To quickly check the names of your groups, try calling `levels(df_full$Group)`.
# 
# facet_by: A string specifying which variable you would like to facet by.
#   This variable must be contained within the `vars_of_interest` vector defined 
#   at the beginning of the file (excluding the name for sample ID).
# 
# sample_cluster_df: This is a data frame specifying information about all
#   metacluster-sample pairs. It is created during the for loop that performs the
#   analysis, and should be one of the data frames in the list `df_full`.
abundance_boxplot = function(groups_of_interest, facet_by, sample_cluster_df) {
  ind = which(sample_cluster_df$Group %in% groups_of_interest)
  
  sample_cluster_df = sample_cluster_df[ind, ]
  
  p = ggplot(sample_cluster_df, aes(x = sample_cluster_df[[facet_by]],
                                    y = 100*(count_row/sum(count_row)),
                                    fill = sample_cluster_df[[facet_by]])) +
    geom_boxplot() +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(title = paste("Relative Abundances in Cell Subtypes by", facet_by),
         x = "",
         y = "Proportion (%)",
         fill = facet_by) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  print(p)
}

abundance_boxplot(c("female_MIBC_NAC", "female_MIBC_No.NAC", "male_MIBC_NAC", "male_MIBC_No.NAC"), 
                  "Sex", df_full[[2]])

# Function to plot expression by group. The following is an example of how it 
# might be called:
# 
# expression_boxplot(groups_of_interest = c("Female_BBN_N", 
#                                           "Female_BBN_Y", 
#                                           "Male_BBN_N", 
#                                           "Male_BBN_Y"), 
#                    marker_of_interest = "PHA-L"
#                    group_by = "Sex", 
#                    sample_cluster_df = df_full)
#
# The parameters are as follows:
# 
# groups_of_interest: A vector defining the groups whose data you would like
#   to use for plotting. The names of these groups must match those found in
#   the "Group" column of the `df_full` object created during the limma analysis.
#   To quickly check the names of your groups, try calling `levels(df_full$Group)`.
# 
# marker_of_interest: A string specifying the marker you would like to plot
#   expression for. It should appear in the variable `markers_of_interest`, 
#   defined above at the beginning of the limma section.
# 
# group_by: A string specifying which variable you would like to group plots by.
#   This variable must be contained within the `vars_of_interest` vector defined 
#   at the beginning of the file (excluding the name for sample ID).
# 
# sample_cluster_df: This is a data frame specifying information about all
#   metacluster-sample pairs. It is created during the for loop that performs the
#   analysis, and should be one of the data frames in the list `df_full`.
expression_boxplot = function(groups_of_interest, marker_of_interest, group_by,
                              sample_cluster_df) {
  ind = which(sample_cluster_df$Group %in% groups_of_interest)
  
  sample_cluster_df = sample_cluster_df[ind, ]
  
  p = ggplot(sample_cluster_df, aes(x = cell_type,
                                    y = sample_cluster_df[[marker_of_interest]],
                                    fill = sample_cluster_df[[group_by]])) +
    geom_boxplot() +
    labs(title = paste(marker_of_interest, "in Cell Subtypes by", group_by),
         x = "Cell Type",
         y = "Median Expression",
         fill = group_by) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  print(p)
}

expression_boxplot(c("female_MIBC_NAC", "female_MIBC_No.NAC", 
                     "male_MIBC_NAC", "male_MIBC_No.NAC"),
                   "IL10R", "Sex", df_full[[2]])

expression_boxplot(c("female_MIBC_NAC", "female_MIBC_No.NAC",
                     "male_MIBC_NAC", "male_MIBC_No.NAC"),
                   "PHA-L", "NAC", df_full[[1]])


# edit_fsom = function(fsom, col, diff_mat) {
#   for (i in 1:nrow(diff_mat)) {
#     ind1 = fsom$data[, "File"] == i
#     for (j in 1:ncol(diff_mat)) {
#       ind2 = GetMetaclusters(fsom) == colnames(diff_mat)[j]
#       inds = which(ind1 & ind2)
#       print(length(inds))
#       fsom$data[inds, col] = fsom$data[inds, col] - diff_mat[i, j]
#     }
#   }
#   return(fsom)
# }

# Edit to take a certain subset rather than edit the whole file.
edit_agg = function(fsom, agg, col, diff_mat) {
  diff_mat = get_sample_mat(df_full[[1]], 3, prepr_files, sample_df, sample_col, meta_names)
  
  for (i in 1:nrow(diff_mat)) {
    ind1 = exprs(agg)[, "File"] == i
    for (j in 1:ncol(diff_mat)) {
      ind2 = GetMetaclusters(fsom) == colnames(diff_mat)[j]
      inds = which(ind1 & ind2)
      
      exprs(agg)[inds, col] = exprs(agg)[inds, col] - diff_mat[i, j]
    }
    print(paste0("Analyzing ", rownames(diff_mat)[i], "..."))
  }
  return(agg)
}

diff_mat = get_sample_mat(df_full[[1]], 3, prepr_files, sample_df, sample_col, meta_names)
df_trans = get_sample_mat(df_full[[1]], 2, prepr_files, sample_df, sample_col, meta_names)
write.csv(df_trans, "Analysis Results/full_26_merged_pha_dMFI_linear_table.csv")

rand_inds = sample(seq(1, nrow(agg)), 1000000)
agg_inv = agg
agg_inv = transform(agg_inv, inv)
exprs(agg_inv) = exprs(agg_inv)[rand_inds, ]
new_fsom = FlowSOMSubset(fSOM, rand_inds)
new_agg = edit_agg(new_fsom, agg_inv, "FITC-A", diff_mat)
agg_log = transform(new_agg, trans)
new_fsom$data = exprs(agg_log)

##############
# More Plots #
##############

# Helper function for calculating standard error.
se_fun = function(x) {
  sd = sd(x)
  se = sd/sqrt(length(x))
  return(se)
}

# Function to plot bar plot of MFIs or dMFIs. You must pass a .csv file (created
# by the `get_sample_mat()` function defined above), the `factors` object created
# above, and a list defining which groups you'd like to plot. It should look
# like the following: 
# 
# grps_of_interest = list("Male Ctrl" = c("male_Ctrl_X"),
#                         "Female Ctrl" = c("female_Ctrl_X"),
#                         "Male MIBC" = c("male_MIBC_No.NAC", "male_MIBC_NAC"),
#                         "Female MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC"))
# 
# The names and order of the groups should appear in the same way you would like
# them to appear in the bar plot. Make sure that the group names you supply are
# the same as they appear in the `factors` object. Finally, if necessary, you
# may adjust the y-axis upper limit with the `upper_lim` parameter.
bar_plot_fun = function(csv_file, factors, grps_of_interest, upper_lim = NULL) {
  orig_df = read.csv(csv_file, check.names = FALSE)
  orig_df = cbind(orig_df, factor_group = factors$group)
  
  df = orig_df[,-1]
  
  for (col in colnames(df)) {
    if (col != "factor_group") {
      df[, col] = as.numeric(df[, col])
    }
  }
  
  # Ensure "group" column is treated as a factor.
  df = df %>%
    mutate(factor_group = as.factor(factor_group))  
  
  df$Group = factor(seq(1:nrow(df)), levels = names(grps_of_interest))
  for (i in 1:nrow(df)) {
    group = as.character(df$factor_group[i])
    for (j in 1:length(grps_of_interest)) {
      if (group %in% grps_of_interest[[j]]) {
        df$Group[i] = names(grps_of_interest)[j]
      }
    }
  }
  
  df = df[, -(which(colnames(df) == "factor_group"))]
  point_df = pivot_longer(df, cols = -Group)
  
  median_df = df %>%
    group_by(Group) %>%
    summarise_if(is.numeric, median, na.rm = TRUE) %>%
    ungroup()
  
  se_df = df %>%
    group_by(Group) %>%
    summarise_if(is.numeric, se_fun) %>%
    ungroup()
  
  med_dat_long = pivot_longer(median_df, cols = -Group)
  se_dat_long = pivot_longer(se_df, cols = -Group)
  
  if (is.null(upper_lim)) {
    upper_lim = max(med_dat_long$value)
    upper_lim = upper_lim + (upper_lim/2)
  }
  
  ggplot(med_dat_long, aes(x = name, y = value, fill = Group)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black") + 
    geom_point(data = point_df,
               mapping = aes(name, value, shape = Group, fill = Group),
               position = position_jitterdodge(jitter.width = 0.05, dodge.width = 1, seed = 42),
               size = 1,
               inherit.aes = FALSE) +
    geom_errorbar(mapping = aes(x = name, ymin = value-se_dat_long$value, ymax = value+se_dat_long$value),
                  position = position_dodge(width = 0.9, preserve = "single"),
                  width = 0.5,
                  inherit.aes = TRUE) +
    scale_fill_viridis(discrete = TRUE, option = "G", begin = 0.2) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8), aspect.ratio = 0.6) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, upper_lim))
}


#########
# UMAPs #
#########

# Create UMAP colored by metacluster. (how is it different from `plot_umap()`?)
plot_group_umap = function(fsom, clustered_markers, factors, num_cells, seed = 42) {
  set.seed(seed)
  
  groups = levels(factors$group)
  
  all_ind = c()
  
  for (group in groups) {
    grp_samples = which(factors[, "group"] == group)
    inds = which(fsom$data[, "File"] %in% grp_samples)
    samp = sample(inds, num_cells)
    all_ind = c(all_ind, samp)
  }
  
  dat = fsom$data[all_ind, GetChannels(fsom, clustered_markers)]
  
  umap = umap(dat)
  
  meta_vec = GetMetaclusters(fsom)[all_ind]
  # meta_vec = factor(meta_vec, levels = c("CD4 T cells", "CD8 T cells", "gdT cells",
  #                                        "B cells", "NK cells", "NK T cells", "Monocytes",
  #                                        "cDC", "pDC", "Undefined"))
  umap_df = data.frame(umap$layout, Metacluster = meta_vec, Indices = all_ind)
  
  print(ggplot(umap_df) + scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = Metacluster), 
                                                        pointsize = 2) + 
          theme_void())
  
  return(umap_df)
}