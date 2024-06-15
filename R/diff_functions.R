#' prepareSampleInfo
#'
#' Prepare sample information for DE and DA analysis via a .csv file.
#'
#' @param filepath A filepath to a .csv file containing sample info.
#' @param name_col The name of the column containing the experiment's sample names.
#' @param filename_col The name of the column containing the experiment's .fcs file names.
#' @param cols_to_use A character vector specifying the column names other than
#' sample name and filename that are relevant for the experiment.
#' @param samples_to_remove An integer vector defining the rows of any samples
#' that should be excluded from the analysis, either because they were discarded
#' during preprocessing or are not of interest. Default is none (\code{NULL}).
#'
#' @details
#' If you have a .csv file specifying sample information, it can be read in and
#' prepared by this function. It MUST have columns for file names and sample
#' names, and any additional columns may specify parameters of interest,
#' like NAC vs. No NAC. Each column should be named after what it represents.
#' This file will be read in as a data frame and will be used to construct the
#' design matrix used for analysis.
#'
#' *REMOVE* differences from orig. script: all columns other than filename and sample names
#' have make.names() applied. 
#'
#' @return A data frame containing sample information.
#'
#' @export
prepareSampleInfo <- function(filepath, name_col, filename_col, cols_to_use,
                              samples_to_remove = NULL) {
  # Read in .csv file.
  sample_df <- utils::read.csv(file = filepath,
                               header = TRUE)
  
  # Get sample names to remove, if necessary
  if (is.numeric(samples_to_remove)) {
    samples_to_remove <- sample_df[[name_col]][samples_to_remove]
  }
  

  # add check that filenames in .csv file exist in the preprocessed dir?

  # Check for typos.
  input_cols <- c(name_col, filename_col, cols_to_use)
  invalid_cols <- which(!(input_cols %in% colnames(sample_df)))
  if (length(invalid_cols) > 0) {
    stop(paste(input_cols[invalid_cols], collapse = ", "), " are not valid column names. ",
         "Please check your input for typos. The following column names are valid: ",
         paste(colnames(sample_df), collapse = ", "))
  }

  # Move sample name and filename columns to the first two positions of the data frame.
  sample_df <- sample_df[, c(name_col,
                             filename_col,
                             setdiff(names(sample_df), c(name_col, filename_col)))]

  # Edit column names and row entries to be R friendly, excluding sample and filenames.
  sample_df[, -c(1,2)] <- sapply(sample_df[, -c(1,2)], make.names)
  
  # Remove samples that were excluded from the analysis from the data frame 
  # `sample_df`, and reorder its rows such that they are the same as 
  # the file order in `dir_prepr()`.
  prepr_files <- list.files(path = dir_prepr())
  matched_ind <- match(prepr_files, sample_df[[filename_col]])
  sample_df <- sample_df[matched_ind, ]
  rownames(sample_df) <- seq(1, length(prepr_files))
  
  # Remove any samples that are not of interest.
  if (!is.null(samples_to_remove)) {
    sample_df <- sample_df[-(which(sample_df[[name_col]] %in% samples_to_remove)), ]
  }

  # Reorder rows alphabetically
  # sample_df <- sample_df[order(sample_df[[filename_col]]), ] # relevant when dir_prepr() not used

  return(sample_df)
}

#' makeFactorDF
#'
#' ! consider merging with prepareSampleInfo()
#'
#' @param sample_df A data frame containing sample information, generated either
#' manually or by [prepareSampleInfo()].
#' @param comparisons A named list of named lists, defining the groups to be
#' compared during analysis. See example for how this variable should be defined.
#'
#' @return A data frame of factors, with an added column for group.
#'
#' @export
makeFactorDF <- function (sample_df, comparisons) {
  # Create factors data frame, which will be used to create our design matrix.
  factors = data.frame(row.names = sample_df[, 1], check.names = FALSE)
  for (i in 1:length(sample_df[, 1])) {
    idx = grep(paste(sample_df[i, 1], collapse = "|"), rownames(factors))
    for (a in colnames(sample_df)[-2]) { # originally `vars_of_interest`
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

  # add choice to relevel?

  return(factors)
}

#' makeDesignMatrix
#'
#' Generate a design matrix. (edit example)
#'
#' @param factors A data frame, generated either manually or by [makeFactorDF()].
#' @param comparisons A named list of named lists, defining the groups to be
#' compared during analysis. See example for how this variable should be defined.
#'
#' @return A matrix where each column is a group of interest, and each row is a
#' sample.
#'
#' @export
#'
#' @examples
#' samples <- prepareSampleInfo("filepath", "Sample.Name", "File.Name", c("Sex", "Disease"))
#'
#' comparisons <- list(
#'   male_vs_female = list(Sex = list("male", "female")),
#'   male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female"))
#'   )
#'
#' factors <- makeFactorDF(sample_df, comparisons)
#'
#' design_matrix <- makeDesignMatrix(factors, comparisons)
makeDesignMatrix <- function(factors, comparisons) {
  design_matrix <- stats::model.matrix(~ 0 + factors$group)
  colnames(design_matrix) <- gsub("factors$group", "", colnames(design_matrix), fixed = TRUE)

  return(design_matrix)
}

#' makeContrastsMatrix
#'
#' Generate a contrasts matrix. (edit example)
#'
#' @inheritParams makeDesignMatrix
#'
#' @return A matrix, where each column corresponds to a comparison, and each row
#' corresponds to a group.
#'
#' @export
#'
#' @examples
#' samples <- prepareSampleInfo("filepath", "Sample.Name", "File.Name", c("Sex", "Disease"))
#'
#' comparisons <- list(
#'   male_vs_female = list(Sex = list("male", "female")),
#'   male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female"))
#'   )
#'
#' factors <- makeFactorDF(sample_df, comparisons)
#'
#' contrasts <- makeContrastsMatrix(factors, comparisons)
makeContrastsMatrix <- function(factors, comparisons) {

  # Create design matrix.
  design_matrix <- makeDesignMatrix(factors, comparisons)

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

  return(Contrasts)
}

#' makeCountMatrix
#'
#' Generate matrix of sample/metacluster cell counts.
#'
#' @param fsom A FlowSOM object.
#' @param factors A data frame containing sample information, generated by
#' [makeFactorDF()].
#' @param meta_names A vector of metacluster names of interest. By default, all
#' metacluster names are used.
#' @param min_cells An integer, the minimum number of cells a metacluster should
#' have in a specified number of samples to be included in the analysis.
#' @param min_samples An integer, the minimum number of samples a metacluster
#' should have at least \code{min_cells} events in to be included in the analysis.
#' By default, this is half the total number of samples.
#'
#' @return A matrix, where each column represents a sample, and each row
#' represents a metacluster.
#'
#' @export
makeCountMatrix <- function(fsom, factors, meta_names = levels(fsom$metaclustering),
                            min_cells = 3, min_samples = NULL) {
  if (nrow(factors) > 0) {

    # Set `min_samples` if it is NULL
    if (is.null(min_samples)) {
      min_samples <- nrow(factors)/2
    }

    count_mat <- c()

    # For each sample
    for (i in 1:nrow(factors)) {
      # Get all cells belonging to the current sample
      ind <- which(fsom$data[, "File"] == i)

      # Get metacluster assignments for current cells
      meta_assignments <- FlowSOM::GetMetaclusters(fsom)[ind]

      meta_counts <- c()

      # For each metacluster
      for (j in meta_names) {
        count <- length(which(meta_assignments == j))
        meta_counts <- c(meta_counts, count)
      }

      # Bind metacluster counts for current sample to matrix
      count_mat <- cbind(count_mat, meta_counts)
      num_col <- ncol(count_mat)

    }

    # Rename columns and rows.
    colnames(count_mat) <- rownames(factors)
    rownames(count_mat) <- meta_names

    # Filter count matrix based on given minimum number of cells and samples
    ind <- count_mat >= min_cells
    meta_to_keep <- apply(ind, 1, function(i) {
      sum(i) >= min_samples
    })
    count_mat <- count_mat[meta_to_keep, , drop = FALSE]
  }

  return(count_mat)
}

#' doDAAnalysis
#'
#' Do differential abundance analysis on the data.
#'
#' @inheritParams makeContrastsMatrix
#' @param design_matrix A design matrix as generated by [makeDesignMatrix()].
#' @param count_mat A count matrix as generated by [makeCountMatrix()].
#' @param contrasts A contrasts matrix as generated by [makeContrastsMatrix()].
#' @param norm_method The normalization method to be used and passed to
#' [edgeR::calcNormFactors()], if any.
#' @param dir_tables The name of a subdirectory to save the resulting tables in.
#' Default is \code{NULL}, in which case the files will be saved in \code{dir_edger}.
#'
#' @details
#' This function performs tests for differential abundances between groups of
#' interest defined by a contrasts matrix. Models are fit and compared via
#' likelihood ratio tests using functions from edgeR. Results are saved as a
#' .csv file for each comparison in the directory defined by \code{dir_edger()} at
#' runtime. If running an edgeR analysis on multiple FlowSOM objects, or running the
#' script multiple times for any reason, the user may wish to create
#' sub-directories within \code{dir_edger} to stay organized. To do this,
#' simply set the \code{dir_tables} parameter to a preferred sub-directory name.
#'
#' @seealso [edgeR::calcNormFactors()], [edgeR::DGEList()], [edgeR::estimateDisp()],
#' [edgeR::glmFit()], [edgeR::glmLRT()], [edgeR::topTags()].
#'
#' @export
doDAAnalysis <- function(design_matrix, count_mat, contrasts, factors,
                         norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                         dir_tables = NULL) {

  # Get sample info from factors object
  sample_info <- factors[, -c(1, ncol(factors))] # remove?

  # Create DGEList object
  if (!is.null(norm_method)) { # if user wants to normalize data
    norm_factors <- edgeR::calcNormFactors(count_mat, method = norm_method)
    dge_list <- edgeR::DGEList(counts = count_mat,
                               samples = sample_info,
                               group = factors$group,
                               norm.factors = norm_factors,
                               remove.zeros = TRUE)
  } else {
    dge_list <- edgeR::DGEList(counts = count_mat,
                               samples = sample_info,
                               group = factors$group,
                               remove.zeros = TRUE)
  }

  # Estimate dispersions. The option `trend_method = "none"` is used, as the
  # dispersion-mean relationship of flow data typically does not resemble that of
  # RNAseq data.
  disp <- edgeR::estimateDisp(dge_list, design = design_matrix, trend.method = "none")

  # Fit GLM models
  glm_fit <- edgeR::glmFit(disp, design = design_matrix)

  if (!is.null(dir_tables)) { # save to user specified subdirectory
    dir_write <- file.path(dir_edger(), dir_tables)
  } else { # save directly to edgeR directory
    dir_write <- dir_edger()
  }

  # Check if directory exists and create if not
  if (!dir.exists(dir_write)) {
    dir.create(dir_write, recursive = TRUE)
  }

  # Perform likelihood ratio tests, and save results as a .csv file for each
  # comparison in the directory `dir_edger`.
  #
  # Note: in edgeR user guide, only one contrast is tested at a time. The same is
  # done here.
  for (i in 1:ncol(contrasts)) {
    glm_lrt <- edgeR::glmLRT(glm_fit, contrast = contrasts[, i])
    top <- edgeR::topTags(glm_lrt, adjust.method = "BH", sort.by = "none")

    utils::write.csv(top, file = file.path(dir_write,
                                           paste0(colnames(contrasts)[i],
                                                  "_", dir_tables, "_edger.csv")))
  }
}



#' clusterControls
#'
#' Prepare isotype or FMO controls for calculating delta MFI.
#'
#' @param sample_df A data frame of sample info, as generated by [prepareSampleInfo()].
#' @param ctrl_col The name of the column containing control filenames.
#' @param dir_prepr_ctrl The directory containing preprocessed control files.
#' @param marker The marker that this control corresponds to.
#' @param fsom A FlowSOM object, to be used for differential expression analysis.
#' @param subsetted_meta A character vector specifying the metaclusters that were
#' used for metaclusters, if \code{fsom} is a reclustering. Default is \code{NULL}.
#' @param dir_clustr_ctrl A filepath to the directory containing clustered control
#' files, if \code{fsom} is a reclustering. Default is \code{NULL}.
#' @param parent_ctrl_fsom The parent control FlowSOM whose cells are to be
#' reclustered, if \code{fsom} is a reclustering. Default is \code{NULL}.
#'
#' @details
#' Additional details...
#' *reduce parameters if possible, edit documentation*
#'
#' @return A FlowSOM object.
#'
#' @export
clusterControls <- function(sample_df, ctrl_col, dir_prepr_ctrl, marker, fsom,
                            subsetted_meta = NULL, dir_clustr_ctrl = NULL,
                            parent_ctrl_fsom = NULL
                            ) { # parent_ctrl_fsom unnecessary if you add metadata

  if (!is.null(subsetted_meta)) {
    # frames <- lapply(list.files(dir_clustr_ctrl, full.names = TRUE), flowCore::read.FCS)
    files <- file.path(dir_clustr_ctrl, gsub(".fcs", "_FlowSOM.fcs", sample_df[[ctrl_col]]))
    frames <- lapply(files, flowCore::read.FCS)

    # Get indices of metaclusters used for the reclustering
    meta_names <- levels(parent_ctrl_fsom$metaclustering)
    meta_ind <- which(meta_names %in% subsetted_meta)

    # Create flowSet of cells only belonging to given metaclusters
    frames <- lapply(frames, filterFCS, metaclusters = meta_ind)
    names(frames) <- file.path(dir_prepr_ctrl, sample_df[[ctrl_col]]) # sample_df unnecessary?
    fs <- methods::as(frames, "flowSet")

    # Make new aggregate
    agg_ctrl <- FlowSOM::AggregateFlowFrames(fs,
                                            cTotal = nrow(fsom$data),
                                            writeOutput = TRUE,
                                            outputFile = file.path(dir_agg(),
                                                                   paste0("aggregate_", ctrl_col, ".fcs"))
                                            )

    # Map aggregate data to FlowSOM object of interest
    ctrl_fsom <- FlowSOM::NewData(fsom, agg_ctrl)
    saveRDS(ctrl_fsom, paste0(dir_rds_edited(), "fsom_", ctrl_col, ".rds"))

    # Not a reclustering, and no pre-existing clustering given
  } else {
    # Create aggregate files from cells in control files.
    ctrl_files <- file.path(dir_prepr_ctrl, sample_df[[ctrl_col]])
    agg_ctrl <- FlowSOM::AggregateFlowFrames(ctrl_files,
                                            cTotal = nrow(fsom$data),
                                            writeOutput = TRUE,
                                            keepOrder = TRUE,
                                            outputFile = file.path(dir_agg(),
                                                                   paste0("aggregate_", ctrl_col, ".fcs"))
                                            )

    # Map above aggregate file to FlowSOM object of interest.
    ctrl_fsom <- FlowSOM::NewData(fsom, agg_ctrl)

    # Save control FlowSOM object
    saveRDS(ctrl_fsom, paste0(dir_rds_edited(), ctrl_col, "_fsom2.rds"))

    # Save clustered control files
    dir_save <- file.path("Data", paste0("Clustered ", ctrl_col, " Files2"))
    if (!dir.exists(dir_save)) {
      dir.create(dir_save, recursive = TRUE)
    }
    FlowSOM::SaveClustersToFCS(fsom = ctrl_fsom,
                               originalFiles = list.files(dir_clustr_ctrl, full.names = TRUE),
                               outputDir = dir_save)
  }

  return(ctrl_fsom)
}

#' getSampleMFIs
#'
#' @param input A FlowSOM object, flowFrame, or matrix of expression values. 
#' Must have a \code{"File"} column in the data matrix.
#'
#' @return A data frame, where each row is a sample, and each column is a channel/marker.
#' 
#' @export
getSampleMFIs <- function(input) {
  if (methods::is(input, "FlowSOM")) {
    input <- input$data
  } else if (methods::is(input, "flowFrame")) {
    input <- flowCore::exprs(input)
  }
  
  medians <- data.frame(input,
                        sample = input[, "File"],
                        check.names = FALSE) %>%
    dplyr::group_by(sample, .drop = FALSE) %>%
    dplyr::summarise_all(stats::median) %>%
    dplyr::select(-sample) %>%
    data.frame(row.names = levels(sample),
               check.names = FALSE)
}

#' prepareControlInfo
#'
#' @param markers ...
#' @param ctrl_cols ...
#' @param ctrl_prepr_dirs ...
#' @param ctrl_clustr_dirs ...
#' @param parent_ctrl_fsom ...
#'
#' @return A data frame with FMO/Isotype control information.
#' 
#' @export
prepareControlInfo <- function(markers, ctrl_cols, ctrl_prepr_dirs, 
                               ctrl_clustr_dirs = NULL, parent_ctrl_fsom = NULL) {
  ctrl_df <- data.frame("Column.Name" = ctrl_cols,
                        "Prepr.Dir" = ctrl_prepr_dirs,
                        row.names = markers)
  
  if (!is.null(ctrl_clustr_dirs)) {
      ctrl_df <- data.frame(ctrl_df, 
                            "Clustered.Parent.Dir" = ctrl_clustr_dirs,
                            "Parent.fsom" = parent_ctrl_fsom)
  }
  
  return(ctrl_df)
}


doDEAnalysis <- function(fsom, sample_df, design_matrix, contrasts, count_mat, 
                         markers_of_interest, meta_names, prepr_transform, 
                         controls_df = NULL, ctrl_fsom_names = NULL,
                         subsetted_meta = NULL) {
  # Set a seed for reproducibility.
  set.seed(42)
  
  de_res <- list("tests" = list(),
                 "data" = list())
  
  # Create empty data frames for each comparison that will be made. These will
  # store the results of the analysis and be used to create the .csv files.
  pval_dfs <- lapply(1:ncol(contrasts), function(i) {data.frame()})
  
  # Master data frame that will contain data about all metaclusters.
  # Maybe change to list of data frames, one for each marker?
  df_full <- lapply(1:length(markers_of_interest), function(i) {data.frame()})
  
  # Samples excluded from analysis
  removed_samples <- which(!(levels(factor(fsom$data[, "File"])) %in% rownames(sample_df)))
  
  for (m in markers_of_interest) {
    if (!is.null(controls_df) && m %in% rownames(controls_df)) {
      print("Preparing controls...")
      
      row <- which(rownames(controls_df) == m)
      
      # Pre-existing clustering given
      if (!is.null(ctrl_fsom_names)) {
        print("Using the FlowSOM objects defined in `ctrl_fsom_names` object...")
        
        ctrl_fsom <- readRDS(paste0(dir_rds_edited(), ctrl_fsom_names[row]))
      } else {
        ctrl_fsom <- clusterControls(sample_df, 
                                     ctrl_col = controls_df[row, 1], 
                                     dir_prepr_ctrl = controls_df[row, 2], 
                                     marker = m, 
                                     fsom = fsom,
                                     subsetted_meta = subsetted_meta,
                                     dir_clustr_ctrl = controls_df[row, 3],
                                     parent_ctrl_fsom = controls_df[row, 4]) ###
      }
      agg_ctrl <- flowCore::read.FCS(file.path(dir_agg(), names(ctrl_fsom$metaData)[[1]])) # !!! check NewData adds metaData
    }
    for (k in 1:length(meta_names)) {
      print(paste0("Analyzing ", meta_names[k], " ..."))
      
      # Get cells belonging to each sample in current metacluster
      if (length(removed_samples) > 0) { 
        idx1 <- FlowSOM::GetMetaclusters(fsom) == meta_names[k]                
        idx2 <- !(fsom$data[, "File"] %in% removed_samples)                      
        inds <- which(idx1 & idx2)
        fsom_sub <- FlowSOM::FlowSOMSubset(fsom, inds)
      } else {
        inds <- which(FlowSOM::GetMetaclusters(fsom) == meta_names[k])
        fsom_sub <- FlowSOM::FlowSOMSubset(fsom, inds)
      }
      
      # Get cluster-sample medians.
      sample_medians <- getSampleMFIs(fsom_sub)
      
      # Get transformation used for preprocessing and its inverse
      trans <- readRDS(file.path(dir_rds_edited(), prepr_transform))
      inv <- flowCore::inverseLogicleTransform(trans)
      
      # Get aggregate file for currect FlowSOM object
      agg <- flowCore::read.FCS(file.path(dir_agg(), names(fsom$metaData)[[1]]))
      # Only columns used for clustering
      cols <- which(colnames(sample_medians) %in% colnames(agg))
      
      # Transform medians to linear scale
      trans_sample_medians <- methods::new("flowFrame", exprs = as.matrix(sample_medians[, cols]), 
                                           parameters = agg@parameters, 
                                           description = agg@description)
      
      trans_sample_medians <- transform(trans_sample_medians, inv)
      trans_sample_medians <- flowCore::exprs(trans_sample_medians)
      
      # Subtract FMO or ISO if necessary
      if (m %in% names(controls_df)) {
        sub_inds <- which(FlowSOM::GetMetaclusters(ctrl_fsom) == meta_names[k]) # ensure that excluded samples
        ctrl_fsom_sub <- FlowSOM::FlowSOMSubset(ctrl_fsom, sub_inds)            # aren't analyzed
        
        agg_ctrl_sub <- agg_ctrl[sub_inds, ] # same as ctrl_fsom$data?
        agg_ctrl_sub_trans <- transform(agg_ctrl_sub, inv)
        
        fmo_sample_medians <- getSampleMFIs(agg_ctrl_sub_trans)
        
        chan <- FlowSOM::GetChannels(fsom, rownames(controls_df)[row])
        all_diff <- c()
        for (i in 1:nrow(sample_df)) { # nrow(sample_df) was length(fmo_num)
          diff <- fmo_sample_medians[rownames(sample_df)[i], chan]
          all_diff <- c(all_diff, diff)
          trans_sample_medians[i, chan] <- trans_sample_medians[i, chan] - diff
        }
      }
      
      # Checks for samples with too few cells
      # Note: samples not of interest were removed from local variables earlier, 
      # before finding `sample_medians`
      missing_samples <- which(!(rownames(sample_df) %in% sample_medians$File))
      
      # Get vector of medians for current marker
      expr_matrix <- t(rbind(sample_medians[, FlowSOM::GetChannels(fsom, m)]))
      expr_matrix <- as.matrix(t(expr_matrix))
      rownames(expr_matrix) <- m
      
      # Get vector of transformed medians for current marker
      trans_matrix <- t(rbind(trans_sample_medians[, FlowSOM::GetChannels(fsom, m)]))
      trans_matrix <- as.matrix(t(trans_matrix))
      rownames(trans_matrix) <- "Transformed Exp."
      
      # Remove samples with too few cells from count matrix, if necessary
      # Note: `count_mat` is found via `factors`, so samples not of interest are already removed 
      if (length(missing_samples) > 0) {
        count_row <- count_mat[k, -missing_samples] # does the same need to be done for design?
      } else {                                      # what if a sample has enough cells in one meta/sample
        count_row <- count_mat[k, ]                 # pair for one marker but not the other?
      }
      
      # Make data frame to store info about this metacluster.
      df <- data.frame(t(expr_matrix),
                       t(trans_matrix),
                       diff = all_diff,
                       counts = count_row,
                       cell_type = meta_names[k],
                       Group = factors$group, #######
                       check.names = FALSE)
      
      # Append to greater data frame.
      df_ind <- which(markers_of_interest == m)
      df_full[[df_ind]] <- rbind(df_full[[df_ind]], df)
      
    }
    colnames(df_full[[df_ind]])[1] <- markers_of_interest[df_ind]
  }
  
  # Create expression matrix for testing.
  expr_matrix <- data.frame()
  for (i in 1:length(df_full)) {
    new_df <- as.matrix(t(getSampleMetaMatrix(df_full[[i]], 1, sample_df, meta_names)))
    expr_matrix <- as.matrix(rbind(expr_matrix, new_df))
  }
  
  # Create linear models.
  # NOTE: weights questionable
  lm_model <- limma::lmFit(object = expr_matrix, 
                           design = design_matrix)
  
  # Perform statistical tests.
  contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
  limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)
  
  # Create tables containing results of our statistical tests, and add them
  # to the data frame corresponding to the relevant comparison.
  for (i in 1:ncol(contrasts)) {
    table <- limma::topTable(limma_ebayes,
                             colnames(contrasts)[i],
                             sort.by = "p",
                             number = 12,
                             p.value = 0.20)
    
    if (nrow(table) != 0) {
      table$marker <- rep(1, nrow(table))
      
      for (j in 1:nrow(table)) {
        curr <- rownames(table)[j]
        ix <- ceiling(as.integer(curr)/length(meta_names))
        table$marker[j] <- markers_of_interest[ix]
      }
      
      temp_ind <- nrow(pval_dfs[[i]]) + 1
      pval_dfs[[i]] <- rbind(pval_dfs[[i]], table)
    }
  }
  
  de_res$tests <- pval_dfs
  de_res$data <- df_full
  
  return(de_res)
}

# Function for creating matrices of interest from `df_full` object.
getSampleMetaMatrix = function(df_full, col_to_use, sample_df, meta_names) {
  vec <- df_full[, col_to_use]
  mat <- matrix(vec, ncol = length(vec)/nrow(sample_df))
  
  df <- as.data.frame(mat, row.names = sample_df[[1]])
  colnames(df) <- meta_names # levels(df[["cell_type"]]) ?
  
  return(df)
}
