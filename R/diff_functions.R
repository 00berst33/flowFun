#' prepareSampleInfo
#'
#' Prepare sample information for DE and DA analysis via a .csv file.
#'
#' @param filepath A filepath to a .csv file containing sample info.
#' @param name_col The name of the column containing the experiment's sample names.
#' @param filename_col The name of the column containing the experiment's .fcs file names.
#' @param comparisons A list of comparisons, where each comparison is further
#' defined by a list of factor levels. Each comparison and group of factor
#' levels should be named. See examples for further detail on how this parameter
#' should be defined.
#' @param samples_to_remove An integer vector defining the rows of any samples
#' that should be excluded from the analysis, either because they were discarded
#' during preprocessing or are not of interest. Default is none (\code{NULL}).
#'
#' @details
#' If you have a .csv file specifying sample information, it can be read in and
#' prepared by this function. It MUST have columns for file names and sample
#' names, and any additional columns should specify parameters of interest,
#' like NAC vs. No NAC. Each column should be named after what it represents.
#' This file will be read in as a data frame and will be used for comparative
#' analysis.
#'
#' *REMOVE* differences from orig. script: all columns other than filename and sample names
#' have make.names() applied.
#' *note* need to make edits to better incorporate controls
#'
#' @return A data frame containing sample information.
#'
#' @export
#'
#' @examples
#' # It is important that `comparisons` is defined properly.
#'
#' # The first level should be named lists, defining each comparison of interest.
#' # The lists within each of these comparisons should include at least two levels
#' # of a group factor to be compared. To compare two groups within another
#' # group (i.e. NAC vs. No NAC in females only), specify the level(s)
#' # of interest for both of these factors. For example:
#' comparisons <- list(
#'   male_vs_female = list(Sex = list("male", "female")),
#'   male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female"))
#' )
#'
#' # Here, at the first level of the object, there are two comparisons, named
#' # "male_vs_female", and "male_vs_female_mibc".
#'
#' # They are further defined by lists of factor levels;
#' # `list(Sex = list("male", "female"))` and
#' # `list(Disease = "MIBC", Sex = list("male", "female"))`, respectively.
#'
#' # So, the first element defines the comparison male vs. female,
#' # and the second element defines the comparison MIBC male vs. MIBC female.
#'
#' # Defining a comparison such as `list(Sex = "male", Disease = "MIBC")`,
#' # where only one level is given for each factor, is useless.
#'
#' # Every element in the list defining a comparison must be named, and correspond
#' # to a column in the file found at `filepath`. So, in this example, the given
#' # file must have a column named "Sex", with entries "male" or "female", and
#' # a column named "Disease", with entries "MIBC".
#'
#' # Finally, note that if a column has many different possible entries - for example,
#' # if the "Disease" column had entries "MIBC", "NMIBC", and "Ctrl" - it is only
#' # necessary to list the groups of interest. So in this case, the above example
#' # is still appropriate.
#'
#' # Read in .csv file
#' filepath <- system.file("extdata", "sample_information.csv", package = "flowFun")
#' sample_df <- prepareSampleInfo(filepath,
#'                                name_col = "Sample_ID",
#'                                filename_col = "Filenames",
#'                                comparisons = comparisons)
#'
#' # Note the added "group" column
#' print(head(sample_df))
#'
#' # Read in the same file, this time removing a couple samples
#' sample_df.rm <- prepareSampleInfo(filepath,
#'                                   name_col = "Sample_ID",
#'                                   filename_col = "Filenames",
#'                                   comparisons = comparisons,
#'                                   samples_to_remove = c(2, 4))
#'
#' # Row names do not change when samples are removed
#' print(head(sample_df.rm))
prepareSampleInfo <- function(filepath, name_col, filename_col, comparisons,
                              samples_to_remove = NULL) {
  # Read in .csv file.
  sample_df <- utils::read.csv(file = filepath,
                               header = TRUE)

  # Get sample names to remove, if necessary
  if (is.numeric(samples_to_remove)) {
    samples_to_remove <- sample_df[[name_col]][samples_to_remove]
  }

  # add check that filenames in .csv file exist in the preprocessed dir?

  # Get names of columns relevant for comparative analysis
  comp_factors <- c()
  for (i in 1:length(comparisons)) {
    new_atts <- attributes(comparisons[[i]])$names
    new_atts <- new_atts[!(new_atts %in% comp_factors)]
    comp_factors <- append(comp_factors, new_atts)
  }

  # Check for typos.
  input_cols <- c(name_col, filename_col, comp_factors)
  invalid_cols <- which(!(input_cols %in% colnames(sample_df)))
  if (length(invalid_cols) > 0) {
    stop(paste(input_cols[invalid_cols], collapse = ", "), " are not valid column names. ",
         "Please check your input for typos. The following column names are valid: ",
         paste(colnames(sample_df), collapse = ", "))
  }

  # Move sample name, filename, and `comp_factors` columns to the front of the data frame
  sample_df <- sample_df[, c(name_col,
                             filename_col,
                             comp_factors,
                             setdiff(names(sample_df), c(name_col, filename_col, comp_factors)))]

  # Edit column names and row entries to be R friendly, excluding sample and filenames.
  sample_df[, -c(1,2)] <- sapply(sample_df[, -c(1,2)], make.names)

  # Remove samples that were excluded from the analysis from the data frame
  # `sample_df`, and reorder its rows such that they are the same as
  # the file order in `dir_prepr()`.
  prepr_files <- sample_df[[filename_col]][-c(30, 40)]#list.files(path = dir_prepr())
  matched_ind <- match(prepr_files, sample_df[[filename_col]])
  sample_df <- sample_df[matched_ind, ]
  rownames(sample_df) <- seq(1, length(prepr_files))

  # Remove any samples that are not of interest.
  if (!is.null(samples_to_remove)) {
    sample_df <- sample_df[-(which(sample_df[[name_col]] %in% samples_to_remove)), ]
  }

  # Reorder rows alphabetically
  # sample_df <- sample_df[order(sample_df[[filename_col]]), ] # relevant when dir_prepr() not used

  # Create a new "group" column that concatenates the factors for comparisons.
  sample_df$group <- do.call(paste, c(sample_df[comp_factors], sep = "_"))

  # Convert the columns in the "factors" data frame to factors.
  sample_df[] <- lapply(sample_df, as.factor)

  # add choice to relevel?

  return(sample_df)
}

#' prepareControlInfo
#'
#' use in clusterControls ?
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

#' makeFactorDF
#'
#' !!! original function now mostly merged with prepareSampleInfo(), only really
#' necessary if sampleinfo is given as a list
#'
#' @param sample_info A list of lists containing sample information
#' @param comparisons A named list of named lists, defining the groups to be
#' compared during analysis. See example for how this variable should be defined.
#'
#' @return A data frame of factors, with an added column for group.
#'
#' @export
makeFactorDF <- function (sample_info, comparisons) {
  # Create factors data frame, which will be used to create our design matrix.
  sample_names <- c()
  for (i in 1:length(sample_info)) {
    sample_names <- c(sample_names, sample_info[[i]][[1]])
  }
  factors <- data.frame(row.names = sample_names, check.names = FALSE) # consider using file names read from dir_prepr() instead
  for (i in 1:length(sample_info)) {
    idx <- grep(paste(sample_info[[i]][[1]], collapse = "|"), rownames(factors))
    for (a in attributes(sample_info[[i]])$names[-1]) {
      if (!(a %in% attributes(factors)$names)) {
        factors[[a]] <- NA
      }
      factors[[a]][idx] <- sample_info[[i]][[a]]
    }
  }

  # Determine the factors that should be used to create a "group" factor that
  # combines individual factors.
  comp_factors <- c()
  for (i in 1:length(comparisons)) {
    new_atts <- attributes(comparisons[[i]])$names
    new_atts <- new_atts[!(new_atts %in% comp_factors)]
    comp_factors <- append(comp_factors, new_atts)
  }

  # Create a new "group" column that concatenates the factors for comparisons.
  factors$group <- do.call(paste, c(factors[comp_factors], sep = "_"))

  # Convert the columns in the "factors" data frame to factors.
  factors[] <- lapply(factors, as.factor)

  # add choice to relevel?

  return(factors)
}

#' makeDesignMatrix
#'
#' Generate a design matrix. (edit example)
#'
#' @param sample_df A data frame, generated either manually or by [prepareSampleInfo()].
#'
#' @return A matrix where each column is a group of interest, and each row is a
#' sample.
#'
#' @export
#'
#' @examples
#' filepath <- system.file("extdata", "sample_information.csv", package = "flowFun")
#' samples <- prepareSampleInfo(filepath, "Sample.Name", "File.Name", comparisons)
#'
#' design <- makeDesignMatrix(samples)
makeDesignMatrix <- function(sample_df) {
  design <- stats::model.matrix(~ 0 + sample_df$group)
  colnames(design) <- gsub("sample_df$group", "", colnames(design), fixed = TRUE)
  rownames(design) <- rownames(sample_df)

  return(design)
}

#' makeContrastsMatrix
#'
#' Generate a contrasts matrix. (edit example)
#'
#' @inheritParams makeDesignMatrix
#' @inheritParams prepareSampleInfo
#'
#' @return A matrix, where each column corresponds to a comparison, and each row
#' corresponds to a group.
#'
#' @export
#'
#' @examples
#' comparisons <- list(
#'   male_vs_female = list(Sex = list("male", "female")),
#'   male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female"))
#' )
#'
#' filepath <- system.file("extdata", "sample_information.csv", package = "flowFun")
#' samples <- prepareSampleInfo(filepath, "Sample.Name", "File.Name", comparisons)
#'
#' contrasts <- makeContrastsMatrix(samples, comparisons)
makeContrastsMatrix <- function(sample_df, comparisons) {

  # Create design matrix.
  design <- makeDesignMatrix(sample_df)

  # Define comparisons by groups.
  grp_comps = list()
  i = 0
  for (comp in comparisons) {
    i = i+1
    factors_idx1 = rep(TRUE, nrow(sample_df))
    factors_idx2 = rep(TRUE, nrow(sample_df))
    for (a in attributes(comp)$names) {
      # If an attribute in comp has only one element, this implies it is the same
      # for both factor levels to be compared.
      if (length(comp[[a]])==1) {comp[[a]] = list(comp[[a]], comp[[a]])}
      # Get the logical row indices of the sample_df dataframe that correspond to the
      # groups to be compared.
      factors_idx1_a = rep(FALSE, nrow(sample_df))
      for (grp_var in comp[[a]][[1]]) {
        factors_idx1_a = factors_idx1_a | (sample_df[[a]]==grp_var)
      }
      factors_idx1 = factors_idx1 & factors_idx1_a
      factors_idx2_a = rep(FALSE, nrow(sample_df))
      for (grp_var in comp[[a]][[2]]) {
        factors_idx2_a = factors_idx2_a | (sample_df[[a]]==grp_var)
      }
      factors_idx2 = factors_idx2 & factors_idx2_a
    }
    # Store the names of groups for each comparison in a list of vectors.
    ### Without prepending group factor names by "group":
    grp1 = unique(sample_df$group[factors_idx1])
    grp2 = unique(sample_df$group[factors_idx2])
    ### With prepending group factor names by "group":
    # grp1 = paste0("group", unique(sample_df$group[factors_idx1]))
    # grp2 = paste0("group", unique(sample_df$group[factors_idx2]))
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
  Contrasts = limma::makeContrasts(contrasts = Contrasts, levels = design)
  colnames(Contrasts) = comparison_names

  return(Contrasts)
}

#' makeCountMatrix
#'
#' Generate matrix of sample/metacluster cell counts.
#'
#' @param input A FlowSOM object or data frame.
#' @param meta_names A vector of metacluster names of interest. By default, all
#' metacluster names are used.
#' @param min_cells An integer, the minimum number of cells a metacluster should
#' have in a specified number of samples to be included in the analysis.
#' @param min_samples An integer, the minimum number of samples a metacluster
#' should have at least \code{min_cells} events in to be included in the analysis.
#' By default, this is half the total number of samples.
#' @param ... other arguments:
#'
#' \code{sample_df}: If \code{input} is a FlowSOM object, a data frame from
#' [prepareSampleInfo()].
#'
#' @return A matrix, where each column represents a sample, and each row
#' represents a metacluster.
#'
#' @seealso [makeCountMatrix.FlowSOM()]
#'
#' @export
makeCountMatrix <- function(input, meta_names = NULL,
                            min_cells = 3, min_samples = NULL, ...) {
  result <- UseMethod("makeCountMatrix")
  return(result)
}

#' makeCountMatrix.FlowSOM
#'
#' Generate matrix of sample/metacluster cell counts from FlowSOM object.
#'
#' @param sample_df A data frame containing sample information, generated by
#' [prepareSampleInfo()].
#'
#' @keywords internal
#' @export
makeCountMatrix.FlowSOM <- function(input, sample_df, meta_names = NULL,
                                    min_cells = 3, min_samples = NULL) {
  if (nrow(sample_df) > 0) {
    # Set `meta_names` if it is NULL
    if (is.null(meta_names)) {
      meta_names <- levels(input$metaclustering)
    }

    # Set `min_samples` if it is NULL
    if (is.null(min_samples)) {
      min_samples <- nrow(sample_df)/2
    }

    counts <- c()

    # For each sample
    for (i in 1:nrow(sample_df)) {
      # Get all cells belonging to the current sample
      ind <- which(input$data[, "File"] == i)

      # Get metacluster assignments for current cells
      meta_assignments <- FlowSOM::GetMetaclusters(input)[ind]

      meta_counts <- c()

      # For each metacluster
      for (j in meta_names) {
        count <- length(which(meta_assignments == j))
        meta_counts <- c(meta_counts, count)
      }

      # Bind metacluster counts for current sample to matrix
      counts <- cbind(counts, meta_counts)
      num_col <- ncol(counts)

    }

    # Rename columns and rows.
    colnames(counts) <- sample_df[, 1]
    rownames(counts) <- meta_names

    # Filter count matrix based on given minimum number of cells and samples
    ind <- counts >= min_cells
    meta_to_keep <- apply(ind, 1, function(i) {
      sum(i) >= min_samples
    })
    counts <- counts[meta_to_keep, , drop = FALSE]
  }

  return(counts)
}

#' makeCountMatrix.data.frame
#'
#' Generate matrix of sample/metacluster cell counts.
#'
#' @keywords internal
#' @export
makeCountMatrix.data.frame <- function(input, meta_names = NULL,
                                       min_cells = 3, min_samples = NULL) {
  .id <- Metacluster <- n <- cell_count <- NULL
  # Set `min_samples` if it is NULL
  if (is.null(min_samples)) {
    min_samples <- input %>%
      tidytable::pull(.id) %>%
      unique() %>%
      length()
  }

  # Set `meta_names` if necessary
  if (is.null(meta_names)) {
    meta_names <- input %>%
      tidytable::pull(Metacluster) %>%
      unique()
  }

  # Get sample-metacluster count matrix
  counts <- input %>%
    tidytable::filter(Metacluster %in% meta_names) %>%
    tidytable::summarise(cell_count = n(),
                         .by = c(.id, Metacluster),
                         .groups = "drop") %>%
    tidytable::pivot_wider(names_from = .id,
                           values_from = cell_count,
                           values_fill = list(cell_count = 0)) %>%
    as.matrix(rownames = TRUE)

  # Rename columns
  colnames(counts) <- basename(colnames(counts))

  # Filter count matrix based on given minimum number of cells and samples
  ind <- counts >= min_cells
  meta_to_keep <- apply(ind, 1, function(i) {
    sum(i) >= min_samples
  })
  counts <- counts[meta_to_keep, , drop = FALSE]

  rem <- names(meta_to_keep)[which(!meta_to_keep)]
  if (length(rem) > 0) {
    print(paste("Metacluster(s)", paste(rem, collapse = ", "), "were discarded due to too few cells."))
  }

  return(counts)
}

#' doDAAnalysis
#'
#' Do differential abundance analysis on the data.
#'
#' @inheritParams makeContrastsMatrix
#' @param design A design matrix as generated by [makeDesignMatrix()].
#' @param counts A count matrix as generated by [makeCountMatrix()].
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
doDAAnalysis <- function(design, counts, contrasts, sample_df,
                         norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                         dir_tables = NULL) {

  # Get sample info from sample_df object
  sample_info <- sample_df[, -c(1, ncol(sample_df))] # remove?

  # Create DGEList object
  if (!is.null(norm_method)) { # if user wants to normalize data
    norm_factors <- edgeR::calcNormFactors(counts, method = norm_method)
    dge_list <- edgeR::DGEList(counts = counts,
                               samples = sample_info,
                               group = sample_df$group,
                               norm.factors = norm_factors,
                               remove.zeros = TRUE)
  } else {
    dge_list <- edgeR::DGEList(counts = counts,
                               samples = sample_info,
                               group = sample_df$group,
                               remove.zeros = TRUE)
  }

  # Estimate dispersions. The option `trend_method = "none"` is used, as the
  # dispersion-mean relationship of flow data typically does not resemble that of
  # RNAseq data.
  disp <- edgeR::estimateDisp(dge_list, design = design, trend.method = "none")

  # Fit GLM models
  glm_fit <- edgeR::glmFit(disp, design = design)

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
    top <- edgeR::topTags(glm_lrt, adjust.method = "BH", sort.by = "PValue")

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


#' clusterControlsTable
#'
#' Cluster FMO and/or isotype controls to facilitate calculation of delta MFIs.
#'
#' @param input A data.table, as generated by [getTableFromFCS()] or [doPreprocessing()].
#' @param fsom The original FlowSOM object, which has not been reclustered.
#' @param num_cells The number of cells to use in the resulting clustering.
#' @param subsetted_fsom If the FlowSOM object of interest is a reclustering,
#' the reclustered FlowSOM object.
#' @param subsetted_meta If the FlowSOM object of interest is a reclustering,
#' the name of the metacluster that was used.
#' @param save_fsom Optional, an .rds filename to be used to name the resulting
#' FlowSOM object. Default is \code{NULL}, which does not save the object.
#'
#' @return A data table with metacluster and cluster assignments for each cell.
#'
#' @export
clusterControlsTable <- function(input, fsom, subsetted_fsom = NULL, subsetted_meta = NULL,
                                 num_cells = nrow(fsom$data), save_fsom = NULL) {
  Metacluster <- Cluster <- cell_id <- NULL

  # Convert input to flowFrame for use with FlowSOM
  mat <- as.matrix(input, rownames = TRUE)
  ctrl_ff <- flowCore::flowFrame(mat)

  # Map input to FlowSOM object of interest.
  ctrl_fsom <- FlowSOM::NewData(fsom, ctrl_ff)

  # Get metacluster and cluster labels
  meta_labels <- FlowSOM::GetMetaclusters(ctrl_fsom)
  clust_labels <- factor(FlowSOM::GetClusters(ctrl_fsom))

  # Append labels as columns to the data.table
  ctrl_dt <- input %>%
    tidytable::mutate(Cluster = clust_labels,
                      Metacluster = meta_labels,
                      .keep = "all")

  # If the FlowSOM object of interest is a reclustering
  if (!is.null(subsetted_meta)) {
    # Get only cells belonging to subsetted metacluster
    ctrl_dt <- ctrl_dt %>%
      tidytable::filter(Metacluster == subsetted_meta) %>%
      tidytable::select(-c(Cluster, Metacluster))

    # Sample cells and subset data
    print(ifelse(num_cells > nrow(ctrl_dt), nrow(ctrl_dt), num_cells))
    print(nrow(ctrl_dt))
    inds <- sample(nrow(ctrl_dt), ifelse(num_cells > nrow(ctrl_dt), nrow(ctrl_dt), num_cells))
    ctrl_dt <- ctrl_dt[inds, ]

    # Create flowFrame of cells only belonging to given metaclusters
    ctrl_ff <- flowCore::flowFrame(ctrl_dt)

    # Map aggregate data to FlowSOM object of interest
    ctrl_fsom <- FlowSOM::NewData(subsetted_fsom, ctrl_ff)

    # Get metacluster and cluster labels
    meta_labels <- FlowSOM::GetMetaclusters(ctrl_fsom)
    clust_labels <- factor(FlowSOM::GetClusters(ctrl_fsom))

    ctrl_dt <- ctrl_dt %>%
      tidytable::mutate(Cluster = clust_labels,
                        Metacluster = meta_labels,
                        .keep = "all")

    # Not a reclustering, and no pre-existing clustering given
  } else {
    # Sample cells and subset data
    inds <- sample(nrow(ctrl_dt), ifelse(num_cells > nrow(ctrl_dt), nrow(ctrl_dt), num_cells))
    ctrl_dt <- ctrl_dt[inds, ]
  }

  # Save clustered object if desired
  if (!is.null(save_fsom)) {
    saveRDS(ctrl_fsom, paste0(dir_rds_edited(), save_fsom)) ###
  }

  # Sort by cell ID
  ctrl_dt <- ctrl_dt %>%
    tidytable::arrange(cell_id)

  return(ctrl_dt)
}

#' getSampleMFIs
#'
#' A helper function to get sample MFIs.
#' add functionality for tidytable
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

# method for tidytable
getSampleMFIsTable <- function(input) {
  .id <- NULL

  fcs_cols <- attr(input, "cols_from_fcs")

  medians <- input %>%
    tidytable::group_by(.id) %>%
    tidytable::summarise(tidytable::across(.cols = fcs_cols,
                                           .fns = stats::median,
                                           .drop = "keep")) %>%
    tidytable::select(c(.id, fcs_cols)) %>%
    data.frame(row.names = ".id",
               check.names = FALSE)

  rownames(medians) <- basename(rownames(medians))

  return(medians)
}

#' getSampleMetaMatrix
#'
#' Helper function for creating matrices of interest from \code{df_full}.
#'
#' @param df_full A data frame with sample/metacluster information for the current
#' marker.
#' @param col_to_use Which column of the data frame should be used to create the
#' matrix.
#'
#' @keywords internal
#'
#' @return A matrix, where each column is a metacluster and each row is a sample.
#' The entries may be median expression on either a linear or logicle-transformed
#' scale.
#'
#' @export
getSampleMetaMatrix = function(df_full, col_to_use) {
  cell_type <- NULL

  if (is.numeric(col_to_use)) {
    col_to_use <- colnames(df_full)[col_to_use]
  }

  df <- df_full %>%
    dplyr::select(c(col_to_use, sample, cell_type)) %>%
    tidyr::pivot_wider(values_from = col_to_use, names_from = cell_type) %>%
    data.frame(check.names = FALSE) %>%
    dplyr::select(-sample)

  rownames(df) <- unique(df_full$sample)

  return(df)
}

#' doDEAnalysis
#'
#' Perform differential expression analysis on markers of interest.
#'
#' @param input A FlowSOM object or data frame whose cells have been clustered.
#' @param sample_df A data frame of sample info, as generated by [prepareSampleInfo()].
#' @param design A design matrix, as generated by [makeDesignMatrix()].
#' @param contrasts A contrasts matrix, as generated by [makeContrastsMatrix()].
#' @param counts A count matrix, as generated by [makeContrastsMatrix()].
#' @param markers_of_interest A character vector of markers to test for differential
#' expression.
#' @param meta_names A character vector of metacluster names to test. Default is
#' all metaclusters.
#' @param ... other arguments:
#'
#' - prepr_transform: The transformation applied during preprocessing.
#' - controls_df: A data frame as generated by [prepareControlInfo()],
#' containing information about FMO or isotype controls to be used to calculate
#' delta MFIs.
#' - ctrl_fsom_names: Optional, a character string defining names of FlowSOM
#' .rds files previously generated from [clusterControls()].
#' - subsetted_meta: If \code{fsom} is a reclustering, and \code{controls_df}
#' is not \code{NULL}, a string or list of strings defining the metacluster(s)
#' \code{fsom} is a reclustering of.
#'
#' @details
#' Additional details...
#'
#' !!! counts parameter may be useless, Group column in `df_full` no longer needed
#'
#' @return A list of lists containing test results and data matrices.
#'
#' @export
doDEAnalysis <- function(input, sample_df, design, contrasts, counts,
                         markers_of_interest, meta_names = NULL, ...) {
  results <- UseMethod("doDEAnalysis")
  return(results)
}

#' doDEAnalysis.FlowSOM
#'
#' Perform differential expression analysis on markers of interest.
#'
#' @param prepr_transform A filepath to the .rds file for the logicle
#' transformation applied during preprocessing. Default is
#' \code{file.path("RDS", "logicle_transformation.rds")}.
#' @param controls_df Optional, a data frame as generated by [prepareControlInfo()],
#' containing information about FMO or isotype controls to be used to calculate
#' delta MFIs.
#' @param ctrl_fsom_names Optional, a character string defining names of FlowSOM
#' .rds files previously generated from [clusterControls()].
#' @param subsetted_meta If \code{fsom} is a reclustering, and \code{controls_df}
#' is not \code{NULL}, a string or list of strings defining the metacluster(s)
#' \code{fsom} is a reclustering of.
#'
#' @export
doDEAnalysis.FlowSOM <- function(input, sample_df, design, contrasts, counts,
                                 markers_of_interest, meta_names = NULL,
                                 prepr_transform = NULL, controls_df = NULL,
                                 ctrl_fsom_names = NULL, subsetted_meta = NULL) {
  # Set a seed for reproducibility
  set.seed(42)

  # Get metacluster names if none were given
  if (is.null(meta_names)) {
    meta_names <- levels(input$metaclustering)
  }

  # Determine samples excluded from analysis, if any
  removed_samples <- which(!(levels(factor(input$data[, "File"])) %in% rownames(sample_df)))

  # Create empty data frames for each comparison that will be made. These will
  # store the results of the analysis and be used to create the .csv files.
  pval_dfs <- lapply(1:ncol(contrasts), function(i) {data.frame()})

  # Master data frame that will contain data about all metaclusters.
  df_full <- lapply(1:length(markers_of_interest), function(i) {data.frame()})

  # Initialize object that will be returned
  de_res <- list("tests" = list(), # `pval_dfs`
                 "data" = list()) # `df_full`

  # Prepare data for each marker for testing
  for (m in markers_of_interest) {
    print(paste0("Analyzing ", m, " ..."))
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
                                     fsom = input,
                                     subsetted_meta = subsetted_meta,
                                     dir_clustr_ctrl = controls_df[row, 3],
                                     parent_ctrl_fsom = controls_df[row, 4]) ###
      }
      agg_ctrl <- flowCore::read.FCS(file.path(dir_agg(), names(ctrl_fsom$metaData)[[1]])) # !!! check NewData adds metaData
    }
    for (k in 1:length(meta_names)) {
      print(paste0("Preparing ", meta_names[k], " ..."))

      # Get cells belonging to each sample in current metacluster
      if (length(removed_samples) > 0) {
        idx1 <- FlowSOM::GetMetaclusters(input) == meta_names[k]
        idx2 <- !(input$data[, "File"] %in% removed_samples)
        inds <- which(idx1 & idx2)
        fsom_sub <- FlowSOM::FlowSOMSubset(input, inds)
      } else {
        inds <- which(FlowSOM::GetMetaclusters(input) == meta_names[k])
        fsom_sub <- FlowSOM::FlowSOMSubset(input, inds)
      }

      # Get cluster-sample medians.
      sample_medians <- getSampleMFIs(fsom_sub)

      # Get transformation used for preprocessing and its inverse
      if (is.null(prepr_transform)) {
        prepr_transform <- file.path("RDS", "logicle_transformation.rds")
      }

      trans <- readRDS(prepr_transform)
      inv <- flowCore::inverseLogicleTransform(trans)

      # Get aggregate file for currect FlowSOM object
      agg <- flowCore::read.FCS(file.path(dir_agg(), names(input$metaData)[[1]]))
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

        chan <- FlowSOM::GetChannels(input, rownames(controls_df)[row])
        # all_diff <- c()
        for (i in 1:nrow(sample_df)) { # nrow(sample_df) was length(fmo_num)
          diff <- fmo_sample_medians[rownames(sample_df)[i], chan]
          # all_diff <- c(all_diff, diff)
          trans_sample_medians[i, chan] <- trans_sample_medians[i, chan] - diff
        }
      }

      # Checks for samples with too few cells
      # Note: samples not of interest were removed from local variables earlier,
      # before finding `sample_medians`
      missing_samples <- which(!(rownames(sample_df) %in% sample_medians$File))

      # Get vector of medians for current marker
      expr_matrix <- t(rbind(sample_medians[, FlowSOM::GetChannels(input, m)]))
      expr_matrix <- as.matrix(t(expr_matrix))
      rownames(expr_matrix) <- m

      # Get vector of transformed medians for current marker
      trans_matrix <- t(rbind(trans_sample_medians[, FlowSOM::GetChannels(, m)]))
      trans_matrix <- as.matrix(t(trans_matrix))
      rownames(trans_matrix) <- "Transformed Exp."

      # Remove samples with too few cells from count matrix, if necessary
      # Note: `counts` is found via `sample_df`, so samples not of interest are already removed
      # if (length(missing_samples) > 0) {
      #   count_row <- counts[k, -missing_samples] # does the same need to be done for design?
      # } else {                                      # what if a sample has enough cells in one meta/sample
      #   count_row <- counts[k, ]                 # pair for one marker but not the other?
      # }

      # Make data frame to store info about this metacluster.
      df <- data.frame(t(expr_matrix),
                       t(trans_matrix),
                       # counts = count_row,
                       sample = sample_df[rownames(sample_medians), 1],
                       cell_type = meta_names[k],
                       # Group = sample_df$group,
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
    new_df <- as.matrix(t(getSampleMetaMatrix(df_full[[i]], 1)))
    expr_matrix <- as.matrix(rbind(expr_matrix, new_df))
  }

  # Create linear models.
  # NOTE: weights questionable
  lm_model <- limma::lmFit(object = expr_matrix,
                           design = design)

  # Perform statistical tests.
  contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
  limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)

  # Create tables containing results of our statistical tests, and add them
  # to the data frame corresponding to the relevant comparison.
  for (i in 1:ncol(contrasts)) {
    table <- limma::topTable(limma_ebayes,
                             coef = colnames(contrasts)[i],
                             sort.by = "p",
                             p.value = 0.20)

    # Add column for marker ID
    if (nrow(table) != 0) {
      table$marker <- rep(markers_of_interest[1], nrow(table))

      if (length(markers_of_interest) > 1) { # if more than one marker is being tested
        for (j in 1:nrow(table)) {
          curr <- rownames(table)[j]
          ix <- ceiling(as.integer(curr)/length(meta_names))
          table$marker[j] <- markers_of_interest[ix]
        }
        rownames(table) <- NULL
      }

      # Append table of results to returned list
      temp_ind <- nrow(pval_dfs[[i]]) + 1
      pval_dfs[[i]] <- rbind(pval_dfs[[i]], table)
    }
  }

  de_res$tests <- pval_dfs
  de_res$data <- df_full

  return(de_res)
}

#' doDEAnalysis.data.frame
#'
#' @keywords internal
#' @export
doDEAnalysis.data.frame <- function(input, sample_df, design, contrasts, counts,
                                    markers_of_interest, meta_names = NULL,
                                    ctrl_list = NULL) {
  Metacluster <- .id <- NULL

  # Set a seed for reproducibility
  set.seed(42)

  # Get metacluster names if none were given
  # if (is.null(meta_names)) {
  #   meta_names <- levels(fsom$metaclustering)
  # }
  if (is.null(meta_names)) {
    meta_names <- input %>%
      tidytable::pull(Metacluster) %>%
      unique()
  }

  filenames <- input %>%
    tidytable::pull(.id) %>%
    unique() %>%
    basename()

  # Determine samples excluded from analysis, if any
  removed_samples <- which(!(filenames %in% sample_df[, 2]))

  # Create empty data frames for each comparison that will be made. These will
  # store the results of the analysis and be used to create the .csv files.
  pval_dfs <- lapply(1:ncol(contrasts), function(i) {data.frame()})

  # Master data frame that will contain data about all metaclusters.
  df_full <- lapply(1:length(markers_of_interest), function(i) {data.frame()})

  # Initialize object that will be returned
  de_res <- list("tests" = list(), # `pval_dfs`
                 "data" = list()) # `df_full`

  # Prepare data for each marker for testing
  for (m in markers_of_interest) {
    print(paste0("Analyzing ", m, " ..."))

    # Get the channel corresponding to the current control
    channel <- NULL
    fcs_cols <- attr(input, "cols_from_fcs")
    for (i in seq_along(fcs_cols)) {
      col <- fcs_cols[i]
      if (!is.na(attr(fsom_dt[[col]], "marker")) & attr(input[[col]], "marker") == m) { # change so there are no NAs in preprocessing
        channel <- col
        break
      }
    }
    if (is.null(channel)) {
      stop(paste0("No corresponding channel was found for marker ", m,
                  " check parameter `markers_of_interest` for typos."))
    }

    # if (!is.null(controls_df) && m %in% rownames(controls_df)) {
    if (!is.null(ctrl_list) && m %in% names(ctrl_list)) {
      print("Preparing controls...")

      # row <- which(rownames(controls_df) == m)

      ctrl_dt <- ctrl_list[[m]]
      # Pre-existing clustering given
      # if (!is.null(ctrl_fsom_names)) {
      #   # print("Using the FlowSOM objects defined in `ctrl_fsom_names` object...")
      #   # ctrl_fsom <- readRDS(paste0(dir_rds_edited(), ctrl_fsom_names[row]))
      #   stop() # temp
      # } else {
      #   ctrl_dt <- doPreprocessing(controls_df[row, 2],
      #                              ld_channel = ,
      #                              compensation = ,
      #                              transformation = attr(input, "transformation"),
      #                              transformation_type = ,
      #                              debris_gate = ,
      #                              live_gate = ,
      #                              nmad = ,
      #                              pctg_live = 0,
      #                              pctg_qc = 0)
      #   ctrl_fsom <- clusterControlsTable(ctrl_dt,
      #                                     fsom = fsom,
      #                                     subsetted_meta = subsetted_meta)
      # }
    }
    for (k in 1:length(meta_names)) {
      print(paste0("Preparing ", meta_names[k], " ..."))

      # Get cells belonging to each sample in current metacluster
      input_sub <- input %>%
        tidytable::filter(Metacluster == meta_names[k])

      # Get cluster-sample medians
      sample_medians <- getSampleMFIsTable(input_sub) ###
      #row_matches <- match(rownames(sample_medians), sample_info$File.Name)###
      #rownames(sample_medians) <- sample_info[row_matches, "Sample.Name"]###

      # Get transformation
      transform <- attr(input, "transformation")
      if (!is.null(transform)) {
        transform_type <- methods::slot(transform, "transformationId")

        # Reverse transformation if necessary
        if (transform_type == "logicle") { # make switch()?
          inv <- flowCore::inverseLogicleTransform(transform)
        } else if (transform_type == "arcsinh") {
          stop()
        } else if (transform_type == "other") {
          stop()
        } else if (transform_type == "none") {
          stop()
        }

        # Transform sample medians to linear scale
        trans_sample_medians <- flowCore::flowFrame(as.matrix(input_sub, rownames = TRUE)) # loses attributes, consider copying transformed values onto old table?
        trans_sample_medians <- flowCore::transform(trans_sample_medians, inv)
        trans_sample_medians <- flowCore::exprs(trans_sample_medians)

        # Subtract FMO or ISO if necessary
        if (m %in% names(ctrl_list)) {
          # Get only cells belonging to current metacluster
          ctrl_dt_sub <- ctrl_dt %>%
            tidytable::filter(Metacluster == meta_names[k]) # ensure that excluded samples aren't analyzed

          # Transform control values to linear scale
          trans_ctrl_dt_sub <- flowCore::flowFrame(as.matrix(ctrl_dt_sub, rownames = TRUE))
          trans_ctrl_dt_sub <- flowCore::transform(trans_ctrl_dt_sub, inv)
          trans_ctrl_dt_sub <- flowCore::exprs(trans_ctrl_dt_sub) # attributes lost

          # Get cluster-sample medians
          fmo_sample_medians <- getSampleMFIs(trans_ctrl_dt_sub)

          # Calculate delta MFIs
          for (i in 1:nrow(sample_df)) {
            diff <- fmo_sample_medians[rownames(sample_df)[i], channel]
            trans_sample_medians[i, channel] <- trans_sample_medians[i, channel] - diff
          }
        }

        # Get vector of transformed medians for current marker
        trans_matrix <- t(rbind(trans_sample_medians[, channel]))
        trans_matrix <- as.matrix(t(trans_matrix))
        rownames(trans_matrix) <- "Transformed Exp."
      }

      # Checks for samples with too few cells
      # Note: samples not of interest were removed from local variables earlier,
      # before finding `sample_medians`
      missing_samples <- which(!(rownames(sample_df) %in% sample_medians$File))

      # Get vector of medians for current marker
      expr_matrix <- t(rbind(sample_medians[, channel]))
      expr_matrix <- as.matrix(t(expr_matrix))
      rownames(expr_matrix) <- m

      # Remove samples with too few cells from count matrix, if necessary
      # Note: `counts` is found via `sample_df`, so samples not of interest are already removed
      # if (length(missing_samples) > 0) {
      #   count_row <- counts[k, -missing_samples] # does the same need to be done for design?
      # } else {                                      # what if a sample has enough cells in one meta/sample
      #   count_row <- counts[k, ]                 # pair for one marker but not the other?
      # }

      # Make data frame to store info about this metacluster.
      df <- data.frame(t(expr_matrix),
                       t(expr_matrix), # was trans_matrix
                       # counts = count_row,
                       sample = sample_df[which(rownames(sample_medians) %in% sample_df$File.Name), 1],
                       cell_type = meta_names[k],
                       # Group = sample_df$group,
                       check.names = FALSE)

      # print(head(df))

      # Append to greater data frame.
      df_ind <- which(markers_of_interest == m)
      df_full[[df_ind]] <- rbind(df_full[[df_ind]], df)

    }
    colnames(df_full[[df_ind]])[1] <- markers_of_interest[df_ind]
  }

  # Create expression matrix for testing.
  expr_matrix <- data.frame()
  for (i in 1:length(df_full)) {
    print(head(df_full[[i]]))
    new_df <- as.matrix(t(getSampleMetaMatrix(df_full[[i]], 1)))
    expr_matrix <- as.matrix(rbind(expr_matrix, new_df))
  }

  # Create linear models.
  # NOTE: weights questionable
  lm_model <- limma::lmFit(object = expr_matrix,
                           design = design)

  # Perform statistical tests.
  contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
  limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)

  # Create tables containing results of our statistical tests, and add them
  # to the data frame corresponding to the relevant comparison.
  for (i in 1:ncol(contrasts)) {
    table <- limma::topTable(limma_ebayes,
                             coef = colnames(contrasts)[i],
                             sort.by = "p",
                             p.value = 0.20)

    # Add column for marker ID
    if (nrow(table) != 0) {
      table$marker <- rep(markers_of_interest[1], nrow(table))

      if (length(markers_of_interest) > 1) { # if more than one marker is being tested
        for (j in 1:nrow(table)) {
          curr <- rownames(table)[j]
          ix <- ceiling(as.integer(curr)/length(meta_names))
          table$marker[j] <- markers_of_interest[ix]
        }
        rownames(table) <- NULL
      }

      # Append table of results to returned list
      temp_ind <- nrow(pval_dfs[[i]]) + 1
      pval_dfs[[i]] <- rbind(pval_dfs[[i]], table)
    }
  }

  de_res$tests <- pval_dfs
  de_res$data <- df_full

  return(de_res)
}


#' calculateSE
#'
#' Helper function for calculating standard error.
#'
#' @param x A numeric vector.
#'
#' @keywords internal
#'
#' @return Standard error.
#'
#' @export
calculateSE <- function(x) {
  sd <- sd(x)
  se <- sd/sqrt(length(x))
  return(se)
}

#' plotGroupMFIBars
#'
#' Function to draw bar plot of MFIs or dMFIs by group.
#'
#' @param input A matrix, data frame, or path to a .csv file, where each row is
#' a sample and each column is a metacluster.
#' @param sample_df A factor data frame as generated by [prepareSampleInfo()].
#' @param grps_of_interest A list defining which groups to plot.
#' @param upper_lim The y-axis upper limit.
#'
#' @details
#' The names and order of the groups defined in \code{grps_of_interest} should
#' appear in the same way you would like them to appear in the bar plot. Make
#' sure that the group names you supply are the same as they appear in the
#' \code{sample_df} object.
#'
#' !!! should be able to better use tidyr and dplyr here
#'
#' @return A bar plot drawn with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
#'
#' @examples
#' # Define groups of interest
#' grps_of_interest <- list("Male Ctrl" = c("male_Ctrl_X"),
#'                          "Female Ctrl" = c("female_Ctrl_X"),
#'                          "Male MIBC" = c("male_MIBC_No.NAC", "male_MIBC_NAC"),
#'                          "Female MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC"))
plotGroupMFIBars <- function(input, sample_df, grps_of_interest, upper_lim = NULL) {
  Group <- value <- name <- NULL

  # Get input as data frame
  if (is.character(input)) {
    orig_df <- utils::read.csv(input, check.names = FALSE)
  } else if (is.matrix(input)) {
    orig_df <- as.data.frame(input, check.names = FALSE)
  }

  # Add group column to input
  orig_df <- cbind(orig_df, factor_group = sample_df$group) # !!! this column is later removed,
                                                          # binding it shouldn't be necessary
  # Remove sample ID column
  df <- orig_df[,-1]

  # Convert entries to numeric
  for (col in colnames(df)) {
    if (col != "factor_group") {
      df[, col] <- as.numeric(df[, col])
    }
  }

  # Ensure "group" column is treated as a factor.
  df <- df %>%
    dplyr::mutate(factor_group <- as.factor(factor_group))

  # Create new column defining which group of interest each sample belongs to
  df$Group <- factor(seq(1:nrow(df)), levels = names(grps_of_interest))
  for (i in 1:nrow(df)) {
    group <- as.character(df$factor_group[i])
    for (j in 1:length(grps_of_interest)) {
      if (group %in% grps_of_interest[[j]]) {
        df$Group[i] <- names(grps_of_interest)[j]
      }
    }
  }

  # Remove group column not of interest
  df <- df[, -(which(colnames(df) == "factor_group"))]

  # Data frame to plot a data point for each sample
  point_df <- tidyr::pivot_longer(df, cols = -Group)

  # Calculate medians for each group
  median_df <- df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise_if(is.numeric, stats::median, na.rm = TRUE) %>%
    dplyr::ungroup()

  # Calculate standard errors for each group
  se_df <- df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise_if(is.numeric, calculateSE) %>%
    dplyr::ungroup()

  # Pivot data to long format
  med_dat_long <- tidyr::pivot_longer(median_df, cols = -Group)
  se_dat_long <- tidyr::pivot_longer(se_df, cols = -Group)

  # If `upper_lim` is null, approximate appropriate value
  if (is.null(upper_lim)) {
    upper_lim <- max(med_dat_long$value)
    upper_lim <- upper_lim + (upper_lim/2)
  }

  # Draw plot
  ggplot2::ggplot(med_dat_long, ggplot2::aes(x = name, y = value, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
    ggplot2::geom_point(data = point_df,
                        mapping = ggplot2::aes(name, value, shape = Group, fill = Group),
                        position = ggplot2::position_jitterdodge(jitter.width = 0.05, dodge.width = 1, seed = 42),
                        size = 1,
                        inherit.aes = FALSE) +
    ggplot2::geom_errorbar(mapping = ggplot2::aes(x = name,
                                                  ymin = value-se_dat_long$value,
                                                  ymax = value+se_dat_long$value),
                           position = ggplot2::position_dodge(width = 0.9, preserve = "single"),
                           width = 0.5,
                           inherit.aes = TRUE) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "G", begin = 0.2) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 8),
                   aspect.ratio = 0.6) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, upper_lim))
}

#' getDensity
#'
#' Helper for plotting by density
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param ... Additional parameters to pass to [MASS::kde2d()].
#'
#' @keywords internal
#'
#' @return A matrix of the estimated density.
#'
#' @export
getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' plotGroupUMAPs
#'
#' Plot colored UMAPs for each group of interest.
#'
#' @param fsom A FlowSOM object.
#' @param sample_df A factor data frame as generated by [prepareSampleInfo()].
#' @param grps_of_interest A list specifying the groups of interest, in terms
#' of file names.
#' @param umap A UMAP plot as generated by [plotUMAP()]. If \code{NULL}
#' (default), a new UMAP will be generated for the plots.
#' @param color_by A string specifying how you would like the UMAP plot to be colored.
#' Your options are "density" or a marker or channel of interest. Default is "density".
#' @param num_cells The number of cells you would like to be sampled for each
#' group. Default is 5000.
#' @param seed Optional, a seed for reproducibility.
#'
#' @details
#' If any group has fewer than \code{num_cells}, then the number of cells to use
#' for plotting will be taken from the group with the lowest cell count.
#'
#' !!! add check for total cells being greater than number used to make fsom
#' when umap = NULL
#'
#' @return Plots faceted by group drawn with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
plotGroupUMAPs <- function(fsom, sample_df, grps_of_interest, umap = NULL,
                           color_by = "density", num_cells = 5000, seed = NULL) {
  # Set seed if desired
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # If no UMAP given
  if (is.null(umap)) {
    fsom_inds <- c()
    group_vec <- rep(1, length(grps_of_interest) * num_cells)
    for (group in names(grps_of_interest)) {
      # Get sample indices belonging to current group
      # Note that `rownames(sample_df)` corresponds to the sample's number in `fsom$data[, "File"]`
      grp_samples <- rownames(sample_df)[which(sample_df$group %in% grps_of_interest[[group]])]
      # Get cell indices of `fsom` belonging to current group
      inds <- which(fsom$data[, "File"] %in% grp_samples)
      # Sample `num_cells` for current group of interest
      samp <- sample(inds, num_cells)
      fsom_inds <- c(fsom_inds, samp)

      # Assign group name for each cell
      upper_bound <- which(names(grps_of_interest) == group) * num_cells
      grp_inds <- seq((upper_bound - num_cells + 1), upper_bound)
      group_vec[grp_inds] <- group
    }

    # Subset data and generate parent UMAP
    fsom <- FlowSOM::FlowSOMSubset(fsom, fsom_inds)
    umap <- plotUMAP(fsom, length(fsom_inds))

    umap_df <- umap$data

  } else { # if a parent UMAP was given
    umap_df <- umap$data

    # Get group name for each cell
    group_vec <- fsom$data[umap_df$Indices, "File"]
    for (group in names(grps_of_interest)) {
      grp_samples <- rownames(sample_df)[which(sample_df$group %in% grps_of_interest[[group]])] # samples belonging to current group
      inds <- which(group_vec %in% grp_samples) # cell indices of `dat` belonging to current group
      group_vec[inds] <- group

    }

    # Get group counts
    tab <- table(factor(group_vec, levels = names(grps_of_interest)))
    grp_inds <- c()

    # If there is a group with no cells, stop
    if (any(tab == 0)) {
      stop(paste("There are no cells in group", names(tab)[which(tab == 0)],
                 "in the given UMAP. Generate a new one by either setting `umap = NULL`,",
                 "or using the function `plotUMAP()`."))

      # If all groups have more cells than `num_cells`
    } else if (all(tab > num_cells)) {
      subtra <- num_cells
      grp_inds <- 1:length(grps_of_interest)

      # If groups have disproportionate cell counts
    } else if (any(tab > min(tab))) {
      subtra <- min(tab)
      grp_inds <- which(tab > min(tab))

      print(paste("There are some groups with less than", num_cells,
                  "cells. Plotting UMAPs with", subtra, "cells each instead."))
    }
    # Resample each group's cells either according to the `num_cells` argument,
    # or the group with the lowest cell count
    for (grp_ind in grp_inds) {
      diff <- tab[grp_ind] - subtra
      old_inds <- which(group_vec == names(tab[grp_ind]))
      samp <- sample(old_inds, diff) # indices of points to remove

      # Remove sampled cells
      umap_df <- umap_df[-samp, ]
      group_vec <- group_vec[-samp]
    }
  }

  # Initialize parameters
  X1 <- X2 <- values <- NULL

  # Get point values and legend name
  if (color_by == "density") {
    # Calculate density values
    Density <- rep(1, length(group_vec))
    for (group in names(grps_of_interest)) {
      inds <- which(group_vec == group)
      Density[inds] <- getDensity(umap_df[inds, 1], umap_df[inds, 2], n = 100)
    }

    # Append density column to data frame
    umap_df <- data.frame(umap_df, group = group_vec, values = Density/max(Density))
    legend_name <- "Density"

  } else {
    if (color_by %in% colnames(fsom$data)) {
      channel <- color_by
    } else if (color_by %in% FlowSOM::GetMarkers(fsom, colnames(fsom$data))) {
      channel <- FlowSOM::GetChannels(fsom, color_by)
    } else {
      stop("The value given to parameter `color_by` is invalid.")
    }

    # Get expression data and append column to data frame
    marker_vec <- fsom$data[umap_df$Indices, channel]
    umap_df <- data.frame(umap_df, group = group_vec, values = marker_vec)
    legend_name <- channel

  }

  # Draw plot
  p <- ggplot2::ggplot(umap_df) +
    scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = values),
                                  pointsize = 2) +
    ggplot2::facet_wrap(~group) +
    viridis::scale_color_viridis(option = "H", name = legend_name) +
    ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio = 1, strip.text = ggplot2::element_text(size = 12))
  print(p)

  return(p)
}
