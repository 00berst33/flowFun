#' prepareSampleInfo
#'
#' Prepare sample information for DE and DA analysis via a .csv file.
#'
#' @param filepath A filepath to a .csv file.
#' @param name_col The name of the column containing the experiment's sample names.
#' @param filename_col The name of the column containing the experiment's .fcs file names.
#' @param cols_to_use A character vector specifying the column names other than
#' sample name and filename that are relevant for the experiment.
#' @param samples_to_remove A vector defining any samples that should be excluded
#' from the analysis, either because they were discarded during preprocessing
#' or are not of interest. Default is none (\code{NULL}).
#'
#' @details
#' If you have a .csv file specifying sample information, it can be read in and
#' prepared by this function. It MUST have columns for file names and sample
#' names, and any additional columns may specify parameters of interest,
#' like NAC vs. No NAC. Each column should be named after what it represents.
#' This file will be read in as a data frame and will be used to construct the
#' design matrix used for analysis.
#'
#' *REMOVE* differences from orig. script: columns not of interest are removed,
#' rather than keeping all columns and storing that info as a variable. filename
#' remains in the data frame. all columns other than filename and sample names
#' have make.names() applied. instead of ordering and removing rows with
#' prepr_files, they are ordered alphabetically according to filename; this is
#' the same as how list.files() reads in filenames. *add option to override?*
#'
#' @return A data frame containing sample information.
#'
#' @export
prepareSampleInfo <- function(filepath, name_col, filename_col, cols_to_use,
                              samples_to_remove = NULL) {
  # Read in .csv file.
  sample_df <- utils::read.csv(file = filepath,
                               header = TRUE)

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

  # Get only the columns of "sample_df" that are relevant for analysis.
  vars_of_interest <- c(name_col, filename_col, cols_to_use)
  sample_df <- sample_df[, vars_of_interest]

  # Edit column names and row entries to be R friendly, excluding sample and filenames.
  sample_df[, -c(1,2)] <- sapply(sample_df[, -c(1,2)], make.names)

  # Remove any samples that are not of interest.
  if (!is.null(samples_to_remove)) {
    sample_df <- sample_df[-samples_to_remove, ]
  }

  # Reorder rows alphabetically
  sample_df <- sample_df[order(sample_df[[filename_col]]), ]

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
  sample_info <- factors[, -c(1, ncol(factors))]

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



#' prepareControls
#'
#' Prepare isotype or FMO controls for calculating delta MFI.
#'
#' @details
#' Additional details...
#'
#' @export
prepareControls <- function(sample_df, control_col, marker, orig_fsom, fsom_name = NULL,
                            agg_name = NULL, subsetted_meta = NULL) { # agg_name unnecessary if you add metadata
  cols_of_interest <- c(colnames(sample_df)[1], control_cols)
  sample_df <- sample_df[, cols_of_interest]

  ind <- which(names(control_cols) == marker)

  col <- control_col[ind]

  # change the order of the if statements
  if (!is.null(fsom_name) && !(is.null(agg_name))) {
    print("Using the FlowSOM object defined in `fsom_name`...")
    agg_fmo <- flowCore::read.FCS(agg_names[ind])
    fmo_fsom <- readRDS(paste0(dir_rds_edited(), fsom_name[ind]))

  } else if (!is.null(subsetted_meta)) {
    frames <- lapply(list.files(fmo_clustr_dir[ind], full.names = TRUE), flowCore::read.FCS)

    meta_name <- readRDS(paste0(dir_rds_edited(), "full_meta_names.rds")) # !!! edit
    meta_ind <- which(meta_name %in% subsetted_meta)

    frames <- lapply(frames, filterFCS, metaclusters = meta_ind)

    names(frames) <- marker_list[[names(col)]]

    fs <- methods::as(frames, "flowSet")

    agg_fmo <- FlowSOM::AggregateFlowFrames(fs,
                                            cTotal = nrow(orig_fsom$data),
                                            writeOutput = TRUE,
                                            outputFile = paste0(dir_agg(), "aggregate_",
                                                                col, "_", str, ".fcs"))

    fmo_fsom <- FlowSOM::NewData(orig_fsom, agg_fmo)

    saveRDS(fmo_fsom, paste0(dir_rds_edited(), col, "_", str, ".rds"))

  } else {
    # Create aggregate files from cells in FMO/Iso files.
    agg_fmo <- FlowSOM::AggregateFlowFrames(marker_list[[marker]],
                                            cTotal = nrow(orig_fsom$data),
                                            writeOutput = TRUE,
                                            keepOrder = TRUE,
                                            outputFile = paste0(dir_agg(), "aggregate_", col, ".fcs"))

    # Map above aggregate file to FlowSOM object of interest.
    fmo_fsom <- FlowSOM::NewData(orig_fsom, agg_fmo)

    saveRDS(fmo_fsom, paste0(dir_rds_edited(), col, "_fsom.rds"))
    saveRDS(meta_names, paste0(dir_rds_edited(), "full_meta_names.rds"))

    dir_save <- paste0("Data/Clustered ", col, " Files/") # !!!

    if(!dir.exists(dir_save)) {
      dir.create(dir_save)
    }

    FlowSOM::SaveClustersToFCS(fsom = fmo_fsom,
                               originalFiles = marker_list[[names(col)]],
                               outputDir = dir_save)
  }

  # Convert the file name to numbers and store them as vectors. This is necessary
  # because FlowSOM objects represent each file as a number rather than a name.
  fmo_num = match(fmo_info[[col]], unique(fmo_info[[col]]))
}
