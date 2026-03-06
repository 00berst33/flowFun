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

  # Ensure column names for sample and file names are those used by other functions
  sample_df <- sample_df %>%
    dplyr::rename(filename = !!filename_col, sample.name = !!name_col)

  # Remove any samples that are not of interest.
  if (!is.null(samples_to_remove)) {
    sample_df <- sample_df[-(which(sample_df[[name_col]] %in% samples_to_remove)), ]
  }

  # Create a new "group" column that concatenates the factors for comparisons.
  sample_df$group <- do.call(paste, c(sample_df[comp_factors], sep = "_"))

  # Convert the columns in the "factors" data frame to factors.
  sample_df[] <- lapply(sample_df, as.factor)

  return(sample_df)
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
#' @param sample_df If \code{input} is a FlowSOM object, a data frame from
#' [prepareSampleInfo()].
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
#' @seealso [makeCountMatrix.FlowSOM()]
#'
#' @export
makeCountMatrix <- function(input, sample_df = NULL, meta_names = NULL,
                            min_cells = 3, min_samples = NULL) {
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
makeCountMatrix.data.frame <- function(input, sample_df = NULL, meta_names = NULL,
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
                         norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {

  # Create DGEList object
  if (!is.null(norm_method)) { # if user wants to normalize data
    norm_factors <- edgeR::calcNormFactors(counts, method = norm_method)
    dge_list <- edgeR::DGEList(counts = counts,
                               samples = sample_df,
                               group = sample_df$group,
                               norm.factors = norm_factors,
                               remove.zeros = TRUE)
  } else {
    dge_list <- edgeR::DGEList(counts = counts,
                               samples = sample_df,
                               group = sample_df$group,
                               remove.zeros = TRUE)
  }

  # Estimate dispersions. The option `trend_method = "none"` is used, as the
  # dispersion-mean relationship of flow data typically does not resemble that of
  # RNAseq data.
  disp <- edgeR::estimateDisp(dge_list, design = design, trend.method = "none")

  # Fit GLM models
  glm_fit <- edgeR::glmFit(disp, design = design)

  pval_dfs <- lapply(1:ncol(contrasts), function(i) {data.frame()})

  # Perform likelihood ratio tests
  # Note: in edgeR user guide, only one contrast is tested at a time. The same is
  # done here.
  for (i in 1:ncol(contrasts)) {
    glm_lrt <- edgeR::glmLRT(glm_fit, contrast = contrasts[, i])
    top <- edgeR::topTags(glm_lrt, adjust.method = "BH", sort.by = "PValue")

    pval_dfs[[i]] <- top$table
  }

  names(pval_dfs) <- colnames(contrasts)

  return(pval_dfs)
}


#' getExprMatDE
#'
#' Find expression matrix, where rows are a metacluster/marker pair, and columns
#' are samples. Helper for `doDEAnalysis()`
#'
#' @param fsom_dt A data.table with columns for markers/channels, a column `File`
#' denoting the sample a cell is from, and a column `Metacluster` denoting the
#' metacluster a cell belongs to
#' @param marker_cols A character vector of markers/channels of interest; these
#' should be column names of `fsom_dt`
#' @param sample_col The column in `fsom_dt` corresponding to sample ID. Assumed
#' to be `.id` by default
#'
#' The resulting table is passed onto [limma::lmFit()] in this package's typical
#' workflow, but it may also serve as a helpful summary of the data.
#'
#' @return A data.table where columns are samples and rows are metacluster/marker groups
#' @export
getExprMatDE <- function(fsom_dt, marker_cols, sample_col = .id) {
  marker <- median_expr <- feature <- NULL
  id <- rlang::enquo(sample_col)

  collapsed <- fsom_dt %>%
    dplyr::group_by(!!id, Metacluster) %>%
    dplyr::summarize(dplyr::across(tidyr::all_of(!!marker_cols), median), .groups = "drop") %>%
    # reshape long to wide
    tidyr::pivot_longer(
      cols = tidyr::all_of(!!marker_cols),
      names_to = "marker",
      values_to = "median_expr"
    ) %>%
    tidyr::unite("feature", Metacluster, marker, sep = ".") %>%   # e.g., C1.markerA
    tidyr::pivot_wider(
      names_from = !!id,
      values_from = median_expr
    ) %>%
    dplyr::arrange(feature)

  return(collapsed)
}

#' gs_makeMFIMatrix
#'
#' @param gs A `GatingSet`
#' @param cols A `character` vector of the columns in the expression matrix
#' to calculate MFIs for
#' @param subpopulations A `character` vector of subpopulations in `gs`
#' @param inverse `boolean`, whether or not the data should be transformed back
#' to its raw scale. If `TRUE`, a transformation with an associated inverse must be attached
#' to the `GatingSet`
#'
#' @return A data.frame of MFIs
#' @export
gs_makeMFIMatrix <- function(gs, cols, subpopulations, inverse = FALSE) {
  var <- rlang::enquo(cols)

  # Get MFIs from GatingSet
  mfis <- flowWorkspace::gs_pop_get_stats(gs, nodes = subpopulations, type = pop.MFI, inverse = inverse) # may set inverse_transform = TRUE for raw data scale
  mfis <- mfis %>%
    dplyr::select(c(sample, pop, !!var)) %>%
    tidyr::pivot_longer(cols = !c(pop,sample), names_to = "feature") %>%
    tidyr::pivot_wider(names_from = sample, values_from = value) %>%
    dplyr::mutate(feature = paste0(feature, pop), .keep = "unused")

  # Ensure columns are in order corresponding to design matrix
  mfis <- mfis %>%
    dplyr::select(feature, dplyr::all_of(samples))

  return(mfis)
}


#' doDEAnalysis
#'
#' @param input A `GatingSet` or `data.table`
#' @param marker_cols A `character` vector specifying which markers/channels
#' should be tested for differential expression. Each element should be a column
#' in the data's expression matrix.
#' @param design A design matrix from [makeDesignMatrix()]
#' @param contrasts A contrasts matrix from [makeContrastsMatrix()]
#' @param subpopulations Only valid if `input` is a `GatingSet`. A `character` vector
#' specifying which populations in the `GatingSet` to include in tests.
#' @param inverse A `boolean`, whether or not data should be transformed back to
#' a linear scale before calculating MFIs and performing tests. Only valid if `input`
#' is a `GatingSet` with a transformation and its inverse attached.
#'
#' @return An `MArrayLM` object resulting from [limma::eBayes()]
#' @export
doDEAnalysis <- function(input, marker_cols, design, contrasts, subpopulations = NULL, inverse = FALSE) {
  # Find expression matrix: metacluster.marker by sample
  if (is(input, "data.frame")) {
    collapsed <- getExprMatDE(input, marker_cols)
  } else if (is(input, "GatingSet")) {
    # Get MFIs and prepare for limma
    collapsed <- gs_makeMFIMatrix(input,
                                  cols = marker_cols,
                                  subpopulations = subpopulations,
                                  inverse = inverse) # may set to TRUE to get MFIs on scale of raw data

    # Ensure columns are in order corresponding to design matrix
    collapsed <- collapsed %>%
      dplyr::select(feature, dplyr::all_of(samples))
  }

  # Create linear models.
  # NOTE: weights questionable
  lm_model <- limma::lmFit(object = collapsed,
                           design = design)

  # Perform statistical tests.
  contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
  limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)

  return(limma_ebayes)
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


#' @keywords internal
getGroups <- function(comparison, sample_df) {
  # Get the names of all factors
  factor_names <- names(comparison)

  # Generate all combinations of factors
  combinations <- expand.grid(comparison, stringsAsFactors = FALSE)

  # Paste each combination together
  groups <- apply(combinations, 1, function(row) paste(row, collapse = " "))

  # Initialize list to store groups
  sorted_groups <- list()

  # Iterate over all possible found groups
  for (i in seq_len(nrow(combinations))) { # for each combination
    combination <- combinations[i, ]
    samples <- sample_df
    for (j in seq_len(ncol(combinations))) { # for each factor
      var <- colnames(combinations)[j]
      level <- combinations[i, j]

      # Get samples belonging to current factor's level
      samples <- samples %>%
        tidytable::filter(which(samples[[var]] == level))
    }

    # Get files belonging to group
    samples <- samples %>%
      tidytable::pull(2) %>%
      list()

    # Name group and add to list
    names(samples) <- groups[i]
    sorted_groups <- append(sorted_groups, samples)
  }

  return(sorted_groups)
}

#' plotGroupMFIBars
#'
#' Function to draw bar plot of MFIs or dMFIs by group.
#'
#' @param input A matrix or data.frame where rows are samples, columns are clusters,
#' and values are MFIs for a channel
#' @param sample_df A factor data frame as generated by [prepareSampleInfo()].
#' @param comparison A list defining a comparison.
#' @param meta_to_plot A vector of strings defining the metaclusters to include
#' in the plot. Default is all metaclusters.
#' @param upper_lim The y-axis upper limit.
#'
#' @details
#' The names and order of the groups defined in \code{comparison} should
#' appear in the same way you would like them to appear in the bar plot. Make
#' sure that the group names you supply are the same as they appear in the
#' \code{sample_df} object.
#'
#'
#' @return A bar plot drawn with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
plotGroupMFIBars <- function(input, sample_df, comparison, meta_to_plot = NULL, upper_lim = NULL) {
  Group <- value <- name <- NULL

  #df <- getSampleMetaclusterMFIs(input, !!sym(col), sample_df)

  # Get list of groups by file
  grps_of_interest <- getGroups(comparison, sample_df)

  # Get input as data frame
  if (is.matrix(input)) {
    input <- as.data.frame(input, check.names = FALSE)
  }

  if (!is.null(meta_to_plot)) {
    vars <- rlang::enquos(meta_to_plot)
    input <- input %>%
      tidytable::select(!!!vars)
  }

  # Function to find corresponding list names
  find_list_name <- function(filename, list_of_lists) {
    list_name <- names(list_of_lists)[sapply(list_of_lists, function(lst) filename %in% lst)]
    ifelse(length(list_name) > 0, list_name, "X")
  }

  # Add the new column with mutate()
  input <- input %>%
    dplyr::mutate(Group = purrr::map_chr(sample_df[, 2], ~ find_list_name(.x, grps_of_interest))) %>%
    dplyr::filter(Group != "X")

  # Data frame to plot a data point for each sample
  point_df <- tidyr::pivot_longer(input, cols = -Group)

  # Calculate medians for each group
  median_df <- input %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise_if(is.numeric, stats::median, na.rm = TRUE) %>%
    dplyr::ungroup()

  # Calculate standard errors for each group
  se_df <- input %>%
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
#' @importFrom MASS kde2d
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
#' @param input A data.table or FlowSOM object.
#' @param sample_df A factor data frame as generated by [prepareSampleInfo()].
#' @param comparison A list specifying the groups of interest, in terms
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
#'
#' @return Plots faceted by group drawn with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
# Get list of groups by file
plotGroupUMAPs <- function(input, sample_df, comparison, umap = NULL,
                           color_by = "density", num_cells = 5000, seed = NULL)  {
  result <- UseMethod("plotGroupUMAPs")
  return(result)
}

#' @keywords internal
#' @export
plotGroupUMAPs.FlowSOM <- function(input, sample_df, comparison, umap = NULL,
                                   color_by = "density", num_cells = 5000, seed = NULL) {
  Group <- NULL

  # Set seed if desired
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get groups
  grps_of_interest <- getGroups(comparison, sample_df)

  # Function to find corresponding list names
  find_list_name <- function(filename, list_of_lists) {
    list_name <- names(list_of_lists)[sapply(list_of_lists, function(lst) filename %in% lst)]
    ifelse(length(list_name) > 0, list_name, "X")
  }

  # Add the new column with mutate()
  df <- sample_df %>% # to what data frame?
    dplyr::mutate(Group = purrr::map_chr(sample_df[, 2], ~ find_list_name(.x, grps_of_interest))) %>%
    dplyr::filter(Group != "X")

  # If no UMAP given
  if (is.null(umap)) {
    fsom_inds <- c()
    group_vec <- rep(1, length(grps_of_interest) * num_cells)
    for (group in names(grps_of_interest)) { # for each group of interest
      # Get sample indices belonging to current group
      # Note that `rownames(sample_df)` corresponds to the sample's number in `fsom$data[, "File"]`
      grp_samples <- rownames(sample_df)[which(sample_df$File.Name %in% grps_of_interest[[group]])]
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
      grp_samples <- rownames(sample_df)[which(sample_df$File.Name %in% grps_of_interest[[group]])] # samples belonging to current group
      inds <- which(group_vec %in% grp_samples) # cell indices of `dat` belonging to current group
      group_vec[inds] <- group
    }

    # Remove indices that don't belong to any groups of interest
    not_in_groups <- which(!(group_vec %in% names(grps_of_interest)))
    if (length(not_in_groups) > 0) {
      group_vec <- group_vec[-not_in_groups]
      umap_df <- umap_df[-not_in_groups, ]
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

  return(p)
}

#' @keywords internal
#' @export
plotGroupUMAPs.data.table <- function(input, sample_df, comparison, umap = NULL,
                                   color_by = "density", num_cells = 5000, seed = NULL) {
  Group <- NULL

  # Set seed if desired
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get groups
  grps_of_interest <- getGroups(comparison, sample_df)

  # Function to find corresponding list names
  find_list_name <- function(filename, list_of_lists) {
    list_name <- names(list_of_lists)[sapply(list_of_lists, function(lst) filename %in% lst)]
    ifelse(length(list_name) > 0, list_name, "X")
  }

  # Add the new column with mutate()
  df <- sample_df %>% # to what data frame?
    dplyr::mutate(Group = purrr::map_chr(sample_df[, 2], ~ find_list_name(.x, grps_of_interest))) %>%
    dplyr::filter(Group != "X") %>%
    dplyr::select(2, Group)

  # Add column to input
  input <- input %>%
    dplyr::left_join(df, by = dplyr::join_by(.id == filename))

  # If no UMAP given
  if (is.null(umap)) {
    # Sample cells for each group
    res <- input %>%
      dplyr::group_by(Group) %>%
      dplyr::slice_sample(n = num_cells)
    attr(res, "clustered") <- attr(input, "clustered")

    # Generate umap
    umap <- plotUMAP(res, num_cells = nrow(res))

    umap_df <- umap$data

    group_vec <- res %>%
      dplyr::pull(Group)

  } else { # if a parent UMAP was given
    umap_df <- umap$data

    # Get group name for each cell
    group_vec <- input[umap_df$Indices, ] %>%
      dplyr::pull(.id)
    # group_vec <- fsom$data[umap_df$Indices, "File"]
    for (group in names(grps_of_interest)) {
      grp_samples <- sample_df$filename[which(sample_df[, 2] %in% grps_of_interest[[group]])] # samples belonging to current group
      inds <- which(group_vec %in% grp_samples) # cell indices of `dat` belonging to current group
      group_vec[inds] <- group
    }

    # Remove indices that don't belong to any groups of interest
    not_in_groups <- which(!(group_vec %in% names(grps_of_interest)))
    if (length(not_in_groups) > 0) {
      group_vec <- group_vec[-not_in_groups]
      umap_df <- umap_df[-not_in_groups, ]
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
    if (color_by %in% colnames(input)) {
      # Get expression data and append column to data frame
      marker_vec <- input[umap_df$Indices, ] %>%
        dplyr::pull(color_by)
      umap_df <- data.frame(umap_df, group = group_vec, values = marker_vec)
      legend_name <- color_by
    } else {
      stop("The value given to parameter `color_by` is invalid.")
    }
  }

  # Draw plot
  p <- ggplot2::ggplot(umap_df) +
    scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2, color = values),
                                  pointsize = 2) +
    ggplot2::facet_wrap(~group) +
    viridis::scale_color_viridis(option = "H", name = legend_name) +
    ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio = 1, strip.text = ggplot2::element_text(size = 12))

  return(p)
}
