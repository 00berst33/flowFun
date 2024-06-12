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
#' from the analysis. Default is none (\code{NULL}).
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
#' have make.names() applied.
#' 
#' *NOTE* dir_prepr() must be set properly for this function to work
#'
#' @return A data frame containing sample information.
#' 
#' @export
prepareSampleInfo <- function(filepath, name_col, filename_col, cols_to_use,
                              samples_to_remove = NULL) {
  # Read in .csv file.
  sample_df <- utils::read.csv(file = filepath,
                               header = TRUE)
  
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
  
  # Get the names of all preprocessed .fcs files. Although file names are included 
  # in the data frame defined above, this variable is included in case any files 
  # were discarded during preprocessing due to issues with the sample, like low 
  # cell viability. 
  # 
  # It is also necessary to include this variable because it reads 
  # the preprocessed files in in the same order as the function we used to create 
  # our aggregate files. 
  prepr_files <- list.files(path = dir_prepr())
  
  # add check that filenames in .csv file exist in the preprocessed dir
  
  # Remove samples that were excluded from the analysis from the data frame 
  # `sample_df`, and reorder its rows such that they are the same as 
  # the file order contained in the `prepr_files` variable.
  if (length(prepr_files) > 0) {
    matched_ind <- match(prepr_files, sample_df[[filename_col]])
    sample_df <- sample_df[matched_ind, ]
    rownames(sample_df) <- seq(1, length(prepr_files))
  } else {
    # stop (edit)
  }
  
  # Remove any samples that are not of interest.
  if (!is.null(samples_to_remove)) {
    sample_df <- sample_df[-samples_to_remove, ]
  }
  
  return(sample_df)
}

#' makeContrastsMatrix
#' 
#' Generate contrasts matrix. (edit example)
#'
#' @param sample_df A data frame, generated either manually or by [prepareSampleInfo()].
#' @param comparisons A named list of named lists, defining the groups to be 
#' compared during analysis. See example for how this variable should be defined.
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
#' contrasts <- makeContrastsMatrix(samples, comparisons)
makeContrastsMatrix <- function(sample_df, comparisons) {
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
  
  # Create design matrix. 
  design_matrix = stats::model.matrix(~ 0 + factors$group)
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
  
  return(Contrasts)
}

#' makeCountMatrix
#' 
#' Generate matrix of sample/metacluster cell counts.
#'
#' @param fsom A FlowSOM object.
#' @param sample_df A data frame containing sample information, generated either
#' manually or by [prepareSampleInfo()].
#' @param meta_names A vector of metacluster names of interest.
#'
#' @return A matrix, where each column represents a sample, and each row
#' represents a metacluster.
#' 
#' @export
makeCountMatrix <- function(fsom, sample_df, meta_names) {
  if (nrow(sample_df) > 0) {
    count_mat <- c()
    
    # For each sample
    for (i in 1:nrow(sample_df)) {
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
    colnames(count_mat) <- sample_df[, 1]
    rownames(count_mat) <- meta_names
  }
}