#' The MislabelSolver class
#' 
#' The MislabelSolver class stores the sample metadata required form mislabel detection and correction.
#' 
#' @slot sample_metadata A data.frame containing sample metadata, with one row per sample
#'                       Must include columns for Sample_ID, Subject_ID, and Genotype_Group_ID.
#' @slot swap_cats (Optional) A data.frame with one row per sample specifying the SwapCat_ID,
#'                 where by experimental design only samples with the same SwapCat_ID may be 
#'                 swapped for one another. For example, assay type or batch ID information 
#'                 can be used to categorize Sample_ID(s) into SwapCat_ID(s)
#' @slot anchor_samples (Optional) A character vector of Sample_ID(s) where the label is known to be correct
#' @slot .solve_state A purely internal slot, used to keep track of sample relabels
#' 
#' @return NULL
#' 
#' @import methods
#' 
#' @export
#' 
setClass("MislabelSolver",
         representation(
             sample_metadata = "data.frame",
             swap_cats = "data.frame",
             anchor_samples = "character",
             .solve_state = "list"
         ),
         prototype(
             swap_cats = NULL,
             anchor_samples = character(0)
         )
)


#' Constructor for the MislabelSolver class
#' 
#' @param sample_metadata A data.frame containing sample metadata, with one row per sample. 
#'                       Must include columns for Sample_ID, Subject_ID, and Genotype_Group_ID.
#' @param swap_cats (Optional) A data.frame with one row per sample specifying the SwapCat_ID,
#'                 where by experimental design only samples with the same SwapCat_ID may be 
#'                 swapped for one another. For example, assay type or batch ID information 
#'                 can be used to categorize Sample_ID(s) into SwapCat_ID(s)
#' @param anchor_samples (Optional) A character vector of Sample_ID(s) where the label is known to be correct
#' 
#' @return A MislabelSolver object
#'
#' @import methods
#' 
#' @export
#' 
MislabelSolver <- function(sample_metadata, swap_cats=NULL, anchor_samples=character(0)) {
    ## Convert and validate inputs
    sample_metadata <- as.data.frame(lapply(sample_metadata, as.character))
    .validate_sample_metadata(sample_metadata)
    
    if (is.null(swap_cats)) {
        swap_cats <- sample_metadata[, "Sample_ID", drop=FALSE]
        swap_cats$SwapCat_ID <- "SwapCat1"
    }
    swap_cats <- as.data.frame(lapply(swap_cats, as.character))
    .validate_swap_cats(sample_metadata, swap_cats)
    
    anchor_samples <- unique(as.character(anchor_samples))
    .validate_anchor_samples(sample_metadata, anchor_samples)
    
    return(methods::new("MislabelSolver", sample_metadata, swap_cats, anchor_samples))
}

#' @import dplyr
#' 
setMethod("initialize", "MislabelSolver",
          function(.Object, sample_metadata, swap_cats=NULL, anchor_samples=character(0)) {
              # Hack to get around the NOTE "no visible binding for global variable"
              Genotype_Group_ID <- Subject_ID <- Sample_ID <- Init_Sample_ID <- NULL
              
              ## Provided there are enough shapes, assign a unique shape to each SwapCat_ID
              all_swap_cat_ids <- names(sort(table(swap_cats$SwapCat_ID), decreasing=TRUE))
              swap_cat_shapes <- data.frame(
                  SwapCat_ID = all_swap_cat_ids,
                  SwapCat_Shape = ifelse(length(all_swap_cat_ids) <= length(VISNETWORK_SWAPCAT_SHAPES),
                                         VISNETWORK_SWAPCAT_SHAPES[seq_len(length(all_swap_cat_ids))], "dot"),
                  vertex_size_scalar = 1
              )
              swap_cats <- swap_cats %>% 
                  dplyr::left_join(swap_cat_shapes, by="SwapCat_ID")
              
              ## Initialize object 'solve_state'
              relabel_data <- sample_metadata %>%
                  dplyr::mutate(
                      Init_Sample_ID = Sample_ID,
                      Init_Subject_ID = Subject_ID,
                      Is_Ghost = is.na(Genotype_Group_ID),
                      Is_Anchor = Init_Sample_ID %in% anchor_samples,
                      Solved = FALSE
                  ) %>%
                  dplyr::left_join(swap_cats, by="Sample_ID")
              unsolved_relabel_data <- relabel_data %>% 
                  dplyr::filter(!is.na(Genotype_Group_ID))
              unsolved_ghost_data <- relabel_data %>% 
                  dplyr::filter(is.na(Genotype_Group_ID))
              putative_subjects <- data.frame(Genotype_Group_ID = character(0),
                                              Subject_ID = character(0))
              ambiguous_subjects <- list()
              solve_state <- list(
                  relabel_data = relabel_data,
                  unsolved_relabel_data = unsolved_relabel_data,
                  unsolved_ghost_data = unsolved_ghost_data,
                  putative_subjects = putative_subjects,
                  ambiguous_subjects = ambiguous_subjects
              )
              
              .Object@sample_metadata <- sample_metadata
              .Object@swap_cats <- swap_cats
              .Object@anchor_samples <- anchor_samples
              .Object@.solve_state <- solve_state
              
              .Object <- .update_solve_state(.Object, initialization=TRUE)
              return(.Object)
          }
)
