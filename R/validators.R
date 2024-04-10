#' @importFrom assertthat assert_that
#' 
.validate_sample_metadata <- function(sample_metadata) {
    required_columns <- c("Sample_ID", "Subject_ID", "Genotype_Group_ID")
    missing_columns <- setdiff(required_columns, colnames(sample_metadata))
    assertthat::assert_that(
        length(missing_columns) == 0,
        msg = paste0("'sample_metadata' is missing required column(s) ", paste(missing_columns, collapse=", "))
    )
    
    duplicated_samples <- sample_metadata$Sample_ID[duplicated(sample_metadata$Sample_ID)]
    assertthat::assert_that(
        length(duplicated_samples) == 0,
        msg = paste0("'sample_metadata' has non-unique 'Sample_ID'(s) ", paste(duplicated_samples, collapse=", "))
    )
    
    assertthat::assert_that(
        all(!is.na(sample_metadata$Sample_ID)),
        msg = "'sample_metadata' missing value(s) in 'Sample_ID' column"
    )
    
    assertthat::assert_that(
        all(!is.na(sample_metadata$Subject_ID)),
        msg = "'sample_metadata' missing value(s) in 'Subject_ID' column"
    )
}

#' @importFrom assertthat assert_that
#' 
.validate_swap_cats <- function(sample_metadata, swap_cats) {
    required_columns <- c("Sample_ID", "SwapCat_ID")
    missing_columns <- setdiff(required_columns, colnames(swap_cats))
    assertthat::assert_that(
        length(missing_columns) == 0,
        msg = paste0("'swap_cats' is missing required column(s) ", paste(missing_columns, collapse=", "))
    )
    
    missing_values <- swap_cats[is.na(swap_cats$SwapCat_ID), "Sample_ID"]
    assertthat::assert_that(
        length(missing_values) == 0,
        msg = paste0("'swap_cats' is missing SwapCat_ID for Sample_ID(s) ", paste(missing_values, collapse=", "))
    )
    
    duplicated_samples <- swap_cats$Sample_ID[duplicated(swap_cats$Sample_ID)]
    assertthat::assert_that(
        length(duplicated_samples) == 0,
        msg = paste0("'swap_cats' has non-unique Sample_ID(s) ", paste(duplicated_samples, collapse=", "))
    )
    
    missing_samples <- setdiff(sample_metadata$Sample_ID, swap_cats$Sample_ID)
    assertthat::assert_that(
        length(missing_samples) == 0,
        msg = paste0("'swap_cats' is missing Sample_ID(s) that exist in 'sample_metadata', check Sample_ID(s) ", paste(missing_samples, collapse=", "))
    )
}

#' @importFrom assertthat assert_that
#' 
#' @import dplyr
#' 
.validate_anchor_samples <- function(sample_metadata, anchor_samples) {
    # Hack to get around the NOTE "no visible binding for global variable"
    Genotype_Group_ID <- Subject_ID <- Sample_ID <- n_Subject_ID <- n_Genotype_Group_ID <- NULL
    
    extra_samples <- setdiff(anchor_samples, sample_metadata$Sample_ID)
    assertthat::assert_that(
        length(extra_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) not in 'sample_metadata', check ", paste(extra_samples, collapse=", "))
    )
    
    ## Check that there are no cases where either
    ## 1. Two or more anchor samples with the different Subject_ID(s) have same Genotype_Group_ID
    ## 2. Two or more anchor samples with the same Subject_ID have different Genotype_Group_ID(s)
    anchor_samples_consistency <- data.frame(Sample_ID = anchor_samples) %>% 
        dplyr::left_join(
            sample_metadata[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")],
            by="Sample_ID"
        ) %>% 
        dplyr::filter(!is.na(Genotype_Group_ID)) %>% 
        dplyr::group_by(Genotype_Group_ID) %>% 
        dplyr::mutate(
            n_Subject_ID = length(unique(Subject_ID))
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(Subject_ID) %>% 
        dplyr::mutate(
            n_Genotype_Group_ID = length(unique(Genotype_Group_ID))
        ) %>% 
        dplyr::ungroup()
    subject_inconsistent_samples <- anchor_samples_consistency %>% 
        dplyr::filter(n_Subject_ID != 1) %>% 
        dplyr::arrange(Sample_ID) %>% 
        dplyr::pull(Sample_ID)
    genotype_inconsistent_samples <- anchor_samples_consistency %>% 
        dplyr::filter(n_Genotype_Group_ID != 1) %>% 
        dplyr::arrange(Sample_ID) %>% 
        dplyr::pull(Sample_ID)
    assertthat::assert_that(
        length(subject_inconsistent_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) that have the same Genotype_Group_ID
                   but different Subject_ID(s), check ", paste(subject_inconsistent_samples, collapse=", "))
    )
    assertthat::assert_that(
        length(genotype_inconsistent_samples) == 0,
        msg = paste0("'anchor_samples' contains Sample_ID(s) that have the same Subject_ID
                   but different Genotype_Group_ID(s), check ", paste(genotype_inconsistent_samples, collapse=", "))
    )
}

# --------------------------------------------------------------------------
# Don't need validators for internally created variables 'putative subjects' and 'relabels'


# .validate_putative_subjects <- function(sample_metadata, putative_subjects) {
#     required_columns <- c("Genotype_Group_ID", "Subject_ID")
#     missing_columns <- setdiff(required_columns, colnames(putative_subjects))
#     assert_that(
#         length(missing_columns) == 0,
#         msg = glue("'putative_subjects' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
#     )
#     
#     extra_genotype_groups <- setdiff(na.omit(putative_subjects$Genotype_Group_ID), sample_metadata$Genotype_Group_ID)
#     assert_that(
#         length(extra_genotype_groups) == 0,
#         msg = glue("'putative_subjects' has 'Genotype_Group_ID'(s) not found in 'sample_metadata', check {paste(extra_genotype_groups, collapse=\", \")}")
#     )
#     extra_subjects <- setdiff(na.omit(putative_subjects$Subject_ID), sample_metadata$Subject_ID)
#     assert_that(
#         length(extra_subjects) == 0,
#         msg = glue("'putative_subjects' has 'Subject_ID'(s) not found in 'sample_metadata', check {paste(extra_subjects, collapse=\", \")}")
#     )
#     
#     duplicated_genotype_groups <- putative_subjects$Genotype_Group_ID[duplicated(na.omit(putative_subjects$Genotype_Group_ID))]
#     assert_that(
#         length(duplicated_genotype_groups) == 0,
#         msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, check 'Genotype_Group'(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
#     )
#     duplicated_subjects <- putative_subjects$Subject_ID[duplicated(na.omit(putative_subjects$Subject_ID))]
#     assert_that(
#         length(duplicated_subjects) == 0,
#         msg = glue("'putative_subjects' does not map 'Subject_ID' to 'Genotype_Group_ID' one-to-one, 'check Subject_ID'(s) {paste(duplicated_samples, collapse=\", \")}")
#     )
# }
# 
# .validate_relabels <- function(sample_metadata, relabels) {
#     required_columns <- c("relabel_from", "relabel_to")
#     missing_columns <- setdiff(required_columns, names(relabels))
#     assert_that(
#         length(missing_columns) == 0,
#         msg = glue("'relabels' is missing required column(s) {paste(missing_columns, collapse=\", \")}")
#     )
#     
#     missing_relabel_from <- setdiff(relabels$relabel_from, sample_metadata$Sample_ID)
#     assert_that(
#         length(missing_relabel_from) == 0,
#         msg = glue("'relabels$relabel_from' has Sample_ID(s) not in 'sample_metadata', check {paste(missing_relabel_from, collapse=\", \")}")
#     )
#     missing_relabel_to <- setdiff(relabels$relabel_to, sample_metadata$Sample_ID)
#     assert_that(
#         length(missing_relabel_to) == 0,
#         msg = glue("'relabels$relabel_to' has Sample_ID(s) not in 'sample_metadata', check {paste(missing_relabel_to, collapse=\", \")}")
#     )
#     
#     duplicated_relabel_from <- relabels$relabel_from[duplicated(relabels$relabel_from)]
#     assert_that(
#         length(duplicated_relabel_from) == 0,
#         msg = glue("'relabels$relabel_from' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_relabel_from, collapse=\", \")}")
#     )
#     duplicated_relabel_to <- relabels$relabel_to[duplicated(relabels$relabel_to)]
#     assert_that(
#         length(duplicated_relabel_to) == 0,
#         msg = glue("'relabels$relabel_to' has non-unique Sample_ID(s), check Sample_ID(s) {paste(duplicated_genotype_groups, collapse=\", \")}")
#     )
#     
#     creating_duplicates <- setdiff(relabels$relabel_to, relabels$relabel_from)
#     assert_that(
#         length(creating_duplicates) == 0,
#         msg = glue("'relabels$relabel_to' contains Sample_ID(s) that don't exist in 'relabels$relabel_from', check Sample_ID(s) {paste(creating_duplicates, collapse=\", \")}")
#     )
# }