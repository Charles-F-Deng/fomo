#' Plot corrections for \code{MislabelSolver} objects
#'
#' Leverages \code{visNetwork} to generate an interactive plot of sample relabels
#'
#' @param object An object of class \code{MislabelSolver}
#' @param path Output directory where solver results are written. Will be created if
#'             directory doesn't already exist. If no directory provided, output will 
#'             be written to the current working directory.
#'
#' @importFrom graphics plot
#' @importFrom openxlsx write.xlsx
#' @import igraph
#' @import dplyr
#' @import visNetwork
#' @import withr
#'
#' @export
#'
writeOutput <- function(object, dir_path=NULL) {
    if (is.null(dir_path)) {
        dir_path <- getwd()
    }
    
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive=FALSE)
        if (dir.exists(dir_path)) {
            print(paste0("Created new directory ", dir_path))
        } else {
            print(paste0("Failed to create new directory ", dir_path, ". Please check if parent directory for the provided path exists."))
        }
    }

    sample_summary <- object@.solve_state$relabel_data %>%
        dplyr::group_by(Genotype_Group_ID, Init_Subject_ID) %>%
        dplyr::mutate(n_agree=n()) %>%
        dplyr::ungroup(Init_Subject_ID) %>%
        dplyr::mutate(n_Genotype_Group=n()) %>%
        dplyr::ungroup(Genotype_Group_ID) %>%
        dplyr::left_join(
            object@.solve_state$putative_subjects %>% 
                filter(!is.na(Genotype_Group_ID)),
            by = "Genotype_Group_ID",
            suffix = c("", "_putative")
        ) %>%
        dplyr::transmute(
            Component_ID = Init_Component_ID,
            Genotype_Group_ID,
            SwapCat_ID,
            Ghost = is.na(Genotype_Group_ID),
            Init_Subject_ID,
            Init_Sample_ID,
            Final_Subject_ID = Subject_ID,
            Final_Sample_ID = Sample_ID,
            Putative_Subject_ID = Subject_ID_putative,
            n_Genotype_Group = ifelse(!Ghost, n_Genotype_Group, NA_integer_),
            Init_Agreement = ifelse(!Ghost, paste0(n_agree-1, " out of ", n_Genotype_Group - 1), NA_character_),
            Relabeled = Init_Sample_ID != Final_Sample_ID,
            Status = case_when(
                Ghost & Relabeled ~ "ghost-relabeled",
                Ghost ~ "ghost",
                is.na(Putative_Subject_ID) ~ "subject_unknown",
                grepl(LABEL_NOT_FOUND, Final_Sample_ID) ~ "deletion_or_duplication",
                Relabeled & (Final_Subject_ID != Putative_Subject_ID | n_Genotype_Group == 1) ~ "relabeled_and_flagged",
                Relabeled ~ "relabeled",
                n_Genotype_Group == 1 ~ "ignored_single-sample",
                n_agree < 2 ~ "ignored",
                TRUE ~ "validated"
            )
        ) %>%
        dplyr::arrange(Component_ID, Genotype_Group_ID, Final_Subject_ID, Final_Sample_ID)

    ## Mislabeled samples in the same genotype group 
    ## with identical swappable categories are
    ## are ambiguities for one another
    ambiguity_summary <- sample_summary %>%
        dplyr::filter(Init_Subject_ID != Putative_Subject_ID) %>%
        dplyr::group_by(Genotype_Group_ID, SwapCat_ID) %>%
        dplyr::mutate(
            n_ambiguities = n() - 1
        ) %>% 
        dplyr::ungroup()
    ambiguity_summary$Ambiguities <- NA_character_
    for (i in seq_len(nrow(ambiguity_summary))) {
        genotype_group_id <- ambiguity_summary$Genotype_Group_ID[i]
        swap_cat_id <- ambiguity_summary$SwapCat_ID[i]
        sample_id <- ambiguity_summary$Final_Sample_ID[i]
        sample_ambiguities <- ambiguity_summary %>%
            filter(
                Genotype_Group_ID == genotype_group_id,
                SwapCat_ID == swap_cat_id,
                Final_Sample_ID != sample_id
            ) %>%
            pull(Final_Sample_ID)
        if (length(sample_ambiguities) > 0) {
            ambiguity_summary[i, "Ambiguities"] <- paste(sample_ambiguities, collapse=", ")
        }
    }
    ambiguity_summary <- ambiguity_summary %>% 
        select(Init_Sample_ID, n_ambiguities, Ambiguities)
    
    sample_summary <- sample_summary %>% 
        left_join(ambiguity_summary, by="Init_Sample_ID")
    
    if (!is.null(object@genotype_matrix)) {
        sample_summary$Neighbor_Connectedness <- NA
        sample_summary$Contamination_Flag <- FALSE
        genotyped_sample_ids <- sample_summary %>% filter(!Ghost) %>% pull(Init_Sample_ID)
        for (sample_id in genotyped_sample_ids) {
            neighbor_samples <- names(which(object@genotype_matrix[sample_id, ] == 1))
            neighbor_matrix <- object@genotype_matrix[neighbor_samples, neighbor_samples]
            existing_edges <- sum(neighbor_matrix[upper.tri(neighbor_matrix)])
            total_edges <- length(neighbor_matrix[upper.tri(neighbor_matrix)])
            sample_summary[sample_summary$Init_Sample_ID == sample_id, "Neighbor_Connectedness"] <-
                paste0(existing_edges, " out of ", total_edges, " pairs of neighbors connected in genotype matrix")
            if (existing_edges < total_edges) {
                sample_summary[sample_summary$Init_Sample_ID == sample_id, "Contamination_Flag"] <- TRUE
            }
        }
    }

    genotype_group_summary <- sample_summary %>%
        group_by(Genotype_Group_ID) %>%
        filter(!is.na(Genotype_Group_ID)) %>% 
        summarize(
            Inferred_Subject_ID = names(sort(table(Final_Subject_ID), decreasing = TRUE)[1]),
            n_Samples_validated = sum(Status == "validated"),
            n_Samples_corrected = sum(Status == "corrected"),
            n_Samples_deletion_or_duplication = sum(Status == "deletion_or_duplication"),
            n_Samples_ignored = sum(str_detect(Status, "^ignored")),
            n_Samples_total = n(),
            Init_Fraction_Match = paste0(n_Samples_validated + n_Samples_ignored, " out of ", n_Samples_total)
        ) %>%
        ungroup() %>%
        select(Genotype_Group_ID, Inferred_Subject_ID, Init_Fraction_Match, everything())
    
    if (!is.null(object@genotype_matrix)) {
        genotype_group_summary$Genotype_Connectedness <- NA
        genotype_group_summary$Genotype_Contamination_Flag <- FALSE
        for (genotype_group_id in genotype_group_summary$Genotype_Group_ID) {
            genotype_group_init_samples <- sample_summary %>%
                filter(Genotype_Group_ID == genotype_group_id) %>%
                pull(Init_Sample_ID)
            genotype_group_matrix <- object@genotype_matrix[genotype_group_init_samples, genotype_group_init_samples]
            existing_edges <- sum(genotype_group_matrix[upper.tri(genotype_group_matrix)])
            total_edges <- length(genotype_group_matrix[upper.tri(genotype_group_matrix)])
            genotype_group_summary[genotype_group_summary$Genotype_Group_ID == genotype_group_id, "Genotype_Connectedness"] <-
                paste0(existing_edges, " out of ", total_edges, " edges")
            if (existing_edges < total_edges) {
                genotype_group_summary[genotype_group_summary$Genotype_Group_ID == genotype_group_id, "Genotype_Contamination_Flag"] <- TRUE
            }
        }
    }

    component_summary <- sample_summary %>%
        group_by(Component_ID) %>%
        summarize(
            n_Genotype_Groups = length(unique(Genotype_Group_ID)),
            n_Subjects = length(unique(Final_Subject_ID)),
            Size_Match = n_Genotype_Groups == n_Subjects,
            Init_Solved = n_Genotype_Groups == 1 & n_Subjects == 1,
            n_Samples_validated = sum(Status == "validated"),
            n_Samples_corrected = sum(Status == "corrected"),
            n_Samples_removed = sum(Status == "removed"),
            n_Samples_ignored = sum(str_detect(Status, "^ignored")),
            n_Samples_total = n()
        )

    dataset_summary <- sample_summary %>%
        summarize(
            n_Components = length(unique(Component_ID)),
            n_Genotype_Groups = length(unique(Genotype_Group_ID)),
            n_Subjects = length(unique(Final_Subject_ID)),
            n_Samples_validated = sum(Status == "validated"),
            n_Samples_corrected = sum(Status == "corrected"),
            n_Samples_removed= sum(Status == "removed"),
            n_Samples_total = n()
        )

    excel_filename <- file.path(dir_path, "corrections_summary.xlsx")
    summary_list <- list(
        "Sample" = sample_summary,
        "Genotype_Group" = genotype_group_summary,
        "Component" = component_summary,
        "Dataset" = dataset_summary)
    openxlsx::write.xlsx(summary_list, file=excel_filename)
    print(paste0("Output successfully written to ", excel_filename))
}

