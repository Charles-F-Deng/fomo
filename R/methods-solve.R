#' Majority-based Sample Relabeling
#' 
#' This heuristic function assigns a Subject_ID to a Genotype_Group_ID if a
#' majority of samples within that Genotype_Group_ID have the same Subject_ID label. 
#' In other words, if most samples in a particular group share the same Subject_ID, 
#' then that Subject_ID is assigned to the entire group.
#' 
#' @param object A MislabelSolver object
#' @param unambiguous_only (Default = FALSE) If true, only correct sample mislabels if they are unambiguous.
#' 
#' @return A MislabelSolver object
#' 
#' @import dplyr
#' @import igraph
#' 
#' @export
#'
solveMajoritySearch <- function(object, unambiguous_only=FALSE) {
    set.seed(1)
    print("Starting majority search")
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        print("0 samples relabeled")
        return(object)
    }
    
    ## 1. Update putative subjects
    votes <- table(object@.solve_state$unsolved_relabel_data$Genotype_Group_ID, 
                   object@.solve_state$unsolved_relabel_data$Subject_ID)
    votes_by_genotype <- data.frame(
        Genotype_Group_ID = rownames(votes),
        Max_Subject_ID = colnames(votes)[apply(votes, 1, which.max)],
        n = rowSums(votes),
        n_Max_Subject_ID = apply(votes, 1, max)
    ) %>%
        dplyr::filter(
            n_Max_Subject_ID >= 2,
            n_Max_Subject_ID > n/2
        ) %>% 
        dplyr::rename(Subject_ID = Max_Subject_ID) %>% 
        dplyr::select(Genotype_Group_ID, Subject_ID)
    votes_by_subject <- data.frame(
        Subject_ID = colnames(votes),
        Max_Genotype_Group_ID = rownames(votes)[apply(votes, 2, which.max)],
        n = colSums(votes),
        n_Max_Genotype_Group_ID = apply(votes, 2, max)
    ) %>% 
        dplyr::filter(
            n_Max_Genotype_Group_ID >= 2,
            n_Max_Genotype_Group_ID > n/2
        ) %>% 
        dplyr::rename(Genotype_Group_ID = Max_Genotype_Group_ID) %>% 
        dplyr::select(Subject_ID, Genotype_Group_ID)
    new_putative_subjects <- dplyr::inner_join(votes_by_subject, votes_by_genotype, by=c("Genotype_Group_ID", "Subject_ID")) %>% 
        dplyr::anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
    object <- .update_putative_subjects(object, new_putative_subjects)
    
    ## 2. Find relabel cycles
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(unsolved_relabel_data, 
                                                            putative_subjects,
                                                            unambiguous_only=unambiguous_only,
                                                            allow_unknowns=FALSE)
    
    ## 3. Relabel samples and update solve state
    object <- .relabel_samples(object, relabels)
    # print(paste(nrow(relabels), "samples relabeled"))
    return(object)
}

#' Comprehensive-based Sample Relabeling
#' 
#' This comprehensive function permutes over all combinations of assigning a 
#' Subject_ID to a Genotype_Group_ID, then picks the assignment that implies the 
#' fewest number of sample mislabels and deletions.
#' 
#' @param object A MislabelSolver object
#' @param max_genotypes (Default = 8) The number of combinations scales in factorial
#'                      with the number of genotypes in the largest connected component.
#'                      The algorithm will skip over all components that exceed this size.
#'                      
#' 
#' @return A MislabelSolver object
#' 
#' @import dplyr
#' @import igraph
#' @import reshape2
#' 
#' @export
#'
solveComprehensiveSearch <- function(object, max_genotypes=8) {
    set.seed(1)
    print("Starting comprehensive search")
    if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
        print("0 samples relabeled")
        return(object)
    }
    
    putative_subjects <- object@.solve_state$putative_subjects
    
    ## 1. Update putative subjects
    component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
    for (component_id in component_ids) {
        cc_unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data %>% 
            dplyr::filter(Component_ID == component_id)
        cc_unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data %>% 
            dplyr::filter(Component_ID == component_id)
        cc_sample_ids <- c(cc_unsolved_relabel_data$Sample_ID, cc_unsolved_ghost_data$Sample_ID)
        cc_swap_cat_ids <- unique(c(cc_unsolved_relabel_data$SwapCat_ID, cc_unsolved_ghost_data$SwapCat_ID))
        cc_genotypes <- unique(cc_unsolved_relabel_data$Genotype_Group_ID)
        cc_subjects <- unique(cc_unsolved_relabel_data$Subject_ID)
        
        ## For now, pass out of components where number of Genotype_Group(s) is greater than number of Subject_ID(s)
        if (length(cc_genotypes) > length(cc_subjects)) {
            # print(component_id)
            next
        }
        
        ## Lock genotypes that already have a putative subject assigned, and find all possible permutations for free genotypes
        locked_genotypes <- intersect(putative_subjects$Genotype_Group_ID, cc_genotypes)
        locked_subjects <- intersect(putative_subjects$Subject_ID, cc_subjects)
        free_genotypes <- setdiff(cc_genotypes, locked_genotypes)
        free_subjects <- setdiff(cc_subjects, locked_subjects)
        
        ## Pass out of components if they have too many genotypes
        if (length(free_genotypes) > max_genotypes | length(free_subjects) > max_genotypes) {
            next
        }
        
        if (length(free_genotypes) > 0 & length(free_subjects) > 0) {
            n <- length(free_subjects)
            r <- length(free_genotypes)
            perm_genotypes <- gtools::permutations(n, r, free_subjects)
            colnames(perm_genotypes) <- sort(free_genotypes)
            n_perm <- nrow(perm_genotypes)
            for (locked_genotype_id in locked_genotypes) {
                locked_subject_id <- putative_subjects[putative_subjects$Genotype_Group_ID == locked_genotype_id, "Subject_ID"][[1]]
                new_perm_col <- matrix(data=locked_subject_id, ncol=1, nrow=n_perm, dimnames=list(NULL, locked_genotype_id))
                perm_genotypes <- cbind(perm_genotypes, new_perm_col)
            }
        } else {
            locked_putative_subjects <- putative_subjects[putative_subjects$Genotype_Group_ID %in% locked_genotypes, ]
            perm_genotypes <- t(locked_putative_subjects$Subject_ID)
            colnames(perm_genotypes) <- locked_putative_subjects$Genotype_Group_ID
        }
        perm_genotypes <- as.matrix(perm_genotypes, dimnames=c("Permutation_ID", "Genotype_Group_ID"))
        n_perms <- nrow(perm_genotypes)
        permutation_ids <- paste0("Permutation", formatC(seq_len(n_perms), width=nchar(n_perms), format="d", flag="0"))
        rownames(perm_genotypes) <- permutation_ids
        
        ## For each Genotype_Group_ID/Subject_ID permutation, determine
        ## 1. The number of existing samples to relabel
        ## 2. The number of ghost samples needed to add
        ## 3. The number of indels required after ghost samples are included
        label_counts <- cc_unsolved_relabel_data %>% 
            dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::group_by(Subject_ID, SwapCat_ID) %>% 
            dplyr::summarize(n_labels = n(), .groups="drop") 
        ghost_label_counts <- cc_unsolved_ghost_data %>% 
            dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::group_by(Subject_ID, SwapCat_ID) %>% 
            dplyr::summarize(n_ghost_labels = n(), .groups="drop") 
        genotype_counts <- cc_unsolved_relabel_data %>% 
            dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::group_by(Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::summarize(n_in_genotype = n(), .groups="drop")
        genotype_subject_concordant_counts <- cc_unsolved_relabel_data %>% 
            dplyr::select(Sample_ID, Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::group_by(Subject_ID, Genotype_Group_ID, SwapCat_ID) %>% 
            dplyr::summarize(n_samples_correct = n(), .groups="drop")
        
        ## Create a "long" version of perm_genotypes
        long_perm_genotypes <- reshape2::melt(perm_genotypes)
        colnames(long_perm_genotypes) <- c("Permutation_ID", "Genotype_Group_ID", "Subject_ID")
        long_perm_genotypes <- as.data.frame(lapply(long_perm_genotypes, as.character))
        
        ## Create an empty matrix of stats for each permutation
        count_cols <- c("n_labels", "n_ghost_labels", "n_in_genotype", "n_samples_correct")
        permutation_stats <- matrix(0, nrow=n_perms, ncol=length(count_cols))
        colnames(permutation_stats) <- count_cols
        rownames(permutation_stats) <- permutation_ids
        for (swap_cat_id in cc_swap_cat_ids) {
            long_perm_genotypes$SwapCat_ID <- swap_cat_id
            
            ## Join in ascending order of size
            merged_long_perm_genotypes <- long_perm_genotypes %>% 
                dplyr::left_join(label_counts, by=c("Subject_ID", "SwapCat_ID")) %>% 
                dplyr::left_join(ghost_label_counts, by=c("Subject_ID", "SwapCat_ID")) %>% 
                dplyr::left_join(genotype_counts, by=c("Genotype_Group_ID", "SwapCat_ID")) %>% 
                dplyr::left_join(genotype_subject_concordant_counts, by=c("Subject_ID", "Genotype_Group_ID", "SwapCat_ID")) %>% 
                dplyr::mutate_at(
                    dplyr::vars(dplyr::all_of(count_cols)),
                    ~dplyr::coalesce(., 0)
                )
            merged_long_perm_genotypes <- merged_long_perm_genotypes %>% 
                dplyr::mutate(
                    ## In each genotype group, the number of samples that can be relabeled to a genotyped sample
                    n_samples_to_relabel = pmin(n_in_genotype, n_labels) - n_samples_correct,
                    ## In each genotype group, the number of mislabeled samples outstanding is
                    ## n_genotypes - n_correct - n_samples_to_relabel. We try to plug this gap with ghost samples
                    n_samples_to_relabel_ghost = pmin(n_in_genotype - n_samples_correct - n_samples_to_relabel, n_ghost_labels),
                    n_label_deletions = pmax(0, n_in_genotype - n_labels - n_ghost_labels),
                    n_genotype_deletions = pmax(0, n_labels - n_in_genotype)
                ) %>% 
                dplyr::arrange(Permutation_ID)
            
            ## Evaluate each permutation
            swap_cat_perm_stats <- merged_long_perm_genotypes %>% 
                dplyr::group_by(Permutation_ID) %>% 
                dplyr::summarize(
                    n_samples_correct = sum(n_samples_correct),
                    n_samples_to_relabel = sum(n_samples_to_relabel),
                    n_samples_to_relabel_ghost = sum(n_samples_to_relabel_ghost),
                    n_genotype_deletions = sum(n_genotype_deletions),
                    n_label_deletions = sum(n_label_deletions),
                    .groups = "drop"               
                ) %>%
                as.data.frame() %>% 
                dplyr::mutate(
                    n_samples_to_relabel = n_samples_to_relabel + pmin(n_genotype_deletions, n_samples_to_relabel_ghost),
                    n_genotype_deletions = pmax(0, n_genotype_deletions - n_samples_to_relabel_ghost),
                    ## The weighting scheme is arbitrary right now
                    perm_score = n_samples_to_relabel + 1.5 * n_samples_to_relabel_ghost + 2 * (n_genotype_deletions + n_label_deletions)
                )
            rownames(swap_cat_perm_stats) <- swap_cat_perm_stats$Permutation_ID
            swap_cat_perm_stats <- swap_cat_perm_stats %>% dplyr::select(-Permutation_ID)
            
            permutation_stats <- permutation_stats + swap_cat_perm_stats
        }
        
        permutation_stats <- permutation_stats %>%
            as.data.frame() %>% 
            mutate(
                Permutation_ID = rownames(.)
            ) %>% 
            arrange(perm_score)
        
        ## To find a single solution, take top row
        best_permutation <- perm_genotypes[permutation_stats$Permutation_ID[1], , drop=FALSE]
        
        ## Record any ties as subject ambiguities
        if (nrow(permutation_stats) > 1) {
            best_score <- permutation_stats$perm_score[1]
            tied_permutation_stats <- permutation_stats %>% 
                filter(perm_score == best_score)
            if (nrow(tied_permutation_stats) > 1) {
                tied_permutations <- perm_genotypes[tied_permutation_stats$Permutation_ID, , drop=FALSE]
                for (curr_genotype_group in colnames(tied_permutations)) {
                    if (!(curr_genotype_group %in% names(object@.solve_state$ambiguous_subjects))) {
                        object@.solve_state$ambiguous_subjects[[curr_genotype_group]] <- tied_permutations[, curr_genotype_group]
                    } else {
                        object@.solve_state$ambiguous_subjects[[curr_genotype_group]] <- 
                            unique(c(object@.solve_state$ambiguous_subjects[[curr_genotype_group]], tied_permutations[curr_genotype_group]))
                    }
                }
            }
        }
        
        new_putative_subjects <- best_permutation %>% 
            t() %>% 
            as.data.frame() %>%
            dplyr::mutate(
                X = rownames(.)
            ) %>% 
            dplyr::relocate(X)
        colnames(new_putative_subjects) <- c("Genotype_Group_ID", "Subject_ID")
        rownames(new_putative_subjects) <- NULL
        
        ## Also update putative_subjects when a Subject_ID in the component 
        ## doesn't have a Genotype_Group_ID, or vice versa
        if (length(cc_genotypes) > length(cc_subjects)) {
            unmatched_genotypes <- setdiff(cc_genotypes, new_putative_subjects$Genotype_Group_ID)
            new_putative_subjects <- rbind(new_putative_subjects, 
                                           data.frame(Genotype_Group_ID = unmatched_genotypes, Subject_ID = NA_character_))
        }
        if (length(cc_genotypes) < length(cc_subjects)) {
            unmatched_subjects <- setdiff(cc_subjects, new_putative_subjects$Subject_ID)
            new_putative_subjects <- rbind(new_putative_subjects, 
                                           data.frame(Genotype_Group_ID = NA_character_, Subject_ID = unmatched_subjects))
        }
        new_putative_subjects <- new_putative_subjects %>%  
            dplyr::anti_join(object@.solve_state$putative_subjects, by=c("Genotype_Group_ID", "Subject_ID"))
        object <- .update_putative_subjects(object, new_putative_subjects)
    }
    
    ## Find relabel cycles
    unsolved_relabel_data <- object@.solve_state$unsolved_relabel_data
    unsolved_ghost_data <- object@.solve_state$unsolved_ghost_data
    putative_subjects <- object@.solve_state$putative_subjects
    relabels <- .find_relabel_cycles_from_putative_subjects(unsolved_relabel_data, putative_subjects, 
                                                            unsolved_ghost_data, allow_unknowns=TRUE)
    
    
    ## Relabel samples and update solve state
    object <- .relabel_samples(object, relabels) 
    # print(paste(nrow(relabels), "samples relabeled"))
    return(object)
}

#' Local Search Sample Relabeling
#' 
#' This search function looks through all possible swaps of 2 samples, and selects
#' the swap that minimizes the sum of within-genotype entropies
#' 
#' @param object A MislabelSolver object
#' @param n_iter (Default = 1) The number of
#' @param include_ghost (Default = FALSE) If TRUE, allow swaps to include ghost samples
#' @param filter_concordant_vertices (Default = FALSE) If TRUE, filter out samples 
#'                                   with at least one concordant edge
#'                      
#' 
#' @return A MislabelSolver object
#' 
#' @import dplyr
#' 
#' @export
#'
solveLocalSearch <- function(object, n_iter=1, include_ghost=FALSE, filter_concordant_vertices=FALSE) {
    set.seed(1)
    print("Starting local search")
    calc_scaled_entropy <- function(x) {
        #n <- sum(x)
        #return(-n*sum(x/n *log(x/n), na.rm=TRUE))
        # n <- sum(x)
        # return(n*sum(log(x/n), na.rm=TRUE))
        return(sum(x*log(x/sum(x)), na.rm=TRUE))
    }
    
    for (i in 1:n_iter) {
        print(paste("Local search iteration (", i, " of ", n_iter, "):: 'include_ghost'=", include_ghost, ", 'filter_concordant_vertices'=", filter_concordant_vertices, sep = ""))
        unsolved_all_data <- rbind(object@.solve_state$unsolved_relabel_data,
                                   object@.solve_state$unsolved_ghost_data)
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            print("0 samples relabeled")
            return(object)
        }
        votes <- table(object@.solve_state$unsolved_relabel_data$Genotype_Group_ID, 
                       object@.solve_state$unsolved_relabel_data$Subject_ID)
        base_entropies <- apply(votes, MARGIN=1, calc_scaled_entropy)
        
        calc_swapped_delta_entropy <- function(swap_from_subject, swap_from_genotype,
                                               swap_to_subject, swap_to_genotype) {
            delta <- 0
            if (!is.na(swap_from_genotype)) {
                genotype_base_entropy <- base_entropies[swap_from_genotype]
                genotype_votes_vec <- votes[swap_from_genotype, ]
                genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] - 1
                genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] + 1
                genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                delta <- delta + genotype_new_entropy - genotype_base_entropy
            }
            if (!is.na(swap_to_genotype)) {
                genotype_base_entropy <- base_entropies[swap_to_genotype]
                genotype_votes_vec <- votes[swap_to_genotype, ]
                genotype_votes_vec[swap_to_subject] <- genotype_votes_vec[swap_to_subject] - 1
                genotype_votes_vec[swap_from_subject] <- genotype_votes_vec[swap_from_subject] + 1
                genotype_new_entropy <- calc_scaled_entropy(genotype_votes_vec)
                delta <- delta + genotype_new_entropy - genotype_base_entropy
            }
            return(delta)
        }
        
        neighbors <- .find_neighbors(object, include_ghost, filter_concordant_vertices) %>% 
            dplyr::left_join(
                unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID")], 
                by=c("Sample_A"="Sample_ID")
            ) %>%
            dplyr::rename(Subject_A = Subject_ID, Genotype_Group_A = Genotype_Group_ID) %>% 
            dplyr::left_join(
                unsolved_all_data[, c("Sample_ID", "Subject_ID", "Genotype_Group_ID", "Component_ID")], 
                by=c("Sample_B"="Sample_ID")
            ) %>% 
            dplyr::rename(Subject_B = Subject_ID, Genotype_Group_B = Genotype_Group_ID) 
        
        print(paste(nrow(neighbors), "candidate swaps being evaluated..."))
        all_component_ids <- sort(unique(object@.solve_state$unsolved_relabel_data$Component_ID))
        relabels <- data.frame(matrix(data=NA, nrow=length(all_component_ids), ncol=2, dimnames=list(c(), c("relabel_from", "relabel_to"))))
        curr_idx <- 1
        for (curr_component_id in all_component_ids) {
            cc_relabel_data <- unsolved_all_data %>% 
                dplyr::filter(Component_ID == curr_component_id)
            cc_neighbors <- neighbors %>% 
                dplyr::filter(Component_ID == curr_component_id)
            
            if (nrow(cc_neighbors) == 0) {next}
            cc_neighbor_objectives <- cc_neighbors %>% 
                dplyr::mutate(
                    delta = mapply(calc_swapped_delta_entropy, 
                                   swap_from_subject=Subject_A, 
                                   swap_from_genotype=Genotype_Group_A,
                                   swap_to_subject=Subject_B,
                                   swap_to_genotype=Genotype_Group_B)
                )
            cc_relabels <- cc_neighbor_objectives %>%
                dplyr::filter(delta > 0, delta == max(delta)) 
            if (nrow(cc_relabels) == 0) {next}
            cc_relabels <- cc_relabels %>% 
                dplyr::sample_n(1) %>%
                dplyr::transmute(
                    relabel_from=Sample_A,
                    relabel_to=Sample_B
                )
            relabels[curr_idx, c("relabel_from", "relabel_to")] <- cc_relabels
            curr_idx <- curr_idx + 1
        }
        
        relabels <- relabels[!is.na(relabels[, 1]) & !is.na(relabels[, 2]), ]
        relabels <- rbind(relabels, data.frame(relabel_from=relabels$relabel_to, relabel_to=relabels$relabel_from))
        object <- .relabel_samples(object, relabels)
        # print(paste(nrow(relabels), "samples relabeled"))
    }
    
    return(object)
}

#' Ensemble Sample Relabeling
#' 
#' This ensemble solver uses a combination of majority-search heuristic, 
#' comprehensive search, and local search to identify and correct mislabels
#' 
#' @param object A MislabelSolver object
#' 
#' @return A MislabelSolver object
#' 
#' @export
#'
solveEnsemble <- function(object) {
    set.seed(1)
    while (TRUE) {
        if (nrow(object@.solve_state$unsolved_relabel_data) == 0) {
            break
        }
        prev_relabel_data <- object@.solve_state$unsolved_relabel_data
        
        object <- solveComprehensiveSearch(object)
        object <- solveMajoritySearch(object)
        object <- solveComprehensiveSearch(object)
        
        comp_relabel_data <- object@.solve_state$unsolved_relabel_data
        object <- solveLocalSearch(object, n_iter=1, include_ghost=TRUE, filter_concordant_vertices=TRUE)
        
        ## If local search found no swaps, try allowing concordant vertices
        if (nrow(comp_relabel_data) == nrow(object@.solve_state$unsolved_relabel_data)) {
            if (identical(comp_relabel_data, object@.solve_state$unsolved_relabel_data)) {
                object <- solveLocalSearch(object, n_iter=1, include_ghost=TRUE, filter_concordant_vertices=FALSE)
            }
        }
        
        if (nrow(prev_relabel_data) == nrow(object@.solve_state$unsolved_relabel_data)) {
            if (identical(prev_relabel_data, object@.solve_state$unsolved_relabel_data)) {
                break
            }
        }
    }
    
    ## After the solve, check if cycles can be broken down
    
    
    return(object)
}