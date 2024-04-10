# This is a workaround for the note "No visible binding for global variable"
# https://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
#' @import utils
#' 
global_vars <- c("Sample_ID", "Subject_ID", "Genotype_Group_ID", "Component_ID",
                 "SwapCat_ID", "Permutation_ID", "n_genotype_deletions", "n_ghost_labels",
                 "n_in_genotype", "n_label_deletions", "n_labels", "n_samples_correct",
                 "n_samples_to_relabel", "n_samples_to_relabel_ghost", "perm_score", ".from",
                 "Init_Component_ID", "Is_Ghost", "Max_Genotype_Group_ID", "Max_Subject_ID",
                 "n_Max_Genotype_Group_ID", "n_Max_Subject_ID", ".", "Col",
                 "Curr_Subject_ID_Genotyped", "Genotype_Group_A", "Genotype_Group_B",
                 "Inferred_Correctly_Labeled", "Invalid_Swap", "Putative_Subject_A",
                 "Putative_Subject_B", "Putative_Subject_ID", "Row", "Sample_A", "Sample_B",
                 "Subject_A", "Subject_B", "SwapCat_A", "SwapCat_B", "X", "delta")
utils::globalVariables(global_vars)