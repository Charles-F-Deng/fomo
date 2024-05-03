#' @import igraph
#' 
.genotype_matrix_to_genotype_df <- function(genotype_matrix) {
    genotype_graph <- igraph::graph_from_adjacency_matrix(genotype_matrix)
    genotype_group_ids <- igraph::components(genotype_graph)$membership
    n_genotype_groups <- length(unique(genotype_group_ids))
    n_digits <- floor(log10(n_genotype_groups)) + 1
    genotype_group_ids <- vapply(genotype_group_ids, \(x) paste0("Genotype_Group", formatC(x, width = n_digits, format = "d", flag = "0")), "character")
    genotype_df <- data.frame(
        Sample_ID = names(genotype_group_ids),
        Genotype_Group_ID = genotype_group_ids
    )
    return(genotype_df)
}

.diagnose_contaminated_genotype_groups <- function(genotype_matrix) {
    
}

