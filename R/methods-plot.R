#' Plot method for \code{MislabelSolver} objects
#'
#' Leverages \code{visNetwork} to generate an interactive plot of sample mislabels
#'
#' @param x An object of class \code{MislabelSolver}.
#' @param y Ignored (included for compatibility with generic \code{plot} function).
#' @param unsolved If \code{TRUE}, plots only samples that are in unsolved components
#' @param collapse_samples If \code{TRUE}, combines samples that are identical in both Subject_ID and Genotype_Group_ID 
#' @param query_by Specifies the field by which to query samples for plotting. Options are: "Init_Component_ID", "Component_ID", "Subject_ID", "Genotype_Group_ID", and "Sample_ID"
#' @param query_val The value to query by. Must not be \code{NULL} if a value is provided in \code{query_by}
#'
#' @importFrom graphics plot
#' @import igraph
#' @import dplyr
#' @import visNetwork
#' @import withr
#'
#' @export
#'
setMethod("plot", signature(x = "MislabelSolver"),
          function(x, 
                   y=NULL, 
                   unsolved=TRUE, 
                   collapse_samples=FALSE,
                   query_by=c("Init_Component_ID", "Component_ID", "Subject_ID", "Genotype_Group_ID", "Sample_ID"),
                   query_val=NULL) {
              
              if (unsolved) {
                  relabel_data <- x@.solve_state$unsolved_relabel_data
                  ghost_data <- x@.solve_state$unsolved_ghost_data
              } else {
                  relabel_data <- x@.solve_state$relabel_data %>% 
                      dplyr::filter(!Is_Ghost)
                  ghost_data <- x@.solve_state$relabel_data %>% 
                      dplyr::filter(Is_Ghost)
              }
              if (!is.null(query_val)) {
                  query_by <- as.character(query_by)
                  query_by <- match.arg(query_by)
                  if (query_by == "Init_Component_ID") {
                      component_id <- rbind(relabel_data, ghost_data) %>% 
                          dplyr::filter(!!sym(query_by) == query_val) %>% 
                          dplyr::pull(Init_Component_ID) %>% 
                          unique()
                  } else {
                      component_id <- rbind(relabel_data, ghost_data) %>% 
                          dplyr::filter(!!sym(query_by) == query_val) %>% 
                          dplyr::pull(Component_ID) %>% 
                          unique()
                  }
                  if (length(component_id) == 0) {
                      warning(paste("No samples found for 'query_by'", paste0("'", query_by, "'"), "and 'query_val'", paste0("'", query_val, "'")))
                      return()
                  }
                  component_id <- component_id[[1]]
                  if (query_by == "Init_Component_ID") {
                      relabel_data <- relabel_data %>% 
                          dplyr::filter(Init_Component_ID == component_id)
                      ghost_data <- ghost_data %>% 
                          dplyr::filter(Init_Component_ID == component_id)
                  } else {
                      relabel_data <- relabel_data %>% 
                          dplyr::filter(Component_ID == component_id)
                      ghost_data <- ghost_data %>% 
                          dplyr::filter(Component_ID == component_id)
                  }
              }
              graph <- .generate_graph(relabel_data, graph_type = "combined", ghost_data, 
                                       populate_plotting_attributes=TRUE, collapse_samples=collapse_samples)
              withr::with_seed(2, {
                  l_mds <- igraph::layout_with_mds(graph)
                  l_drl <- igraph::layout_with_drl(graph, use.seed=TRUE, seed=l_mds)
                  visNetwork::visIgraph(graph, layout = "layout_with_graphopt", start=l_drl) 
              })
          }
)

#' Plot corrections for \code{MislabelSolver} objects
#'
#' Leverages \code{visNetwork} to generate an interactive plot of sample relabels
#'
#' @param object An object of class \code{MislabelSolver}
#'
#' @importFrom graphics plot
#' @import igraph
#' @import dplyr
#' @import visNetwork
#' @import withr
#'
#' @export
#'
plotCorrections <- function(object, init_component_id=NULL) {
    relabels_df <- object@.solve_state$relabel_data %>% 
        filter(Init_Sample_ID != Sample_ID)

    if (!is.null(init_component_id)) {
        relabels_df <- relabels_df %>% 
            filter(Init_Component_ID == init_component_id)
    }
    
    edges_df <- relabels_df %>% select(Init_Sample_ID, Sample_ID)
    corrections_graph <- graph_from_data_frame(edges_df , directed=TRUE)
    
    
    all_samples <- unique(c(relabels_df$Init_Sample_ID, relabels_df$Sample_ID))
    ghost_samples <- relabels_df$
    
    relabels_df <- relabels_df %>% 
        select(Init_Sample_ID, Sample_ID, Is_Ghost, SwapCat_ID, SwapCat_Shape, vertex_size_scalar, is)
    
    
}