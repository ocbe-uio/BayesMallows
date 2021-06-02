# Translation to R of C++ and Python code found here
# https://www.geeksforgeeks.org/all-topological-sorts-of-a-directed-acyclic-graph/
all_topological_sorts <- function(graph, path, discovered, n_items){
  flag <- FALSE

  for(i in seq_len(n_items)){
    if(attr(graph, "indegree")[[i]] == 0 && !discovered[[i]]){
      attr(graph, "indegree")[graph[[i]]] <- attr(graph, "indegree")[graph[[i]]] - 1

      path <- c(path, i)
      discovered[[i]] <- TRUE
      all_topological_sorts(graph, path, discovered, n_items)

      attr(graph, "indegree")[graph[[i]]] <- attr(graph, "indegree")[graph[[i]]] + 1
      path <- path[-length(path)]
      discovered[[i]] <- FALSE

      flag <- TRUE
    }
  }
  if(length(path) == n_items) print(path)
}
