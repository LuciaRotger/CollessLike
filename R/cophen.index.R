#' @title Computes the cophenetic balance index of a phylogenetic tree
#' 
#' @description Given a phylogenetic tree, computes the cophenetic balance index of that phylogenetic tree.
#' 
#' @param tree a single phylogenetic tree. It can be entered as a string the Newick format, as a 'phylo' object (\code{ape} package) or as an 'igraph' object (\code{igraph} package). 
#' @param norm a logical variable that indicates whether the index should be normalized or not.
#' 
#' 
#' @details The cophenetic index is computed as the sum of the depths of the least common ancestor (LCA) of every pair of leaves.
#'
#' @return A numeric value.
#' 
#' @references 
#' A. Mir, F. Rossello, L.Rotger, A new balance index for phylogenetic trees. \emph{Mathematical Biosciences} \bold{241} (2013), 125-136.
#' 
#' 
#' @examples 
#' # Computation of the cophenetic balance index of trees 
#' # entered in newick format:
#' cophen.index("(1,2,3,4,5);")
#' cophen.index("(1,(2,(3,(4,5))));")
#' 
#' # Computation of the cophenetic balance index of trees 
#' # entered as a phylo object:
#' require(ape)
#' random.tree = rtree(5,rooted=TRUE)
#' cophen.index(random.tree)
#' 
#' # Computation of the cophenetic balance index of a tree
#' # entered as a igraph object. The tree is randomly 
#' # generated from all trees with 5 leaves following
#' # the alpha-gamma model with alpha=0.5 and gamma=0.3.
#' a.g.tree = a.g.model(5,0.5,0.3)
#' cophen.index(a.g.tree)
#' 
#' #All of them can be normalized (values between 0 and 1)
#' cophen.index("(1,2,3,4,5);",norm=TRUE)
#' cophen.index("(1,(2,(3,(4,5))));",norm=TRUE)
#' cophen.index(random.tree,norm=TRUE)
#' cophen.index(a.g.tree,norm=TRUE)
#' 
#' @importFrom ape read.tree
#' @importFrom igraph graph.edgelist
#' @importFrom igraph degree
#' @importFrom igraph delete.vertices
#' @importFrom igraph get.shortest.paths 
#' 
#' @author Lucia Rotger
#' 
#' @export
cophen.index <-
  function(tree,norm=FALSE){  
    if(class(tree)=="character") 
      tree=read.tree(text = tree)
    if (class(tree)=="phylo") 
      tree=graph.edgelist(tree$edge, directed=TRUE)  
    if(class(tree)!="igraph")
      stop("Not an igraph object. Please introduce a newick string, an ape tree or an igraph tree.")
    root.node = which(degree(tree,mode="in")==0)
    deg.out = degree(tree,mode="out")
    ########
    if(deg.out[root.node]==1){ #exists a root-edge
      tree = delete.vertices(tree,root.node) 
      deg.out = degree(tree,mode="out") 
      root.node = which(degree(tree,mode="in")==0)
    } 
    leaves = which(deg.out==0)
    root.list = get.shortest.paths(tree,root.node)$vpath
    # COPHENETIC # 
    N = length(leaves)
    COPHEN = 0  
    for(i in 1:(N-1))
      for(j in (i+1):N){
        aux  = length(intersect(root.list[[leaves[i]]],root.list[[leaves[j]]]))-1
        COPHEN = COPHEN + aux    
      } 
    if(norm){ 
      max.c = N*(N-1)*(N-2)/6 
      COPHEN = COPHEN/max.c
    }
    return(COPHEN)
  }
