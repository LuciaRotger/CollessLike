#' @title Computes the cophenetic vector of a phylogenetic tree
#' 
#' @description Given a phylogenetic tree, computes the cophenetic vector of that phylogenetic tree.
#' 
#' @param tree a single phylogenetic tree. It can be entered as a string the Newick format, as a 'phylo' object (\code{ape} package) or as an 'igraph' object (\code{igraph} package). 
#' @param set.of.labels a list containing the set of labels of the tree (in the order of appearance from its newicks form). By default is \code{NULL} due to if \code{tree} is a \code{phylo} object or a newick string it is no necessary to specify it.
#' 
#' @details The cophenetic vector is formed by all the cophenetic values between each pair of leaves. The cophenetic value of two different leaves i and j is the depth of the lowest common ancestor of the nodes in the tree labeled with i and j. If i=j, then the cophenetic value is the depth of that leaf.
#'
#' @return A numeric array.
#' 
#' @references 
#' A. Mir, F. Rossello, L.Rotger, A new balance index for phylogenetic trees. \emph{Mathematical Biosciences} \bold{241} (2013), 125-136.
#' 
#' G. Cardona, A. Mir, F. Rossello, L. Rotger, D. Sanchez, Cophenetic metrics for phylogenetic trees, after Sokal and Rohlf. \emph{BMC Bioinformatics} \bold{14} (2013), 3.
#' 
#' @examples 
#' # Computation of the cophenetic vector of trees 
#' # entered in newick format:
#' cophen.vector("(1,2,3,4,5);")
#' cophen.vector("(1,(2,(3,(4,5))));")
#' 
#' # Computation of the cophenetic vector of trees 
#' # entered as a phylo object:
#' require(ape)
#' random.tree = rtree(5,rooted=TRUE)
#' cophen.vector(random.tree) 
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
cophen.vector <-
  function(tree,set.of.labels=NULL){  
    if(class(tree)=="character")
      tree=read.tree(text = tree)
    if (class(tree)=="phylo"){
      if(is.null(set.of.labels))
        set.of.labels = tree$tip.label 
      tree=graph.edgelist(tree$edge, directed=TRUE)  
    }
    if(class(tree)!="igraph")
      stop("Not an igraph object. Please introduce a newick string, an ape tree or an igraph tree.")
    if(is.null(set.of.labels))    
      stop("Please insert the set of labels or a phylo object")
    root.node = which(degree(tree,mode="in")==0)
    deg.out = degree(tree,mode="out")
    ########
    if(deg.out[root.node]==1){ #exists a root-edge
      tree = delete.vertices(tree,root.node) 
      deg.out = degree(tree,mode="out") 
      root.node = which(degree(tree,mode="in")==0)
    } 
    leaves = which(deg.out==0)
    if(length(set.of.labels)!=length(leaves))
      stop("Please insert the correct set of labels or a phylo object")
    root.list = get.shortest.paths(tree,root.node)$vpath
    # COPHENETIC # 
    N = length(leaves)
    COPHEN = c()
    ordered.leaves=order(set.of.labels)
    for(i in 1:(N-1)){
      leaf.i = ordered.leaves[i]
      COPHEN = c(COPHEN,length(root.list[[leaf.i]])-1)
      for(j in (i+1):N){
        leaf.j = ordered.leaves[j]
        aux  = length(intersect(root.list[[leaves[leaf.i]]],root.list[[leaves[leaf.j]]]))-1
        COPHEN = c(COPHEN,aux)    
      }
    }  
    COPHEN = c(COPHEN,length(root.list[[ordered.leaves[N]]])-1)
    return(COPHEN) 
  }
