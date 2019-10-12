#' @title Computes Colles-like, Sackin and cophenetic balance indices of a phylogenetic tree.
#' 
#' @description Given a phylogenetic tree, computes Colles-like, Sackin and cophenetic balance indices of that tree.
#' 
#' @param tree a single phylogenetic tree. It can be entered as a string in Newick format, as a 'phylo' object (\code{ape} package) or as an 'igraph' object (\code{igraph} package). 
#' @param norm a logical variable that indicates whether the indices should be normalized or not.
#' @param binary.Colless a logical variable FALSE by default. If is TRUE then the classical Colless index is computed (only for binary trees).
#'
#' @details The Colless-like index is the generalization of the Colless' index for non-binary trees (see Mir et al. , 2018).
#' 
#' The Sackin's index is computed as the sum of the number of ancestors for each leave of the tree (see Mir et al. , 2013).
#' 
#' The cophenetic index is computed as the sum of the depths of the least common ancestor (LCA) of every pair of leaves of the tree(see Sackin et al, 1972).
#'
#' @return A numeric vector with the three computed balance indices of the tree: 
#' \code{Colless-like}, \code{Sackin} and \code{Cophenetic} values.
#' 
#' @references  
#' A. Mir, F. Rossello, L.Rotger, Sound Colless-like balance indices for multifurcating phylogenetic trees.\emph{PloS ONE} 13 (2018), e0203401.
#' 
#' A. Mir, F. Rossello, L.Rotger, A new balance index for phylogenetic trees. \emph{Mathematical Biosciences} \bold{241} (2013), 125-136.
#' 
#' M. J. Sackin, "Good" and  "bad" phenograms. \emph{Sys. Zool}, \bold{21} (1972), 225-226.
#' 
#' @examples 
#' # Computation of the Colless-Like, Sackin and Cophenetic 
#' # balance indices of trees entered in newick format:
#' balance.indices("(1,2,3,4,5);")
#' balance.indices("(1,(2,(3,(4,5))));")
#' 
#' # Computation of the Colless-Like, Sackin and Cophenetic
#' # balance indices of a tree entered as a phylo object:
#' require(ape)
#' random.tree = rtree(5,rooted=TRUE)
#' balance.indices(random.tree)
#' 
#' # Computation of the Colless-Like, Sackin and Cophenetic 
#' # balance indices of a tree entered as a igraph object.
#' # The tree is randomly generated from all trees with 5
#' # leaves following the alpha-gamma model with alpha=0.5
#' # and gamma=0.3.
#' a.g.tree = a.g.model(5,0.5,0.3)
#' balance.indices(a.g.tree)
#' 
#' # All of them can be normalized (values between 0 and 1)
#' balance.indices("(1,2,3,4,5);",norm=TRUE)
#' balance.indices("(1,(2,(3,(4,5))));",norm=TRUE)
#' balance.indices(random.tree,norm=TRUE)
#' balance.indices(a.g.tree,norm=TRUE)
#' 
#' @importFrom ape read.tree 
#' @importFrom igraph graph.edgelist
#' @importFrom igraph degree
#' @importFrom igraph delete.vertices get.shortest.paths  neighborhood V
#' @importFrom stats median
#' 
#' @author Lucia Rotger
#' 
#' @export
balance.indices <-
  function(tree,norm=FALSE,binary.Colless=FALSE){  
    if(class(tree)=="character") 
      tree=read.tree(text = tree)
    if (class(tree)=="phylo") 
      tree=graph.edgelist(tree$edge, directed=TRUE)  
    if(class(tree)!="igraph")
      stop("Not an igraph object. Please introduce a newick string, an ape tree or an igraph tree.")
    root.node = which(degree(tree,mode="in")==0)
    deg.out = degree(tree,mode="out")
    
    # COLLESS.MDM.LN
    D.MDM = function(xx) return(sum(abs(xx-median(xx)))/length(xx))
    f.ln   = function(n)  return(log(n+exp(1))) 
    int.nodes = (1:length(V(tree)))[deg.out>0] #nodes interiores
    decendents = neighborhood(tree,1,int.nodes,mode = "out")  
    fun.nodes.deltas = function(nodes){ 
      aux = neighborhood(tree,length(deg.out)-1,nodes,mode = "out")[[1]]  
      return(sum(f.ln(deg.out[aux]))) 
    }
    fun.children = function(children){
      children = children[-1] 
      result =  unlist(lapply(children,fun.nodes.deltas))
      return(result)
    } 
    deltas = lapply(decendents,fun.children)
    Vdis = lapply(deltas, D.MDM)
    COLLESS = sum(unlist(Vdis)) 
    ########
    if(deg.out[root.node]==1){ #exists root-edge
      tree = delete.vertices(tree,root.node) 
      deg.out = degree(tree,mode="out") 
      root.node = which(degree(tree,mode="in")==0)
    } 
    leaves = which(deg.out==0)
    root.list = get.shortest.paths(tree,root.node)$vpath
    # SACKIN #
    depths = unlist(lapply(root.list,function(xx){length(xx)-1}))
    SACKIN=sum(depths[leaves]) 
    # COPHENETIC # 
    N = length(leaves)
    COPHEN = 0  
    for(i in 1:(N-1))
      for(j in (i+1):N){
        aux  = length(intersect(root.list[[leaves[i]]],root.list[[leaves[j]]]))-1
        COPHEN = COPHEN + aux    
      } 
    result = c("Colles-Like"=COLLESS,"Sackin"=SACKIN,"Cophenetic"=COPHEN)
    if(binary.Colless){
      if(sum(!(deg.out %in% c("0","2")))==0) 
        result[1] = result[1]/( (log(0 + exp(1))+log(2 + exp(1)))/2 )
      else warning("The tree introduced is not binary, Colless-like index for multifurcated trees is computed.")
    } 
    else{
      if(norm){
        max.cl = ( log(0+exp(1)) + log(2+exp(1)) )*(N-1)*(N-2)/4
        max.s = N*(N-1)/2 + N-1
        max.c = N*(N-1)*(N-2)/6
        result[1] = result[1]/max.cl
        result[2] = (result[2]-N)/(max.s-N)
        result[3] = result[3]/max.c
      }
    }
    return(result)
  }
