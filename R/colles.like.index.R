#' @title Computes the Colless-like balance index of a phylogenetic tree
#' 
#' @description Given a phylogenetic tree, computes the Colless-like balance index of that phylogenetic tree.
#' 
#' @param tree a single phylogenetic tree. It can be entered as a string in Newick format, as a 'phylo' object (\code{ape} package) or as an 'igraph' object (\code{igraph} package). 
#' @param norm a logical object indicating if the indices should been normalized or not.
#' @param f.size function to compute the f-size of the tree. See (Mir et al. , 2017) for details . Its default value is "ln" for f(n)=ln(n+e). Other value can be "exp" (f(n)=exp(n). It can also be a user-defined function but in this case, the index cannot be normalized
#' @param diss by default the dissimilarity used to compute the balance index. See (Mir et al. , 2017) for details. Its default value is MDM (mean deviation from the median). Other values can be set as "sd" (sample standard deviation) or "var" (sample variance) . It can also be a user-defined function but in this case the index cannot be normalized.
#' 
#' @details 
#' The Colless-Like balance index is the generalization of the Colless balance index (see Colless,1982) for non-binary trees.
#' 
#' Given a function that computes the f-size of a tree and a dissimarity function that computes the difference of the f-sizes of the subtrees rooted at the children of every internal node of the tree, the Colless-Like index is defined as the sum of these dissimilarities for all internal nodes of the tree. (Mir et al. , 2017)
#' 
#' By default, the f-size function is \code{f(n)=exp(n)} and the dissimilarity is the mean deviation from the median (MDM). 
#' It is possible to change them by specifying it with the parameters \code{f.size} and \code{diss}, with "exp" the f-size would be \code{f(n)=exp(n)}, and with "var" (or "sd") the dissimilarity would be the sample variance (or the sample standard deviation).
#' It is also possible to set a new function for both parameters, see "References".
#' 
#' @return A numeric value.
#' 
#' @references 
#' A. Mir, F. Rossello, L.Rotger, A Colless-like balance index for multifurcating phylogenetic trees. 
#' D. H. Colless, Review of "Phylogenetics: the theory and practice of phylogenetic systematics". Sys. Zool, 31 (1982), 100--104.
#' 
#' @examples 
#' # Computation of the Colless-Like balance index of trees 
#' # entered in newick format:
#' colless.like.index("(1,2,3,4,5);")
#' colless.like.index("(1,(2,(3,(4,5))));")
#' 
#' # Computation of the Colless-Like balance index of trees 
#' # entered as a phylo object:
#' require(ape)
#' random.tree = rtree(5,rooted=TRUE)
#' colless.like.index(random.tree)
#' 
#' # Computation of the Colless-Like balance index of a tree
#' # entered as a igraph object. The tree is randomly 
#' # generated from all trees with 5 leaves following
#' # the alpha-gamma model with alpha=0.5 and gamma=0.3.
#' a.g.tree = a.g.model(5,0.5,0.3)
#' colless.like.index(a.g.tree)
#' 
#' # All of them can be normalized (values between 0 and 1)
#' colless.like.index("(1,2,3,4,5);",norm=TRUE)
#' colless.like.index("(1,(2,(3,(4,5))));",norm=TRUE)
#' colless.like.index(random.tree,norm=TRUE)
#' colless.like.index(a.g.tree,norm=TRUE)
#' 
#' # Computation of the Colless-Like balance index of the
#' # previous generated tree with f-size function f(n)=exp(n):
#' colless.like.index(a.g.tree,f.size="exp")
#' 
#' # Computation of the Colless-Like balance index of the 
#' # previous generated tree that sets the sample variance 
#' # and the sample standard deviation as dissimilarity.
#' colless.like.index(a.g.tree,diss="var")
#' colless.like.index(a.g.tree,diss="sd")
#' 
#' # Computation of the Colless-Like balance index of the 
#' # previous generated tree with f-size function f(n)=exp(n)
#' # that sets the sample variance and the sample standard 
#' # deviation as dissimilarity.
#' colless.like.index(a.g.tree,f.size="exp",diss="var")
#' colless.like.index(a.g.tree,f.size="exp",diss="sd")
#' 
#' @importFrom ape read.tree
#' @importFrom igraph graph.edgelist
#' @importFrom igraph degree
#' @importFrom igraph delete.vertices
#' @importFrom igraph get.shortest.paths 
#' @importFrom stats median
#' 
#' @author Lucia Rotger
#' 
#' @export
colless.like.index <-
  function(tree,f.size="ln",diss="MDM",norm=FALSE){  
    if(class(tree)=="character") 
      tree=read.tree(text = tree)
    if (class(tree)=="phylo") 
      tree=graph.edgelist(tree$edge, directed=TRUE)  
    if(class(tree)!="igraph")
      stop("Not an igraph object. Please introduce a newick string, an ape tree or an igraph tree.")
    root.node = which(degree(tree,mode="in")==0)
    deg.out = degree(tree,mode="out")
     
    case.norm = 0
    if(class(f.size)=="character"){
      if(f.size=="ln"){
        f.size = function(nn)  return(log(nn+exp(1)))
        case.norm = 1
      }
      else if((f.size=="exp")||(f.size=="e")){
        f.size = function(nn)  return( exp(nn) )
        case.norm = 4
      }
      else stop("The f-size introduced is not correct.")
    }
    if(class(diss)=="character"){
      if((diss=="MDM")||(diss=="mdm")){
        diss = function(xx) return(sum(abs(xx-median(xx)))/length(xx))
        case.norm = case.norm*1
      }
      else if(diss=="var"){
        diss = function(xx) return(sum( (xx-mean(xx))^2 )/(length(xx)-1))
        case.norm = case.norm*2
      }
      else if(diss=="sd"){
        diss = function(xx) return(sqrt( sum( (xx-mean(xx))^2 )/(length(xx)-1) ))
        case.norm = case.norm*3
      }
      else stop("The dissimilarity introduced is not correct.")
    } 
    int.nodes = (1:length(V(tree)))[deg.out>0] #nodes interiores
    decendents = neighborhood(tree,1,int.nodes,mode = "out")  
    fun.nodes.deltas = function(nodes){ 
      aux = neighborhood(tree,length(deg.out)-1,nodes,mode = "out")[[1]]  
      return(sum(f.size(deg.out[aux]))) 
    }
    fun.children = function(children){
      children = children[-1] 
      result =  unlist(lapply(children,fun.nodes.deltas))
      return(result)
    } 
    deltas = lapply(decendents,fun.children)
    Vdiss = lapply(deltas, diss)
    COLLESS = sum(unlist(Vdiss)) 
    if(norm){
      if(case.norm==0) warning("Indices can not be normalized")
      else{ 
        N = length(which(deg.out==0))
        # ln MDM
        if(case.norm==1)  max.cl = ( f.size(0) + f.size(2) )*(N-1)*(N-2)/4#( log(0+exp(1)) + log(2+exp(1)) )*(N-1)*(N-2)/4 
        # ln var
        if(case.norm==2)  max.cl = ( f.size(0) + f.size(2) )^2*(N-1)*(N-2)*(2*N-3)/12
        # ln sd
        if(case.norm==3)  max.cl = ( f.size(0) + f.size(2) )*(N-1)*(N-2)/(2*sqrt(2))
        # e^n var
        if(case.norm==8)  max.cl = (f.size(N-1)+N-2)^2/2
        if(N==4){
          # e^n MDM
          if(case.norm==4)  max.cl = (f.size(2)+1)*3/2
          # e^n sd
          if(case.norm==12) max.cl =(f.size(2)+1)*3/sqrt(2)
        }
        else{
          # e^n MDM
          if(case.norm==4)  max.cl = (f.size(N-1)+N-2)/2
          # e^n sd
          if(case.norm==12) max.cl = (f.size(N-1)+N-2)/sqrt(2)
        }
        COLLESS  = COLLESS /max.cl 
      }
    }
    return(COLLESS)
  }
