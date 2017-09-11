#' @title Generates a random tree
#' 
#' @description Given alpha, gamma and the number of leaves n, generates a random phylogenetic tree between all trees with n leaves following the alpha-gamma model.
#' 
#' @param n the number of leaves in the tree.
#' @param alpha parametrer of the alpha-gamma model, between 0 and 1.
#' @param gamma parametrer of the alpha-gamma model, between 0 and alpha.
#' 
#' @return An igraph object that represents the generated phylogenetic tree. 
#' 
#' @references  
#' B. Chen, D. Ford, M. Winkel, A new family of Markov branching trees: the alpha-gamma model. \emph{Electr. J. Probab}. \bold{14} (2009), 400-430. 
#' 
#' @examples  
#' 
#' # A phylogenetic tree with 10 leaves and
#' # parameters alpha=0.8 and gamma=0.1
#' tree = a.g.model(10,0.8,0.1) 
#' # plot(tree,layout=layout.reingold.tilford(tree,root=which(degree(tree,mode="in")==0)))
#' 
#' # A phylogenetic tree with 5 leaves and
#' # parameters alpha=0.5 and gamma=0.3
#' tree = a.g.model(5,0.5,0.3)
#' # plot(tree,layout=layout.reingold.tilford(tree,root=which(degree(tree,mode="in")==0)))
#' 
#' @importFrom igraph graph.edgelist
#' @importFrom igraph degree 
#' 
#' @author Lucia Rotger
#' 
#' @export
a.g.model <-
  function(n,alpha,gamma){
    if(n<3) 
      stop("n<3")
    else{
      if((alpha>1)||(alpha<0)||(gamma>1)||(gamma<0))
        stop("alpha and gamma must been between 0 and 1")
      else{
        if(alpha<gamma)
          stop("alpha < gamma")
        else{  
          edge.matrix = matrix(c(1,2,2,3,2,4),byrow=T,ncol=2) 
          n.nodes = 4
          n.edges = 3
          for(n.leaves in 3:n){#Add new leaf
            #Asign probabilities
            probabilities = rep(0,n.nodes+n.edges)
            degrees = rep(0,n.nodes)
            degree.table = table(edge.matrix[,1]) 
            degrees[as.numeric(names(degree.table))] = degree.table
            
            leaves = which(degrees==0)
            leaf.edge =  which(edge.matrix[,2]%in%leaves)
            
            probabilities[1:n.edges + n.nodes] = gamma
            probabilities[leaf.edge+n.nodes] = 1-alpha
            probabilities[which(degrees>1)] = (degrees[degrees>1]-1)*alpha-gamma
            probabilities = probabilities/(n.leaves-alpha)
            #
            random = sample(c(1:n.nodes,1:n.edges+n.nodes),1,prob=probabilities)
            
            if(random<=n.nodes){#a node is selected
              edge.matrix = rbind(edge.matrix,c(random,n.nodes+1)) 
              n.nodes = n.nodes+1
              n.edges = n.edges+1
            }
            else{#an edge is selected
              random = random - n.nodes
              edge.matrix = rbind(edge.matrix,c(edge.matrix[random,1],n.nodes+1)) 
              edge.matrix = rbind(edge.matrix,c(n.nodes+1,edge.matrix[random,2])) 
              edge.matrix = rbind(edge.matrix,c(n.nodes+1,n.nodes+2)) 
              edge.matrix = edge.matrix[-random,]
              n.nodes = n.nodes+2
              n.edges = n.edges+2
            }
          }
          tree = graph.edgelist(edge.matrix) 
          deg.out = degree(tree,mode="out")
          root.node = which(degree(tree,mode="in")==0)
          if(deg.out[root.node]==1){ #Erase the root-edge
            tree = delete.vertices(tree,root.node)  
          } 
          return(tree)
        }
      }
    }
  }
