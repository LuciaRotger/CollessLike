#' @title Generates random trees and compute their balance indices
#' 
#' @description Generates a list of trees according to the introduced parameters for the alpha-gamma model. Then, this 3 balance index are calulated: Colless-like, Sackin and Cophenetic.
#'  
#' @param n the number of leaves of the tree.
#' @param alpha parametrer of the alpha-gamma model, between 0 and 1.
#' @param gamma parametrer of the alpha-gamma model, between 0 and alpha.
#' @param repetitions the number of trees to generate.
#' @param norm a logical object indicating if the indices should been normalized or not. 
#'
#' @details Given a number of leaves, the function generates a tree with that number of leaves and computates the three indeces of balance (Colles-like, Sackin and Cophenetic with function \code{\link{balance.indices}}). This is done as many times as it is set by 'repetitions' parameter, and it generates a 3-column data.frame of indices. 
#' 
#' The trees are generated according to the alpha-gamma model. These parameters can be specified by \code{alpha} and \code{gamma} parameters of the function. The following cases are distinguished:
#' \itemize{
#'  \item{\code{alpha = NA } and \code{ gamma  = NA} :} {All the 66 combinations of \code{alpha} in \{{ 0, 0.1, 0.2, ... ,0.9, 1 \}} and \code{gamma} in \{{ 0, 0.1, ... ,\code{alpha} \}} are done.}
#'  \item{\code{alpha} in [0,1] and \code{gamma  = NA} :} { Since \code{alpha} is fixed, all the combinations with that \code{alpha} and \code{gamma} in \{{ 0, 0.1, ... ,\code{alpha} \}} are done.}
#'  \item{\code{alpha} in [0,1] and \code{gamma} in [0,\code{alpha}] :} { Both parameters are fixed. Then, only that combination is done.}
#' } 
#'
#' @return A 3-column data.frame with the Colless-like, Sackin and Cophenetic balance indices for every generated tree. If more than one data.frame has to be generated, then the returned value is a data.frame list (its names specify which alpha and gamma parameters have generated that data.frame, for instance "a0.5g0.3" indicates alpha=0.5 and gamma=0.3). 
#' 
#' @seealso \code{\link{balance.indices}}
#' 
#' @references  
#' B. Chen, D. Ford, M. Winkel, A new family of Markov branching trees: the alpha-gamma model. \emph{Electr. J. Probab}. \bold{14} (2009), 400-430. 
#' 
#' A. Mir, F. Rossello, L.Rotger, A Colless-like balance index for multifurcating phylogenetic trees.\emph{}
#' 
#' A. Mir, F. Rossello, L.Rotger, A new balance index for phylogenetic trees. \emph{Mathematical Biosciences} \bold{241} (2013), 125-136.
#' 
#' M. J. Sackin, "Good" and  "bad" phenograms. \emph{Sys. Zool}, \bold{21} (1972), 225-226.
#'  
#'  
#' 
#' @examples #('Repetitions' set as 100 for a fast example)
#' \donttest{indices.table = indices.simulation(5,0.5,0.3,repetitions=10)}
#' \donttest{head(indices.table)}
#' 
#' #Normalized indices (between 0 and 1)
#' \donttest{indices.table = indices.simulation(5,0.5,0.3,repetitions=10,norm=TRUE)}
#' \donttest{head(indices.table)}
#' 
#' #Without specifying alpha and gamma
#' \donttest{indices.list = indices.simulation(5,repetitions=100)}
#' #by default alpha=seq(0,1,0.1) and gamma=seq(0,alpha,0.1), thus
#' \donttest{length(indices.list) #=66}
#' #all the elements of the list have a name that identifies its parameters
#' \donttest{indices.list$a0.5g0.3}
#' \donttest{indices.list$a0.7g0.2}
#'  
#' @author Lucia Rotger
#' 
#' @export
indices.simulation <-
function(n,alpha=NA,gamma=NA,repetitions=1000,norm=FALSE){ 
    only.one=FALSE
    if(is.na(alpha)){
        parameters = expand.grid(seq(0,1,0.1),seq(0,1,0.1),n)
        parameters = parameters[which(parameters[,1]>=parameters[,2]),]
    }    
    else{
      if(is.na(gamma)){
        parameters = expand.grid(alpha,seq(0,alpha,0.1),n)
        parameters = parameters[which(parameters[,1]>=parameters[,2]),]
      }
      else{ 
        if((alpha>1)||(alpha<0)||(gamma>1)||(gamma<0))
          stop("alpha and gamma must been between 0 and 1")
        else
          if(alpha<gamma)
            stop("alpha < gamma")
          else{
            parameters = c(alpha,gamma,n)
            only.one = TRUE
          }
      }
    }
    generator = function(idx,n,alpha,gamma){return(a.g.model(n,alpha,gamma))}
    iterate.ford = function(tab){ #tab = [alpha,gamma,n]
        if(tab[1]>=tab[2]){
            alpha = tab[1]
            gamma = tab[2]
            n = tab[3]
            print(paste("n :",n," alpha :",alpha," gamma :",gamma))
            tree.list = lapply(1:repetitions,generator,n,alpha,gamma)  
            result = matrix(unlist(lapply(tree.list,balance.indices,norm=norm)),ncol=3,byrow=T)
            colnames(result) = c("COLLES.MDM.LN","SACKIN","COPHENETIC") 
        }
        return(result)
    } 
    if(only.one){
      result = iterate.ford(parameters) 
    }
    else{
      parameters2=lapply(1:(dim(parameters)[1]), function(i) as.numeric(parameters[i,]))
      result = lapply(parameters2,iterate.ford)
      paste.param = function(tab){return(paste("a",tab[1],"g",tab[2],sep=""))}
      names(result) = apply(parameters,1,paste.param)
    }
    return(result)
}
