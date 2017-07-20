\name{CollessLike-package}
\alias{CollessLike-package}
\alias{CollessLike}
\docType{package}
\title{
\packageTitle{CollessLike}
}
\description{
\packageDescription{CollessLike}
}
\details{

%The DESCRIPTION file:
\packageDESCRIPTION{CollessLike}
\packageIndices{CollessLike}
%~~ An overview of how to use the package, including the most ~~
%~~ important functions ~~
}
\author{
\packageAuthor{CollessLike}

Maintainer: \packageMaintainer{CollessLike}
}
\references{ %Chen, B., Ford, D., Winkel, M., A new family of Markov branching trees: the alpha-gamma model. \emph{Electr. J. Probab}. \bold{14} (2009), 400-430. MR2480547

A. Mir, F. Rossello, L.Rotger, A Colless-like balance index for multifurcating phylogenetic trees.
 
A. Mir, F. Rossello, L.Rotger, A new balance index for phylogenetic trees. Math. Biosc. 241 (2013).
 
M. J. Sackin, "Good" and  "bad" phenograms. Sys. Zool, 21 (1972), 225-226.

\insertRef{alphagamma}{CollessLike}
%~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file ~~
%~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# A random phylogenetic tree of 5 leaves between all trees with 5 leaves 
# following the alpha-gamma model with alpha=0.5 and gamma=0.3
a.g.tree = a.g.model(5,0.5,0.3)
# To compute the percentile of that tree of the Colless-Like,
# Sackin and cophenetic normalized balance indices under the alpha-gamma
# model with alpha=0.5 and gamma=0.3, and a distribution plot.
#distribution(a.g.tree,0.5,0.3,db.path=getwd())
# For a percentile plot set the parameter percentile.plot as TRUE
#distribution(a.g.tree,0.5,0.3,db.path=getwd(),percentile.plot=TRUE)
 
# Computation of a sample of size 100 of the Colless-Like, Sackin and 
# cophenetic balance indices of the distribution of phylogenetic trees 
# with 5 leaves.
  indices.data = indices.simulation(5,0.5,0.3,100)
# Computation of the percentile of the random tree using the previous 
# generated sample  
  distribution(a.g.tree,0.5,0.3,set.indices=indices.data)
}
