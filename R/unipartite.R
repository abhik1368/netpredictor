
#" Adamic-Adar index
#" @description Similarity measure counting common neighbors weighted by their degrees.
#" @param graph: a igraph object
aaSim <- function(graph){
    score <- igraph::similarity.invlogweighted(graph)
    return (score)
}

## Average commute time
## @description This function calculates vertex similarity based on the average number of
## steps, that random walker on the graph needs to get from one vertex to
## another.
## @param graph: a igraph object


actSim <- function(graph){
    L <- igraph::graph.laplacian(graph)
    n <- igraph::vcount(graph)
    m <- igraph::ecount(graph)
    
    L_psinv <- solve(L - 1/n) + 1/n
    score <- 2 * m * (diag(L_psinv) %*% t(rep(1, n)) +
                          rep(1, n) %*% t(diag(L_psinv)) -
                          2 * L_psinv)
    score <- 1 / score
    return (score)
}

## Common neighbors vertex similarity
## @description Similarity measure counting number of common neighbors.
## @param graph: a igraph object

cnSim <- function(graph){
    score <- igraph::cocitation(graph)
    return (score)
}

## Jaccard Index similarity
## @description Similarity measure based on jaccard index.This is a simple wrapper to an \pkg{igraph} function \code{\link[igraph]{similarity.jaccard}},
#  which counts the proportion of neighbours shared by two vertices.
## @param graph: a igraph object

jcSim <- function(graph){
    score <- igraph::similarity.jaccard(graph)
    return (score)
}

## Dice similarity
## @description Similarity measure based on dice similarity. This function measures a relative size of an intersection of neighbors sets
#  of two vertices.
## @param graph: a igraph object

diceSim <- function(graph){
    score <- igraph::similarity.dice(graph)
    return (score)
}


## Katz Index similarity
## @description Similarity measure based on all paths in a graph. This function counts all the paths between given pair of nodes, with shorter
#  paths counting more heavily. Weigths are exponential.
## @param graph: a igraph object


katzSim <- function(graph,beta = 0.001){
    A  <- igraph::get.adjacency(graph)
    I <- diag(nrow(A))
    tmp <- I - beta * A
    score <- solve(tmp)
    diag(score) <- 0
    return (score)
}
## Pseudoinverse of the Laplacian
## @description Similarity measure based solely on the pseudoinverse of the Laplacian matrix. This function counts all the paths between given pair of nodes, with shorter
#  paths counting more heavily. Weigths are exponential.
## @param graph: a igraph object

invLSim <- function(graph){
    L <- igraph::graph.laplacian(graph)
    n <- vcount(graph)
    
    score <- solve(L - 1/n) + 1/n
    return (score)
}

## Geodesic distance vertex similarity
## @description This function calculates similarity score for vertices based on the shortest paths between them.
## @param graph: a igraph object

distSim <- function(graph){
    score <- igraph::shortest.paths(graph)
    score[is.infinite(score)] <- igraph::vcount(graph) + 1
    score <- 1 / score
    diag(score) <- 0
    return (score)
}

## Cosine vertex similarity/ Salton index
## @description This function measures the cosine of the angle between columns of
## the adjacency matrix, corresponding to given nodes.
## @param graph: a igraph object

cosineSim <- function(graph){
    x <- igraph::get.adjacency(graph)
    return (x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2)))))
}

## Preferential attachment vertex similarity
## @description This function calculates similarity score for vertices based on their degrees.
## @param graph: a igraph object
paSim <- function(graph){
    deg <- igraph::degree(graph)
    score <- outer(deg,deg)
    diag(score) <- 0
    return (score)
}

## Local Paths Index
## @description This function counts the number of two-paths and three-paths between nodes,
#  with three-paths weighted by a parameter \eqn{\epsilon}.
## @param graph: igraph object
## @param eps: epsilon parameter
## @return a matrix object of similarity scores.
## @param graph: a igraph object

lpSim <- function(graph,eps = 0.01){
    A <- igraph::get.adjacency(graph)
    score <- A %*% A
    score <- score + score %*% A * eps
    return (score)
}

## Leicht-Holme-Newman Index
## @description This index gives high similarity to node pairs that have many common
## neighbors compared not to the possible maximum, but to the expected number of such neighbors
## @param graph: a igraph object
lhnSim <- function(graph){
    deg <- igraph::degree(graph)   
    score <- igraph::cocitation(graph)
    score <- score / outer(deg, deg)
    return (score)
}
    
## Hub Depressed Index
## @description This measures assigns lower scores to links adjacent to hubs (high degree nodes). It counts common neighbors of two vertices and weigths the result
#  by the higher of those vertices# degrees.
## @param graph: a igraph object

hdiSim <- function(graph){
    deg <- igraph::degree(graph)
    
    score <- igraph::cocitation(graph)
    score <- score / outer(deg, deg, pmin)
    score[is.nan(score)]<-0
    return (score)
}

## Hub promoted Index
## @description This measures assigns higher scores to links adjacent to hubs (high degree nodes). It counts common neighbors of two vertices and weigths the result
#  by the lower of those vertices# degrees.
## @param graph: a igraph object

hpiSim <- function(graph){
    deg <- igraph::degree(graph)
    
    score <- igraph::cocitation(graph)
    score <- score / outer(deg,deg,pmax)
    score[is.nan(score)]<-0
    return (score)
}

## Similarity measure based on resource allocation process.
## @param graph: igraph object 
## @description This function counts the number of common neighbours weighted by the inverse
## of their degrees.
## @param graph: a igraph object

raSim <- function(graph){
    
    n <- igraph::vcount(graph)
    score <- matrix(integer(n^2), nrow = n)   
    neighbors <- igraph::neighborhood(graph, 1)
    neighbors <- lapply(neighbors, function(x) x[-1])
    
    degrees <- igraph::degree(graph)
    for (k in seq(n)){
        tmp <- neighbors[[k]]
        l <- degrees[[k]]
        if (l > 1){
            for (i in 1:(l-1)){
                n1 <- tmp[i]
                for (j in (i+1):l){
                    n2 <- tmp[j]
                    score[n1, n2] <- score[n1, n2] + 1/l
                    score[n2, n1] <- score[n2, n1] + 1/l
                }
            }
        }
    }
    return(score)
}

## function for matching methods
match_method <- function(method){
    method <- match.arg(method, c("act", "aa", "cn", "jc", "dice","invL", "hdi", "hpi", "katz","dist",
                                  "lhn", "lp", "cosine", "pa","ra"))
    paste0(method,"Sim")
}

#' Get similarity using various link prediction methods
#' @title unetSim
#' @description This function calculates vertex proximity (similarity) with selected method
#' and between selected vertices.
#' @param graph an object of class \code{igraph}
#' @param method a method (single string) for calculating similarities , see details for methods.
#' @details Following methods are available:
#' \code{aa} Adamic-Adar index
#' \code{act} average commute time
#' \code{cn} common neighbours
#' \code{jc} cosine similarity
#' \code{dice} cosine similarity on L+
#' \code{invL} graph distance
#' \code{hdi} Hub Depressed Index
#' \code{hpi} Hub Promoted Index
#' \code{dist} Jaccard coefficient
#' \code{katz} Katz index
#' \code{lhn} Leicht-Holme-Newman Index
#' \code{lp} Local Path Index
#' \code{cosine} Matrix Forest Index
#' \code{pa} preferential attachment
#' \code{ra} resource allocation index
#' @return a matrix of similarity scores.
#' @references
#' \itemize{
#'   \item T. Zhou, L. Lv, and Y.-C. Zhang.Predicting missing links via local information, The European Physical Journal B - Condensed Matter and Complex Systems, vol. 71, no. 4, pp. 623-630, Oct. 2009.
#'   }
#' @examples
#' \donttest{
#' # net <- netSim(graph=proj_weighted,method="pa")
#' }
#' 
#' @export


unetSim<- function(graph, method){
    if (class(graph) != "igraph"){
        stop("The function applies to #igraph# object.\n")
    }
    if (igraph::is.directed(graph))
        stop("Graph has to be undirected")
    method <- match_method(method)
    print (method)
    result <- do.call(method, list(graph=graph))
    return (result)
}    



