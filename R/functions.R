## Performing jaccard similarity between two entities
## Performing jaccard similarity between two entities
## Using parallel to make cluster

jaccard.sim <- function(df){
    cl <- makeCluster(detectCores())
    a <- suppressWarnings(snow::parMM(cl,df,t(df)))
    b <- suppressWarnings(snow::parMM(cl,df,(1-t(df))))
    c <- suppressWarnings(snow::parMM(cl,(1-df),t(df)))

    sim <- a/(a + b + c)
    
    # changing diagonals to 1 and removing NaN's and NA
    diag(sim) <- 1
    s = as.matrix(sim)
    s[is.nan(s)] <- 0
    return (s)
    stopCluster(cl)
}

## Normalizing the matrix based on matrix columns

colNorm <- function(PTmatrix){
    
    if (ncol(PTmatrix) > 1){
        col_sum <- apply(PTmatrix, 2, sum)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
        res <- PTmatrix/col_sum_matrix
        res <- as.matrix(res)
        res[is.na(res)] <- 0
        return(res)
    }else {
        res <- PTmatrix[,1]/sum(PTmatrix[,1])
        return (as.matrix(res))
    }
    
    
}


## Get the transition matrix for a Bipartite Graph. Input a sequence similarity matrix (s1),
## chemical similarity matrix (s2) and drug target adjacency matrix (g1) where rows are 
## protein targets and columns as Drugs.
## @export

tMat <- function(g1,s1,s2,normalise="chen"){
    
    cl <- makeCluster(detectCores())
    
    g1 <- t(g1)
    seq<- as.matrix(s1)          ## sequence similairty matrix normalised between 0 and 1 
    drugProt <- as.matrix(g1)           ## drug target matrix   
    csim <- as.matrix(s2)         ## drug similarity matrix normalised between 0 and 1
    
    new.drug_drug <- suppressWarnings(snow::parMM(cl,drugProt,t(drugProt)))
    
    new.prot_prot <- suppressWarnings(snow::parMM(cl,t(drugProt),drugProt))

    #calculate drug-drug similarity based on shared proteins based on jaccard similarity            
    norm_drug <- jaccard.sim(drugProt)
    
    # Jaccard similarity of two proteins based on shared compounds
    norm_prot <- jaccard.sim(t(drugProt))
    
    # Normalizing the matrices with equal weights
    drug.similarity.final <- 0.5*(csim)+0.5*(norm_drug)
    prot.similarity.final <- 0.5*(seq)+0.5*(norm_prot)
    

    
    if(normalise == "laplace"){
        
        D1  <- diag(x=(rowSums(drugProt))^(-0.5))
        D2  <- diag(x=(colSums(drugProt))^(-0.5))   
        
        MTD1 <- suppressWarnings(snow::parMM(cl,D1,g1))
        MTD <- suppressWarnings(snow::parMM(cl,MTD1,D2))
        
        D3  <- diag(x=(rowSums(prot.similarity.final))^(-0.5))
        
        MTT1 <- suppressWarnings(snow::parMM(cl,D3,seq))
        MTT  <- suppressWarnings(snow::parMM(cl,MTT1,D3))
        
        D4  <- diag(x=(rowSums(drug.similarity.final))^(-0.5))
        MDD1  <- suppressWarnings(snow::parMM(cl,D4,csim))
        MDD  <- suppressWarnings(snow::parMM(cl,MDD1,D4))
        
        M1<-cbind(MTT,t(MTD))
        M2<-cbind(MTD,MDD)
        M <- rbind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        # Returning the final matrix 
        return(as.matrix(M))
    }    
    if(normalise == "chen"){
        
        ADT<-colSums(drugProt)
        ATD<-rowSums(drugProt)
        Sd<-rowSums(drug.similarity.final)
        St<-rowSums(prot.similarity.final)
        MDT<-mat.or.vec(nrow(drugProt),nrow(seq))
        MTD<-mat.or.vec(nrow(seq),nrow(drugProt))
        print (dim(MDT))
        print (dim(MTD))
        
        A<-t(drugProt)
        
        D1  <- diag(x=(rowSums(drugProt))^(-0.5))
        D2  <- diag(x=(colSums(drugProt))^(-0.5))   
        
        #MTD1 <- suppressWarnings(snow::parMM(cl,D1,g1))
        ##MTD <- suppressWarnings(snow::parMM(cl,MTD1,D2))
        #MDT <- t(MTD)
        D3  <- diag(x=(rowSums(prot.similarity.final))^(-0.5))
        
        MTT1 <- suppressWarnings(snow::parMM(cl,D3,seq))
        MTT  <- suppressWarnings(snow::parMM(cl,MTT1,D3))
        
        D4  <- diag(x=(rowSums(drug.similarity.final))^(-0.5))
        MDD1  <- suppressWarnings(snow::parMM(cl,D4,csim))
        MDD  <- suppressWarnings(snow::parMM(cl,MDD1,D4))
        
        for (i in 1:nrow(seq)){
            for (j in 1:nrow(seq)){
                if (ATD[i]==0){
                    
                    MTT[i,j]<-prot.similarity.final[i,j]/St[i]
                }
                else{
                    MTT[i,j]<-(0.8*(prot.similarity.final[i,j])/St[i])
                }
            }
        }
        
        for (i in 1:ncol(drugProt)){
            for (j in 1:ncol(drugProt)){
                if (ADT[i]==0){
                    MDD[i,j]<-drug.similarity.final[i,j]/Sd[i]
                }
                else{
                    MDD[i,j]<-(0.8*(drug.similarity.final[i,j])/Sd[i])
                }
            }
        }
        
        for (i in 1:ncol(drugProt)){
            for (j in 1:nrow(seq)){
                if (ADT[i]!=0){
                    MDT[i,j]<- (0.2*A[i,j])/ADT[i]
                }
                else{
                    MDT[i,j]<-0
                }
            }
        }
        for (i in 1:nrow(seq)){
            for (j in 1:ncol(drugProt)){
                if (ATD[i]!=0){
                    MTD[i,j]<- (0.2*A[j,i])/ATD[i]
                }
                else{
                    MTD[i,j]<-0
                }
            }
        }
        
        M1<-cbind(MTT,MTD)
        M2<-cbind(MDT,MDD)
        M <- rbind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        # Returning the final matrix 
        return(as.matrix(M))
        
    }

    
    if(normalise == "none"){
        MDD <- drug.similarity.final
        MTT <- prot.similarity.final
        MTD <- drugProt
        M1<-cbind(MTT,t(MTD))
        M2<-cbind(MTD,MDD)
        M <- rbind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        M <- colNorm(M)
        
        # Returning the final matrix 
        return(as.matrix(M))
    }
    on.exit(stopCluster(cl))
}

## performing random walk with restart both in parallel and non-parallel way.
## It takes the transition matrix W the initial matrix P0 matrix parameter for parallization
## and the number of cores to run parallelisation on. Also a restart parameter is available as r 
## to get better results with different datasets one needs to tune restart parameter r.

rwr <- function(W,P0matrix,r=0.9){
    
    stop_delta <- 1e-07     
    library(foreach)
    library(doParallel)
    registerDoParallel(detectCores())  
    message(paste(c("Executing RWR in parallel way .. \n")))
    PTmatrix <- matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix))
    pb <- txtProgressBar(min = 1, max = ncol(P0matrix), style = 3)
    #flag = dCheckParallel(multicores = 4,verbose = T)
    #if (flag){
    j = 1
    PTmatrix <- foreach::foreach(j=1:ncol(P0matrix), .inorder=T, .combine='cbind') %dopar% {
        setTxtProgressBar(pb, j)
        P0 <- P0matrix[,j]
        ## Initializing variables
        step <- 0
        PT <- P0
        ## Iterative update till convergence (delta<=1e-7)
        while (step<=stop_step){
            PX <- (1-r) * t(W) %*% PT + r * P0
            delta <- sum(abs(PX-PT))
            PT <- PX
            step <- step+1
        }
        as.matrix(PT)
    }
    
    
    return(PTmatrix)
    
}

#' Get the communities in a graph
#' @title Get communities from a graph
#' @name get.Communities
#' @description Get  the communities in a graph. It First decomposes the graph into components and then uses igraph's 
#' various community detection algorithms to detect communities in a graph. 
#' @param g igraph object
#' @param num.nodes Number of nodes to keep in the community.
#' @param calgo The community algorithm to use to find the communities.
#' @usage get.Communities(g,num.nodes = 2,calgo = walktrap.community)
#' @export 

get.Communities <- function(g,num.nodes = 2,calgo = walktrap.community){
    if (class(g)!= "igraph"){
        stop("The function must apply to 'igraph' object.\n")
    }

    ## decompose the graph 
    gps <- decompose.graph(g,min.vertices=num.nodes)
    # Create a list to store list of graphs for community identification
    result <- list()
    nc<- sapply(gps, vcount)
    # Get the indexes which has more than num.nodes
    indx <- which(nc > num.nodes)
    
    
    if (length(indx) < 1){
        stop("No communities found less than num.nodes")
    }
    #g = Z$cgraph
    community.significance.test <- function(graph, vs, ...) {
        if (is.directed(graph)) stop("This method requires an undirected graph")
        graph = g
        vs <- membership(comm)
        subgraph <- induced.subgraph(graph, vs)
        in.degrees <- degree(subgraph)
        out.degrees <- degree(graph, vs) - in.degrees
        wilcox.test(in.degrees, out.degrees, ...)
    }
    
    j <- 1
    for (i in indx){
        comm <- cluster_walktrap(gps[[i]],weights=E(gps[[i]])$weight)
        mem <- as.data.frame(as.matrix(membership(comm)))
        colnames(mem)[1] <- "Members"
        comm.sig <- community.significance.test(g,membership(comm))
        result[[j]] = list(community = comm,
                           community.sig = comm.sig,
                           cgraph = gps[[i]],
                           members = mem) 
        j <- j+1
        
    }
    class(result) <- "community"
    return (result)
}

#' plotting Communities 
#' @title plot_Community
#' @description This uses an object of getCommuntiy class and extracts parameters to plot current communities. 
#' @name plot_Community
#' @param gc receives the results from \code{get.Communities} 
#' @param cols number of columns to show in the plot
#' @examples
#' \dontrun{
#' ## Run the Bipartite Random walk with Restart
#' Q = biNetwalk(g1,S,S1,normalise="laplace",parallel=FALSE,verbose=T)
#' ## Get the significant vertices which are in columns
#' Z = sig.net(data=A,g=g1,Amatrix=Q,num.permutation=100,adjp.cutoff=0.01,p.adjust.method="BH",parallel=FALSE)
#' ## Get the graph for plotting commnuties
#' g <- Z$cgraph
#' gp <- get.Communities(g)
#' ## Total number of communities
#' length(gp)
#' ## Plot the communities with 5 columns
#' plot_Community(gp,cols=5)
#' }
#' @export

plot_Community <- function(gc,cols=3) {
    
    if (class(gc) != "community"){
        stop("The function applies to community object.\n")
    }

    numPlots = length(gc)
    r = ceiling(numPlots/cols)
    if (numPlots==1) {
        
        plot(gc[[1]]$community, gc[[1]]$cgraph)
        
    } else {
        
        par(mfrow = c(ceiling(numPlots/cols), cols))
        par(cex = 0.6)
        par(mar = c(4, 4, 2, 0.5), oma = c(2, 2, 2, 2))
    
        for (i in 1:numPlots){
            plot(gc[[i]]$community, gc[[i]]$cgraph,
                 main = paste(c("Community",i)),
                 layout=layout.fruchterman.reingold,
                 vertex.label.family="sans",
                 vertex.color= rainbow(10, .8, .8, alpha=.8))
            
        }
        
    }
    par(mfrow = c(ceiling(numPlots/cols), cols))
    mtext("Significant Communties", cex=1,font=2, col="Black", outer=TRUE)
}

#' Compute Degree Centrality of the Bipartite Graph
#' @title Degree centrality for Bipartite Graphs.
#' @name get.biDegreeCentrality
#' @description Measures graph degree centrality for bipartite graphs as well as Unipartite Graphs.The degree centrality of a vertex can be defined as fraction of incident edges over the number of all possible edges.
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param A Biparite matrix object
#' @param loops Either TRUE or FALSE keep loops or remove it.
#' @param SM Either TRUE or FALSE . If TRUE returns average centrality for rows and columns. If false returns degree centrality for all the nodes.
#' @return It returns a list of score of degree centrality of nodes at rows and nodes at columns.
#' @usage get.biDegreeCentrality(A,loops=FALSE,SM=TRUE)
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of the degree centrality the nodes in the rows are scored first followed by column vertices.
#' By deafult \code{SM=TRUE} which return the betweeness centrality for two different sets of nodes.
#' @export

get.biDegreeCentrality <-function(A,loops=FALSE,SM=TRUE){
        
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    if (SM==TRUE){
    if (!is.null(V(g)$type)){
        V(g)$degree <- degree(g)
        # determine maximum degrees for each vertex set
        max_d_TRUE <- max(V(g)[type==TRUE]$degree)
        max_d_FALSE <- max(V(g)[type==FALSE]$degree)
        # determine centralization for TRUE vertex subset
        g$sm.degree.centralization.columns <-sum(max_d_TRUE - V(g)[type==TRUE]$degree)/((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1))
        # determine centralization for FALSE vertex subset
        g$sm.degree.centralization.rows <-sum(max_d_FALSE - V(g)[type==FALSE]$degree)/((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1))
        # return both values as list
        return(list("SM.bi.degreeCentrality.columns"=g$sm.degree.centralization.columns,"SM.bi.degreeCentrality.rows"=g$sm.degree.centralization.rows))
    }
    else {
        # boolean vertex attribute 'type' is required
        warning("vertex attribute <type> is missing")
    }
   } else if (!is.null(V(g)$type)){
            
            for (i in V(g)){
                V(g)[i]$degreeCentrality <- degree(g,v=i)/length(V(g)[type==!V(g)[i]$type])
            }
            # return value vector as list item
            return(list("bi.degreeCentrality"=V(g)$degreeCentrality))
        }
        else {
            for (i in V(g)){
                if (!loops){
                    V(g)[i]$degreeCentrality <- degree(g,v=i,loops=FALSE)/(length(V(g))-1)
                }
                else{
                    V(g)[i]$degreeCentrality <- degree(g,v=i,loops=TRUE)/(length(V(g)))
                }
            }
            # return value vector as list item
            return(list("bi.degreeCentrality"=V(g)$degreeCentrality))
        }
}


#' Compute Graph density Bipartite Graph
#' @title Graph Density for Bipartite graphs.
#' @name get.biDensity
#' @description Measures graph density  of bipartite graphs . The density captures the fraction of actual present edges over the number of all possible edges, given that multiple edges are not allowed.  
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param A Biparite matrix object
#' @return It returns a list of score of degree centrality of nodes at rows and nodes at columns.
#' @usage get.biDensity(A)
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of the degree Density the nodes in the rows are scored first followed by column vertices.
#' @export

get.biDensity <-function(A){
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    if (!is.null(V(g)$type)){
        # return value as list item
        if (!is.directed(g)){
            return(list("bi.Density"=length(E(g))/(length(which(V(g)$type))*length(which(!V(g)$type)))))
        }
        else{
            return(list("bi.Density"=length(E(g))/(2*length(which(V(g)$type))*length(which(!V(g)$type)))))
        }
    }
    else {
        # boolean vertex attribute 'type' is required
        warning("vertex attribute <type> is missing")
    }
}

bipartite.closeness.centrality <- function(g){
    if (class(g) != "igraph"){
        stop("The function must apply to igraph object.\n")
        
    }
    if (!is.null(V(g)$type)){
        # determine maximal raw scores for both vertex subsets
        mrs_TRUE <- length(V(g)[type==FALSE]) + 2*length(V(g)[type==TRUE]) - 2
        mrs_FALSE <- length(V(g)[type==TRUE]) + 2*length(V(g)[type==FALSE]) - 2
        # get sum of all geodesic paths for each vertex
        sp <- shortest.paths(g)
        sp[sp==Inf] <- length(V(g))
        rowsums_shortest_paths <- rowSums(sp)
        # "bipartite" normalization of scores
        for (i in V(g)){
            if (V(g)[i]$type==TRUE){
                V(g)[i]$closeness.centrality <- mrs_TRUE/rowsums_shortest_paths[i+1]
            }
            else{
                V(g)[i]$closeness.centrality <- mrs_FALSE/rowsums_shortest_paths[i+1]
            }
        }
        # return value as list
        return(list("Bipartite.Closeness.Centrality"=V(g)$closeness.centrality))
    }
    else {
        # boolean vertex attribute 'type' is required
        cat("vertex attribute <type> is missing")
    }
}
#' Closeness centrality for bipartite graphs
#' @title Closeness centrality for bipartite graphs.
#' @name get.biClosenessCentralit
#' @description Measures Closeness centrality of bipartite graphs . The closeness centrality of a vertex is inversely proportional to the total geodesic distance to all other vertices.  
#' It also has way to calculate single mode betweeness centrality. Single mode centralization for bipartite graphs measures the extent to which vertices in one vertex subset are central relative only to other vertices in the same subset.
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param A Biparite matrix object
#' @param SM Either TRUE or FALSE . If TRUE returns average centrality for rows and columns. If false returns degree centrality for all the nodes.
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of closeness centrality of the nodes.By deafult \code{SM=TRUE} which return the closeness centrality for two different sets of nodes.
#' @export

get.biClosenessCentrality <-function(A,SM=TRUE){
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    if (SM==TRUE){
        
        if (!is.null(V(g)$type)){
            # determine bipartite closeness centrality scores
            d <- bipartite.closeness.centrality(g)[[1]]
            d[is.na(d)] <- 0
            V(g)$Bipartite.Closeness.Centrality <- d
            
            # get maximum scores
            max_c_TRUE <- max(V(g)[type==TRUE]$Bipartite.Closeness.Centrality)
            max_c_FALSE <- max(V(g)[type==FALSE]$Bipartite.Closeness.Centrality)
            
            # determine denominators for both vertex subsets
            if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
                denom_TRUE <- ((length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-2))/(2*length(V(g)[type==TRUE])-3) + ((length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-length(V(g)[type==FALSE])))/(length(V(g)[type==TRUE])+length(V(g)[type==FALSE])-2)
            }
            else{
                denom_TRUE <- ((length(V(g)[type==TRUE])-2)*(length(V(g)[type==TRUE])-1))/(2*length(V(g)[type==TRUE])-3)
            }
            if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
                denom_FALSE <- ((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-2))/(2*length(V(g)[type==FALSE])-3) + ((length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-length(V(g)[type==TRUE])))/(length(V(g)[type==FALSE])+length(V(g)[type==TRUE])-2)
            }
            else{
                denom_FALSE <- ((length(V(g)[type==FALSE])-2)*(length(V(g)[type==FALSE])-1))/(2*length(V(g)[type==FALSE])-3)
            }
            
            # determine centralization for TRUE vertex subset
            g$sm.closeness.centralization.columns <-sum(max_c_TRUE - V(g)[type==TRUE]$Bipartite.Closeness.Centrality)/denom_TRUE
            # determine centralization for FALSE vertex subset
            g$sm.closeness.centralization.rows <-sum(max_c_FALSE - V(g)[type==FALSE]$Bipartite.Closeness.Centrality)/denom_FALSE
            # return both values as list
            return(list("SM.bi.ClosenessCentrality.columns"=g$sm.closeness.centralization.columns,"SM.bi.ClosenessCentrality.rows"=g$single.mode.closeness.centralization.rows))
        }
        else {
            # boolean vertex attribute 'type' is required
            warning("vertex attribute <type> is missing")
        } 
        
    } else if (!is.null(V(g)$type)){
        # determine maximal raw scores for both vertex subsets
        mrs_TRUE <- length(V(g)[type==FALSE]) + 2*length(V(g)[type==TRUE]) - 2
        mrs_FALSE <- length(V(g)[type==TRUE]) + 2*length(V(g)[type==FALSE]) - 2
        # get sum of all geodesic paths for each vertex
        sp <- shortest.paths(g)
        sp[sp==Inf] <- length(V(g))
        rowsums_shortest_paths <- rowSums(sp)
        # "bipartite" normalization of scores
        for (i in V(g)){
            if (V(g)[i]$type==TRUE){
                V(g)[i]$closeness.centrality.columns <- mrs_TRUE/rowsums_shortest_paths[i+1]
            }
            else{
                V(g)[i]$closeness.centrality.rows <- mrs_FALSE/rowsums_shortest_paths[i+1]
            }
        }
        # return values as list
        return(list("bi.ClosenessCentrality.columns"=na.omit(V(g)$closeness.centrality.columns),"bi.ClosenessCentrality.rows"=na.omit(V(g)$closeness.centrality.rows)))
    } else {
        # boolean vertex attribute 'type' is required
        warning("vertex attribute <type> is missing")
    }
}
#' Betweeness centrality for bipartite graphs
#' @title Betweeness centrality for bipartite graphs.
#' @name get.biBetweenessCentrality
#' @description Measures Betweeness centrality of bipartite graphs. The betweenness centrality of a vertex  may be roughly defined as the number of geodesic paths that pass through a given vertex, 
#' weighted inversely by the total number of equivalent paths between the same two vertices, including those that do not pass through the given vertex.  
#' It also has way to calculate single mode betweeness centrality. Single mode centralization for bipartite graphs measures the extent to which vertices in one vertex subset are central 
#' relative only to other vertices in the same subset.This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param A Biparite matrix object
#' @param SM Either TRUE or FALSE . If TRUE returns average centrality for rows and columns. If false returns betweeness centrality for all the nodes.
#' @usage get.biBetweenessCentrality(A,SM=TRUE)
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of closeness centrality of the nodes. By deafult \code{SM=TRUE} which return the betweeness centrality for two different sets of nodes.
#' @export


get.biBetweenessCentrality <- function(A,SM=TRUE){
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    if (SM==TRUE){
        if (!is.null(V(g)$type)){
            # determine denominators for both vertex subsets
            if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
                denom_TRUE <- 2*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
            }
            else{
                denom_TRUE <- (length(V(g)[type==TRUE])-1) * ( 0.5*length(V(g)[type==FALSE])*(length(V(g)[type==FALSE])-1) + 0.5*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-2) + (length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1) )
            }
            if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
                denom_FALSE <- 2*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
            }
            else{
                denom_FALSE <- (length(V(g)[type==FALSE])-1) * ( 0.5*length(V(g)[type==TRUE])*(length(V(g)[type==TRUE])-1) + 0.5*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-2) + (length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1) )
            }
            
            # determine raw betweenness centrality scores from igraph
            V(g)$betweenness.raw <- betweenness(g,directed=FALSE)
            # get maximum scores
            max_c_TRUE <- max(V(g)[type==TRUE]$betweenness.raw)
            max_c_FALSE <- max(V(g)[type==FALSE]$betweenness.raw)
            
            # determine centralization for TRUE vertex subset
            g$sm.betweenness.centralization.columns <-sum(max_c_TRUE - V(g)[type==TRUE]$betweenness.raw)/denom_TRUE
            # determine centralization for FALSE vertex subset
            g$sm.betweenness.centralization.rows <-sum(max_c_FALSE - V(g)[type==FALSE]$betweenness.raw)/denom_FALSE
            # return both values as list
            return(list("SM.bi.BetweenessCentrality.columns"=g$sm.betweenness.centralization.columns,"SM.bi.BetweenessCentrality.rows"=g$sm.betweenness.centralization.rows))
        }
        else {
            # boolean vertex attribute 'type' is required
            warning("vertex attribute <type> is missing")
        }       
        
    } else if (!is.null(V(g)$type)){
        # determine maximal raw scores for both vertex subsets
        if (length(V(g)[type==FALSE])<length(V(g)[type==TRUE])){
            mrs_TRUE <- 2*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
        }
        else{
            mrs_TRUE <- 0.5*length(V(g)[type==FALSE])*(length(V(g)[type==FALSE])-1)+0.5*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==TRUE])-2)+(length(V(g)[type==FALSE])-1)*(length(V(g)[type==TRUE])-1)
        }
        if (length(V(g)[type==TRUE])<length(V(g)[type==FALSE])){
            mrs_FALSE <- 0.5*length(V(g)[type==TRUE])*(length(V(g)[type==TRUE])-1)+0.5*(length(V(g)[type==FALSE])-1)*(length(V(g)[type==FALSE])-2)+(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
            
        }
        else{
            mrs_FALSE <- 2*(length(V(g)[type==TRUE])-1)*(length(V(g)[type==FALSE])-1)
        }
        
        # get raw betweenness centrality scores from igraph
        betweenness_rs <- betweenness(g,directed=FALSE)
        # "bipartite" normalization of scores
        for (i in V(g)){
            if (V(g)[i]$type==TRUE){
                V(g)[i]$betweenness.centrality.columns <- betweenness_rs[i+1]/mrs_TRUE
            }
            else{
                V(g)[i]$betweenness.centrality.rows <- betweenness_rs[i+1]/mrs_FALSE
            }
        }
        # return value as list
        return(list("bi.BetweenessCentrality.columns"=na.omit(V(g)$betweenness.centrality.columns),"bi.BetweenessCentrality.rows"=na.omit(V(g)$betweenness.centrality.rows)))
    }
    else {
        # boolean vertex attribute 'type' is required
        warning("vertex attribute <type> is missing")
    }   
}

#' Projecting Bipartite Networks
#' @title Projecting Bipartite Networks.
#' @name get.biWeightedProjection
#' @description Projecting the Bipartite Network based on the vertex type (\code{TRUE} or \code{FALSE}).The mode \code{"shared-neighbours"} adds the number of shared neighbours in the bipartite graph to the edge linking two neighbours in the monopartite projection. 
#' \code{"newman"} adopts the weighting scheme presented by Newman (2001), that weighs the contributions of shared neighbours by the size of their linkage profiles minus one.  
#' @param A Biparite matrix object
#' @param vertex Either TRUE of FALSE where if FALSE then off-diagonal elements represent edge weights, diagonal elements vertex weights
#' @param mode \code{'shared-neighbours'} and \code{'newman'} 
#' @usage get.biWeightedProjection(A,vertex=FALSE,mode='shared-neighbours')
#' @references
#' \itemize{
#'   \item Newman M (2001) Scientific collaboration networks. II. Shortest paths, weighted networks, and centrality. Physical Review E 64.
#'   }
#' @return Returns a igraph object of the projected matrix with edge weights.
#' @export


get.biWeightedProjection <- function(A,vertex=FALSE,mode='shared-neighbours'){
    if (class(A) != "matrix"){
        stop("The function must apply to matrix object.\n")
        
    }
    g <- graph.incidence(A,mode = 'all')
    
    if (!bipartite.mapping(g)$res){
        stop("The function applies to bipartite graphs.\n")
    }

    if (is.bipartite(g)){
        V(g)$type <- bipartite.mapping(g)$type
        if (vertex) proj <- bipartite.projection(g)[[2]]
        else proj <- bipartite.projection(g)[[1]]
        V(proj)$name <- V(g)[type==vertex]$name
        ## number of shared vertices as edge weights 
        if (mode=='shared-neighbours'){
            incidenceMatrix <- get.incidence(g)	## TODO: might need some more testing if type and vType match
            if (!vertex){
                ## off-diagonal elements represent edge weights, diagonal elements vertex weights (currently unused)
                m <- incidenceMatrix %*% t(incidenceMatrix)
            }
            else {
                m <- t(incidenceMatrix) %*% incidenceMatrix
            }
            ## assign edge weights from matrix m
            if (length(E(proj))>0){
                for (i in 1:(length(E(proj)))){
                    E(proj)[i]$weight <- m[V(proj)[get.edge(proj,i)[[1]]]$name, V(proj)[get.edge(proj,i)[[2]]]$name]
                }
            }
            return(proj)
        }
        ## Newman2001 Formula
        if (mode=='newman'){
            incidenceMatrix <- get.incidence(g)
            if (!vertex){
                incidenceMatrix <- t(incidenceMatrix)
                if (length(c(which(rowSums(incidenceMatrix)==1)))!=0){
                    m2 <- incidenceMatrix[c(which(rowSums(incidenceMatrix)==1))*(-1),]
                }
                else{
                    m2 <- incidenceMatrix
                }
                m <- t(m2) %*% (m2*(1/(rowSums(m2)-1)))
            }
            else {	
                if (length(c(which(rowSums(incidenceMatrix)==1)))!=0){
                    m2 <- incidenceMatrix[c(which(rowSums(incidenceMatrix)==1))*(-1),]
                }
                else{
                    m2 <- incidenceMatrix
                }
                m <- t(m2) %*% (m2*(1/(rowSums(m2)-1)))
            }
            if (length(E(proj))>0){
                for (i in 1:(length(E(proj)))){
                    E(proj)[i]$weight <- m[V(proj)[get.edge(proj,i)[[1]]]$name, V(proj)[get.edge(proj,i)[[2]]]$name]
                }
            }
            return(proj)
        }
    }
}



