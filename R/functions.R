
#' Performing jaccard similarity between two entities
jaccard.sim <- function(df){
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    
    sim <- a/(a + b + c)
    
    # changing diagonals to 1 and removing NaN's and NA
    diag(sim) <- 1
    s = as.matrix(sim)
    s[is.nan(s)] <- 0
    return (s)
}

#' Normalizing the matrix based on matrix columns

colNorm <- function(PTmatrix){
    if (ncol(PTmatrix) > 1){
        col_sum <- apply(PTmatrix, 2, sum)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
        res <- as.matrix(PTmatrix)/col_sum_matrix
        res[is.na(res)] <- 0
        return(res)
    }else {
        res <- PTmatrix[,1]/sum(PTmatrix[,1])
        return (as.matrix(res))
    }
    
    
}


#' Get the transition matrix for a Bipartite Graph. Input a sequence similarity matrix (s1),
#' chemical similarity matrix (s2) and drug target adjacency matrix (g1) where rows are 
#' protein targets and columns as Drugs.
#' @export

tMat <- function(g1,s1,s2,normalise="laplace"){
    g1 <- t(g1)
    seq<-(s1)          ## sequence similairty matrix normalised between 0 and 1 
    drugProt <- as.matrix(g1)           ## drug target matrix   
    csim<-(s2)         ## drug similarity matrix normalised between 0 and 1
    new.drug_drug<- drugProt %*% t(drugProt)    
    
    new.prot_prot<- t(drugProt) %*% drugProt

    #calculate drug-drug similarity based on shared proteins based on jaccard similarity            
    norm_drug <- jaccard.sim(drugProt)

    # Jaccard similarity of two proteins based on shared compounds
    norm_prot <- jaccard.sim(t(drugProt))
      
    # Normalizing the matrices with equal weights
    drug.similarity.final<- as.matrix(0.5*(csim)+0.5*(norm_drug))
    prot.similarity.final<- as.matrix(0.5*(seq)+0.5*(norm_prot))
    #print (prot.similarity.final)
    #print (rowSums(prot.similarity.final)) 
    if(normalise == "laplace"){
        D1  <- diag(x=(rowSums(drugProt))^(-0.5))
        D2  <- diag(x=(colSums(drugProt))^(-0.5))   
        MTD <- D1 %*% g1 %*% D2
        D3  <- diag(x=(rowSums(prot.similarity.final))^(-0.5))
        MTT <- D3 %*% seq %*% D3
        D4  <- diag(x=(rowSums(drug.similarity.final))^(-0.5))
        MDD  <- D4 %*% csim %*% D4 
        print (class(MTD))
        print (class(MTT))
        print (class(MDD))
        M1<-cbind(MTT,t(MTD))
        M2<-cbind(MTD,MDD)
        M <- rbind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        # Returning the final matrix 
        return(M)
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
        return(M)
    }
}

#' performing random walk with restart both in parallel and non-parallel way.
#' It takes the transition matrix W the initial matrix P0 matrix parameter for parallization
#' and the number of cores to run parallelisation on. Also a restart parameter is available as r 
#' to get better results with different datasets one needs to tune restart parameter r.

rwr <- function(W,P0matrix,par=FALSE,r=0.7,multicores=multicores){
    # run on sparse matrix package
    
    flag_parallel <- F
    stop_delta <- 1e-07     
        if(par==TRUE){
            flag_parallel <- dCheckParallel(multicores=multicores, verbose=T)
            message(paste(c("Executing parallel:")))
            print (flag_parallel)
            
            PTmatrix <- matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix))
            if(flag_parallel){
                j <- 1
                pb <- txtProgressBar(min = 1, max = ncol(P0matrix), style = 3)
                PTmatrix <- foreach::`%dopar%` (foreach::foreach(j=1:ncol(P0matrix), .inorder=T, .combine="cbind"), {
                    P0 <- P0matrix[,j]
                    ## Initializing variables
                    delta <- 1
                    PT <- P0
                    setTxtProgressBar(pb, j) 
                    ## Iterative update till convergence (delta<=1e-10)
                    while (delta>stop_delta){
                        PX <- (1-r) * t(W) %*% PT + r * P0
                        # p-norm of v: sum((abs(v).p)^(1/p))
                        delta <- sum(abs(PX-PT))
                        PT <- PX
                    }
                    as.matrix(PT)
                })
                
            }
            print (dim(PTmatrix))
           return(PTmatrix)
        } else if(flag_parallel==F){
            message(paste(c("Executing in non parallel way .. \n")))
            PTmatrix <- matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix))
            delta <- 1
            pb <- txtProgressBar(min = 1, max = ncol(P0matrix), style = 3)
            for(j in 1:ncol(P0matrix)){
                P0 <- P0matrix[,j]
                PT <- P0
                setTxtProgressBar(pb, j)
                while (delta > stop_delta){
                    PX <- (1-r) * t(W) %*% PT + r * P0
                    # p-norm of v: sum((abs(v).p)^(1/p))
                    delta <- sum(abs(PX-PT))
                    #print (delta)
                    PT <- PX
                }
                PTmatrix[,j] <- matrix(PT)
            }
            return(PTmatrix)
        }
}

#' Get the communities in a graph
#' @title Get communities from a graph
#' @description Get  the communities in a graph. It First decomposes the graph into components and then uses igraph's 
#' various community detection algorithms to detect communities in a graph. 
#' @param g: igraph object
#' @param num.nodes: Number of nodes to keep in the community.
#' @param calgo: The community algorithm to use to find the communities.
#' @name get.Communities
#' @export 

get.Communities<- function(g,num.nodes = 3,calgo = walktrap.community){
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
    
    community.significance.test <- function(graph, vs, ...) {
        if (is.directed(graph)) stop("This method requires an undirected graph")
        subgraph <- induced.subgraph(graph, vs)
        in.degrees <- degree(subgraph)
        out.degrees <- degree(graph, vs) - in.degrees
        wilcox.test(in.degrees, out.degrees, ...)
    }
    
    j <- 1
    for (i in indx){
         
        comm <- calgo(gps[[i]],weights=E(gps[[i]])$weight)
        mem <- data.frame(membership(comm))
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
#' @name plot.Community
#' @description This uses an object of getCommuntiy class and extracts parameters to plot current communities. 
#' @export

plot.Community <- function(gc,cols=3) {
    
    if (class(gc) != "community"){
        stop("The function applies to community object.\n")
    }

    numPlots = length(gc)
    r = ceiling(numPlots/cols)
    if (numPlots==1) {
        
        plot(gc[[1]]$community, gc[[1]]$cgraph)
        
    } else {
#         # Set up the page
#         grid.newpage()
#         pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#         
#         # Make each plot, in the correct location
#         for (i in 1:numPlots) {
#             # Get the i,j matrix positions of the regions that contain this subplot
#             matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#             #p <- plot(gc$community[[i]], gc$cgraph[[i]])
#             print(plot(gc[[i]]$community, gc[[i]]$cgraph), vp = viewport(layout.pos.row = matchidx$row,
#                                             layout.pos.col = matchidx$col))
        
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
#' @description Measures graph degree centrality for bipartite graphs as well as Unipartite Graphs.The degree centrality of a vertex can be defined as fraction of incident edges over the number of all possible edges.
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param g: igraph object
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of the degree centrality the nodes in the rows are scored first followed by column vertices.

degreeCentrality <-function(g,loops=FALSE){
    
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
        
    } 
    # if boolean vertex attribute <type> is present, calculate bipartite degree centrality, otherwise monopartite degree centrality
    if (!is.null(V(g)$type)){
        for (i in V(g)){
            V(g)[i]$degreeCentrality <- degree(g,v=i)/length(V(g)[type==!V(g)[i]$type])
        }
        # return value vector as list item
        return(list("Bipartite.Degree.Centrality"=V(g)$degreeCentrality))
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
        return(list("Unipartite.Degree.Centrality"=V(g)$degreeCentrality))
    }
}


#' Compute Graph density Bipartite Graph
#' @title Graph Density for Bipartite graphs.
#' @description Measures graph density  of bipartite graphs . The density captures the fraction of actual present edges over the number of all possible edges, given that multiple edges are not allowed.  
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param g: igraph object
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of the degree centrality the nodes in the rows are scored first followed by column vertices.

bgraphDensity <-function(g){
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
        
    } 
    if (!is.null(V(g)$type)){
        # return value as list item
        if (!is.directed(g)){
            return(list("Density"=length(E(g))/(length(which(V(g)$type))*length(which(!V(g)$type)))))
        }
        else{
            return(list("Density"=length(E(g))/(2*length(which(V(g)$type))*length(which(!V(g)$type)))))
        }
    }
    else {
        # boolean vertex attribute 'type' is required
        warning("vertex attribute <type> is missing")
    }
}
#' Closeness centrality for bipartite graphs
#' @title Closeness centrality for bipartite graphs.
#' @description Measures Closeness centrality of bipartite graphs . TThe closeness centrality of a vertex is inversely proportional to the total geodesic distance to all other vertices.  
#' This function is the adaption to bipartite graphs as presented in Borgatti and Everett (1997).
#' @param g: igraph object
#' @references
#' \itemize{
#'   \item Borgatti, S. P. and Everett, M. G. (1997) Network analysis of 2--mode data. Social Networks \bold{19}, 243--269.
#'   }
#' @return Returns a list scores of closeness centrality of the nodes.
#' 

bcloseCentrality <-function(g){
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
        
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
        warning("vertex attribute <type> is missing")
    }
}

bBetweenCentrality <- function(g){
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
        
    }     
    if (!is.null(V(g)$type)){
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
                V(g)[i]$betweenness.centrality <- betweenness_rs[i+1]/mrs_TRUE
            }
            else{
                V(g)[i]$betweenness.centrality <- betweenness_rs[i+1]/mrs_FALSE
            }
        }
        # return value as list
        return(list("Bipartite.Betweenness.Centrality"=V(g)$betweenness.centrality))
    }
    else {
        # boolean vertex attribute 'type' is required
        cat("vertex attribute <type> is missing")
    }
}







getDrugbanksdf <- function (id, parallel = 2) {
    
    if (is.null(id)){
        stop("id parameter is null")
    }
    
    URL = paste0('http://www.drugbank.ca/structures/structures/small_molecule_drugs/', id, '.sdf')   
    SDF = getURLAsynchronous(url = URL, perform = parallel)   
    return(SDF)
    
}

getDrugBankSmi <- function (id, parallel = 2) {
    
    if (is.null(id)){
        stop("id parameter is null")
    }
    id = c("DB00295","DB00291","DB00545")
    parallel = 2
    URL = paste0('http://www.drugbank.ca/structures/structures/small_molecule_drugs/', id, '.smiles')    
    SMILES = getURLAsynchronous(url = URL, perform = parallel)
    smiles <- gsub("[\r\n]", "", SMILES)
    return(smiles)
    
}


getPubchemSmi <- function (id, parallel = 5) {
    
    if (is.null(id)){
        stop("id parameter is null")
    }
    
    URL = sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES/TXT', id)    
    SMILES = getURLAsynchronous(url = URL, perform = parallel)    
    smiles <- gsub("[\r\n]", "", SMILES)
    return(smles)
    
}

getPubchemSdf <- function(id, parallel = 5){
    if (is.null(id)){
        stop("id parameter is null")
    }
    
    SdfURL = paste0('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=', 
                    id, '&disopt=DisplaySDF')
    
    SdfTxt = getURLAsynchronous(url = SdfURL, perform = parallel)
    
    return(SdfTxt)
}


getChEMBLSmi <- function(id, parallel = 5){
    if (is.null(id)){
        stop("id parameter is null")
    }
    URL = sprintf('https://www.ebi.ac.uk/chemblws/compounds/%s.json', id)    
    jsonTxt = getURLAsynchronous(url = URL, perform = parallel)
    data <- fromJSON(jsonTxt)
    smiles <- data$compound$smiles
    smiles <- gsub("[\r\n]", "", SMILES)
    return(smiles)
}


getChEMBLSdf = function (id, parallel = 5) {
    
    URL = paste0('https://www.ebi.ac.uk/chembldb/compound/inspect/', id)
    MolPageTxt = getURLAsynchronous(url = URL, perform = parallel)
    
    n = length(id)
    tmp1 = rep(NA, n)
    tmp2 = rep(NA, n)
    
    for (i in 1:n) {
        tmp1[i] = strsplit(MolPageTxt, 
                           "<a href='/chembldb/download_helper/getmol/")[[1]][2]
    }
    
    for (i in 1:n) {
        tmp2[i] = strsplit(tmp1[i], "'>Download MolFile")[[1]][1]
    }
    
    MolURL = paste0('https://www.ebi.ac.uk/chembldb/download_helper/getmol/', tmp2)
    sdf = getURLAsynchronous(url = MolURL, perform = parallel)
    
    return(sdf)
    
}
   
### More function related to proteins needs to be added
getSeqUniProt = function (id, parallel = 5) {

    fastaURL = paste0('http://www.uniprot.org/uniprot/', id, '.fasta')    
    fastaTxt = getURLAsynchronous(url = fastaURL, perform = parallel)    
    return(fastaTxt)
    
}


