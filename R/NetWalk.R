#' Random walk on unipartite networks
#' @title Perform Random walk on a Unipartite Network
#' @description Peforms random walk with restart with preferred seed sets. If seed sets are not given then the adjacencny
#' matrix is taken as the input as the input seed sets. THe restart parameter controls the random walk probability . This can be 
#' changed default is set to 0.8. Normalization of the matrix can be done by row,column,laplacian. For faster computation
#' Parallalization is implemented with multicores. Parallization is done using foreach package. 
#' @param ig : igraph object
#' @param normalise : normalise method 
#' @param setSeeds: vector or dataframe
#' @param restart: restart probability parameter
#' @param parallel: to execute in parallel either TRUE or FALSE
#' @param multcores: Number of cores to be used when running in parallel
#' @param Verbose: Verbose output
#' @name uNetwalk
#' @references  
#' \itemize{
#'   \item Kohler S, et al. Walking the Interactome for Prioritization of Candidate Disease Genes. American Journal of Human Genetics. 2008;82:949–958.
#'   \item Can, T., Çamoǧlu, O., and Singh, A.K. (2005). Analysis of protein-protein interaction networks using random walks. In BIOKDD '05: Proceedings of the 5th international workshop on Bioinformatics (New York, USA: Association for Computing Machinery). 61–68
#' }
#' @export
#' @examples
#' \donttest{
#' # Get source compound ids and source information
#' # Using ChEMBL ID and source
#' get.scid.sid("CHEMBL12",1)
#' # Using drugbank id and source
#' get.scid.sid("DB00789",2)
#' }

uNetwalk <- function(ig, normalise=c("row","column","laplacian","none"), setSeeds=NULL, restart=0.8, parallel=TRUE, multicores=NULL, verbose=T) {
    
    startT <- Sys.time()
    # Stopping criteria
    stop_delta <- 1e-7   
    if (class(ig) != "igraph"){
        stop("The function must apply to either 'igraph' or 'matrix' object.\n")
    
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=T)
    }
    
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
        c <- 0.8
    }
    else{
        c <- restart
    }
    normalise <- match.arg(normalise)
    #normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    

    if ("weight" %in% list.edge.attributes(ig)){
        adjM <- get.adjacency(g, type="both", attr="weight", edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using weighted graph!"), appendLF=T)
        }
    }else{
        adjM <- get.adjacency(g, type="both", attr=NULL, edges=F, names=T, sparse=getIgraphOpt("sparsematrices"))
        if(verbose){
            message(sprintf("\tNotes: using unweighted graph!"), appendLF=T)
        }
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Normalising the adjacency matrix using %s normalisation (%s) ...", normalise, as.character(now)), appendLF=T)
    }
  
    A <- adjM!=0
    if(normalise == "row"){
        D <- Matrix::Diagonal(x=(Matrix::rowSums(A))^(-1))
        nadjM <- adjM %*% D
    }else if(normalise == "column"){
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-1))
        nadjM <- D %*% adjM
    }else if(normalise == "laplacian"){
        D <- Matrix::Diagonal(x=(Matrix::colSums(A))^(-0.5))
        nadjM <- D %*% adjM %*% D
    }else{
        nadjM <- adjM
    }
    
    ##  A function to make elements in each steady probability vector is one column normalize
    colNorm<- function(m){
        #res <- t(t(m)/colSums(m))
        
        col_sum <- apply(PTmatrix, 2, sum)
        col_sum_matrix <- matrix(rep(col_sum, nrow(PTmatrix)), ncol=ncol(PTmatrix), nrow=nrow(PTmatrix), byrow =T)
        res <- as.matrix(PTmatrix)/col_sum_matrix
        res[is.na(res)] <- 0
        return(res)
    }
    
    if(is.null(setSeeds)){
        
        P0matrix <- Matrix::Matrix(diag(vcount(ig)), sparse=T)
        rownames(P0matrix) <- V(ig)$name
        colnames(P0matrix) <- V(ig)$name
        
    }else{
        ## check input data
        
        if(is.null(rownames(data))) {
            stop("The function must require the row names of the input setSeeds.\n")
        }else if(any(is.na(rownames(data)))){
            warning("setSeeds with NA as row names will be removed")
            data <- data[!is.na(rownames(data)),]
        }
        
        cnames <- colnames(data)
        if(is.null(cnames)){
            cnames <- seq(1,ncol(data))
        }
        
        if (class(ig) == "igraph"){
            ind <- match(rownames(data), V(ig)$name)
            nodes_mapped <- V(ig)$name[ind[!is.na(ind)]]
            if(length(nodes_mapped)!=vcount(ig)){
                warning("The row names of input setSeeds do not contain all those in the input graph.\n")
            }
            P0matrix <- matrix(0,nrow=nrow(nadjM),ncol=ncol(data))
            P0matrix[ind[!is.na(ind)],] <- as.matrix(data[!is.na(ind),])
        
            ## make sure the sum of elements in each steady probability vector is one
            P0matrix <- colNorm(P0matrix)
        
            ## Assign row and colnames 
            rownames(P0matrix) <- V(ig)$name
            colnames(P0matrix) <- cnames
        } 
    }    
        if(restart==1){
            ## just seeds themselves
            PTmatrix <- P0matrix
        }else{
            ###### Run in parallel
            flag_parallel <- F
            if(parallel==TRUE){
                
                flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
                if(flag_parallel){
                    j <- 1
                    PTmatrix <- foreach::`%dopar%` (foreach::foreach(j=1:ncol(P0matrix), .inorder=T, .combine="cbind"), {
                        P0 <- P0matrix[,j]
                        ## Initializing variables
                        delta <- 1
                        PT <- P0
                        ## Iterative update till convergence (delta<=1e-10)
                        while (delta>stop_delta){
                            PX <- (1-c) * nadjM %*% PT + c * P0
                            # p-norm of v: sum((abs(v).p)^(1/p))
                            delta <- sum(abs(PX-PT))
                            PT <- PX
                        }
                        as.matrix(PT)
                    })
                    
                    PTmatrix[PTmatrix<1e-6] <- 0
                    #PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
                }
            }
            if(flag_parallel==F){
                PTmatrix <- Matrix::Matrix(0, nrow=nrow(P0matrix), ncol=ncol(P0matrix), sparse=T)
                for(j in 1:ncol(P0matrix)){
                    #P0 <- as.matrix(P0matrix[,j],ncol=1)
                    P0 <- P0matrix[,j]
                    
                    ## Initializing variables
                    delta <- 1
                    
                    PT <- P0
                    ## Iterative update till convergence (delta<=1e-10)
                    while (delta>stop_delta){
                        PX <- (1-c) * nadjM %*% PT + c * P0
                        
                        # p-norm of v: sum((abs(v).p)^(1/p))
                        delta <- sum(abs(PX-PT))
                        
                        PT <- PX
                        step <- step+1
                    }
                    #PTmatrix[,j] <- as.matrix(PT, ncol=1)
                    PT[PT<1e-6] <- 0
                    PTmatrix[,j] <- Matrix::Matrix(PT, sparse=T)
                }
            }
        }
        if(verbose){
            now <- Sys.time()
            message(sprintf("Rescaling steady probability vector (%s) ...", as.character(now)), appendLF=T)
        }
        PTmatrix <- colNorm(PTmatrix) 
        PTmatrix <- Matrix::Matrix(PTmatrix, sparse=T)
        rownames(PTmatrix) <- rownames(P0matrix)
        colnames(PTmatrix) <- colnames(P0matrix)
        
        endT <- Sys.time()
        runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
        message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
        
        invisible(PTmatrix)
 
}    


#' Network based inference on Bipartite networks
#' @title Network Based Inference
#' @description Given a bipartite graph , a two phase resource transfer Information from  X(x,y,z) set of nodes gets distributed to Y set of nodes and then again goes back to resource X. 
#' This process allows us to define a technique for the calculation of the weight matrix W.
#' @name nbiNet
#' @param A: Adjacency Matrix
#' @param lamda: lamda parameter
#' @param alpha: alpha parameter
#' @param S: Target Similarity matrix
#' @param S1: Chemical Similarity Matrix
#' @param format: type of file as Adjacnency file
#' @name nbiNet
#' @references 
#' \itemize{
#' \item Cheng F, et al. Prediction of drug-target interactions and drug repositioning via network-based inference. PLoS Comput. Biol. 2012;8:e1002503.
#' \item Zhou T, et al. Solving the apparent diversity-accuracy dilemma of recommender systems. Proc. Natl Acad. Sci. USA 2010;107:4511-4515.
#' \item Zhou T, et al. Bipartite network projection and personal recommendation. Phys. Rev. E Stat. Nonlin. Soft Matter Phys. 2007;76:046115.
#' \item Blog post from Abhik Seal \url{http://data2quest.blogspot.com/2015/02/link-prediction-using-network-based.html}
#' }
#' @export

nbiNet <- function (A, lambda=0.5, alpha=0.5, S=NA, S1=NA,format = c("pairs","matrix")) {
    
    format <- match.arg(format)
    
    if (format == "matrix"){
        
        if(is.data.frame(A)){
            adjM <- as.matrix(A)
        }
    } else if(format == "pairs") {
        d<- graph.data.frame(A) ## only accepts pairs file
        V(d)$type <- V(d)$name %in% A[,1]
        data <-  get.incidence(d)
        adjM <- transpose(data)
    } else {
        stop("The function apply to either 'pairs file type' or 'matrix file type' object.\n")
    }
    
    n = nrow(adjM)
    m = ncol(adjM)
    if (nrow(S1) != m || ncol(S1) != m) {
        stop("The matrix S1 should be an m by m matrix with same number of columns as A.")
    }
    if (nrow(S) != n || ncol(S) != n) {
        stop("The matrix S should be an n by n matrix with same number of rows as A")
    }
    
    Ky <- diag(1/colSums(adjM))
    Ky[is.infinite(Ky) | is.na(Ky)] <- 0
    
    kx <- rowSums(adjM)
    Nx <- 1/(matrix(kx, nrow=n, ncol=n, byrow=TRUE)^(lambda) * 
                 matrix(kx, nrow=n, ncol=n, byrow=FALSE)^(1-lambda))
    Nx[is.infinite(Nx) | is.na(Nx)] <- 0 
    kx[is.infinite(kx) | is.na(kx)] <- 0 
    
    W <- t(adjM %*% Ky)
    W <-  adjM %*% W
    W <- Nx * W
    rownames(W) <- rownames(adjM)
    colnames(W) <- rownames(adjM)
    
    P1 <- adjM %*% S1 %*% t(adjM)
    P2 <- adjM %*% matrix(1, nrow=m, ncol=m) %*% t(adjM)
    S2 <- P1 / P2
    W  <- W * ((alpha * S) + ((1-alpha) * S2))
    
    W[is.nan(W)] <- 0
    rM <- W %*% adjM
    return (rM)
}


#' Randomm walk with restart on Bipartite networks
#' @title Bipartite Random Walk
#' @name biNetwalk
#' @param g1 :igraph object
#' @param s1: Similarity matrix for targets
#' @param s2: Similarity matrix for compounds
#' @param normalise: Normalisation of matrix using laplacian or None(the transition matrix will be column normalized)
#' @param setSeeds: seeds file
#' @param file: Accepts file for seeds
#' @param restart: restart value
#' @param parallel: parallel performance either True or False . Parallelization is implemented using foreach.
#' @param multicore: using multicores
#' @references 
#' \itemize {
#' \item Chen X, et al. Drug–target interaction prediction by random walk on the heterogeneous network. Mol. BioSyst 2012;8:1970-1978.
#' \item Vanunu O, Sharan R. Proceedings of the German Conference on Bioinformatics. Germany: GI; 2008. A propagation-based algorithm for inferring gene-disease assocations; pp. 54–63.
#' }
#' @export

biNetwalk <- function(g1,s1,s2,normalise=c("laplace","none"), setSeeds=NULL, file=NULL,restart=0.8, parallel=TRUE, multicores=NULL, verbose=T,weight=FALSE) {
    
    startT <- Sys.time()
    # Stopping criteria
    
    if (!exists('s1') || !exists('s2')){
        stop("You must submit s1 and s2 matrices.\n")
    }
    
    if (class(g1) != "igraph"){
        stop("The function applies to 'igraph' object.\n")
    }
    
    if (!bipartite.mapping(g1)$res){
        stop("The function applies to bipartite graphs.\n")
    }
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=T)
    }
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
        c <- 0.8
    }
    else{
        c <- restart
    }
    normalise <- match.arg(normalise)
    
    
    if ("weight" %in% list.edge.attributes(g1)){
        adjM <- get.incidence(g1, attr="weight", names=T)
        if(verbose){
            message(sprintf("Notes: using weighted graph!"), appendLF=T)
        }
    }else{
        adjM <- get.incidence(g1, attr=NULL, names=T)
        if(verbose){
            message(sprintf("Notes: using unweighted graph!"), appendLF=T)
        }
    }
    adjM <- as.matrix(adjM)
    # get the transition matrix
    W = tMat(adjM,s1,s2,normalise=normalise)
    
    if(is.null(setSeeds) && is.null(file)){
        
        M<-Matrix(adjM)
        M2<-0.99*M
        d<-Matrix(0.01*diag(nrow(s2)))
        P0matrix<-rBind(M2,d)
        
    }else if (is.null(file)){
        
        ## check input seeds
        if(is.vector(setSeeds)){           
            data <- as.matrix(setSeeds, ncol=1)
            rownames(data)<-setSeeds
        }
        
        cnames <- "V1"    
        ind <- match(rownames(data), rownames(W))
        nodes_mapped <- rownames(W)[ind[!is.na(ind)]]
        if(length(nodes_mapped)!=length(rownames(data))){
            warning("The row names of input setSeeds do not contain all those in the input graph.\n")
        }  
        P0matrix <- matrix(0,nrow=nrow(W),ncol=1)
        P0matrix[ind[!is.na(ind)],] <- 1 
        P0matrix <- colNorm(P0matrix)
    }else{
        
        # part of the section for input file name
        drug.names <- as.character(unique(file$V2))
        print (drug.names)
        P0matrix <- matrix(0,nrow=nrow(W),ncol=length(drug.names))
        for (i in 1:length(drug.names)){
            sub.fr <- file[file$V2==drug.names[i],]
            proteins <- as.character(sub.fr$V1)
            ind1 <- match(proteins, rownames(W))
            ind2 <- match(drug.names[i],rownames(W))
            ind <- append(ind1,ind2)
            nodes_mapped <- rownames(W)[ind[!is.na(ind)]]
            if(length(nodes_mapped)!=length(ind)){
                warning("The row names of input setSeeds do not contain all those in the input graph.\n")
            }
            
            P0matrix[ind[!is.na(ind)],i] <- 1 
        }
        P0matrix <- colNorm(P0matrix)
        
    }
    
    if (exists("W")){
        rmat <- rwr(W,P0matrix,par=parallel,r=c,multicores=multicores)
    } else{
        stop("Transition matrix couldnt be generated..")
    }
    
    if (!exists("rmat")){
        stop("Couldn't return the RWR matrix. \n")
    }else{
        if(verbose){
            now <- Sys.time()
            message(sprintf("Rescaling steady probability vector (%s) ...", as.character(now)), appendLF=T)
        }
        
        rmat[rmat < 1e-06] <- 0
        rmat <- rmat[1:nrow(adjM),]
        
        rmat <- colNorm(as.matrix(rmat))
        rownames(rmat)<- rownames(adjM)
        if(!is.null(file)){
            colnames(rmat)<- drug.names
            invisible(rmat)
        } else {
            colnames(rmat)<- colnames(adjM)
            invisible(rmat)
        }
        endT <- Sys.time()
        runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
        message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)    
        invisible(rmat)
        
    }
    
}

#' Significance of Bipartite networks 
#' @title Significant Network
#' @param data: n x m Adjancency matrix of seeds or dataframe of pairs.
#' @param g: igrah object of bipartite indcident adjacency matrix.
#' @param Amatrix: Affinity matrix computed from biNetWalk or uNetWalk.
#' @param num.permutation: number of permutation of Affinity matrix needed to performed.
#' @param adjp.cutoff: pvalue cutoff 0.05 
#' @param p.adjust.method: Adjusting the pvalue by diiferent method.It uses method from the stats package.  
#' @param parallel: Using parallization either True or False
#' @param multicores: If parallisation is set TRUE number of cores to perform parallel computaion.
#' @param Verbose: For verbose output.
#' @name Significant network
#' @export


sig.net <- function(data, g, Amatrix, num.permutation=10, adjp.cutoff=0.05, p.adjust.method=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr"), parallel=TRUE, multicores=NULL, verbose=T)
{   
    library(igraph)
    startT <- Sys.time()
    
    ####################################################################################
    permutation <- "random"
    p.adjust.method <- match.arg(p.adjust.method)
    ## check input data
    if(is.matrix(data) | is.data.frame(data)){
        
        data <- as.matrix(data)
        if (ncol(data)==2){
            g.data <- graph.data.frame(data) ## only accepts pairs file
            V(g.data)$type <- V(g.data)$name %in% data[,1]
        } else if (ncol(data) > 2){ ## Needs check seeds can be adjacency matrix
            data <- t(data)
            g.data <- graph.incidence(data)
        } else if(ncol(data) < 2){
            stop("The input data must be matrix or dataframe with at least two columns.\n")
        }
        
        
        
    }else if(is.null(data)){
        stop("The input data must be matrix.\n")
    }
    
    ## check input graph
    if (class(g) != "igraph"){
        stop("The function must apply to 'igraph' object.\n")
    }
    
    if (!bipartite.mapping(g.data)$res){
        stop("The function applies to igrah bipartite graphs.\n")
    }
    
    ## check bipartite mapping between input seed data and graph
    p <- bipartite.projection(g.data)
    rnames <- V(p[[1]])$name # rownames of the seed matrix
    cnames <- V(p[[2]])$name # colnames of the seed matrix
    
    # Get the adjacency matrix from igraph object
    g.incident <- get.incidence(g)
    
    ind <- match(cnames, rownames(g.incident))
    nodes_mapped <- rownames(Amatrix)[ind[!is.na(ind)]]
    P0matrix <- matrix(0, nrow=length(rownames(Amatrix)),ncol=length(rnames))
    P0matrix[ind[!is.na(ind)],] <- 1
    rownames(P0matrix) <- rownames(Amatrix)
    colnames(P0matrix) <- rnames
    
    ## check mapping between input Affinity matrix and graph
    ind1 <- match(rownames(Amatrix), rownames(g.incident))
    ind2 <- match(colnames(Amatrix), rnames)    
    if(length(ind1[!is.na(ind1)])!=length(rownames(g.incident)) && length(ind2[!is.na(ind2)])!=length(colnames(g.incident))) {
        stop("The function must require input Affinity matrix (Amatrix) has the same names (both columns and rows) as the input graph.\n")
    } 
    
    
    PTmatrix <- Amatrix[ind1[!is.na(ind1)],ind2]  
    PTmatrix <- colNorm(as.matrix(PTmatrix))
    
    
    ####################################################
    
    obs <- as.matrix(t(PTmatrix) %*% PTmatrix)
    B <- num.permutation
    if(verbose){
        message(sprintf("Third, generate the distribution of contact strength based on %d permutations on nodes respecting %s (%s)...", B, permutation, as.character(Sys.time())), appendLF=T)
    }
    
    
    f <- function(){
        pb <- txtProgressBar(min=1, max=num.permutation-1,style=3)
        count <- 0
        function(...) {
            count <<- count + length(list(...)) - 1
            setTxtProgressBar(pb,count)
            Sys.sleep(0.01)
            flush.console()
            c(...)
        }
    } 
    ###### parallel computing
    flag_parallel <- F
    if(parallel==TRUE){
        
        flag_parallel <- dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            b <- 1
            exp_b <- foreach::`%dopar%` (foreach::foreach(b=1:B, .inorder=T,.combine=f()), {
                PT_random <- PTmatrix[sample(nrow(PTmatrix)),sample(ncol(PTmatrix))]
                ## make sure the sum of elements in each steady probability vector is one
                PT_random <- colNorm(as.matrix(PT_random))
                as.matrix(t(as.matrix(PT_random)) %*% PT_random)
            })
        }
    }
    
    ## non-parallel computing
    if(flag_parallel==F){
        exp_b <- lapply(1:B, function(b){
            PT_random <- PTmatrix[sample(nrow(PTmatrix)),sample(ncol(PTmatrix))]
            ## make sure the sum of elements in each steady probability vector is one
            PT_random <- colNorm(as.matrix(PT_random))
            as.matrix(t(as.matrix(PT_random)) %*% PT_random)
        })
    }
    
    n <- ncol(obs)
    ## for zscore
    exp_mean <- matrix(0, ncol=n, nrow=n)
    exp_square <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        exp <- exp_b[[b]]
        exp_mean <- exp_mean + exp
        exp_square <- exp_square + exp^2
    }
    exp_mean <- exp_mean/B
    exp_square <- exp_square/B
    exp_std <- sqrt(exp_square-exp_mean^2)
    zscore <- (obs-exp_mean)/exp_std
    zscore[is.na(zscore)] <- 0
    zscore[is.infinite(zscore)] <- 0
    
    ## for pvalue
    num <- matrix(0, ncol=n, nrow=n)
    for(b in 1:B){
        num <- num + (obs < exp_b[[b]])
    }
    pval <- num/B
    colnames(pval) <- colnames(obs)
    rownames(pval) <- rownames(obs)
    
    ## for adjusted pvalue
    adjpval <- pval
    ## lower part
    flag_lower <- lower.tri(pval, diag=F)
    adjpval[flag_lower] <- stats::p.adjust(pval[flag_lower], method=p.adjust.method)
    ## upper part
    flag_upper <- upper.tri(pval, diag=F)
    adjpval[flag_upper] <- stats::p.adjust(pval[flag_upper], method=p.adjust.method)
    
    if(verbose){
        message(sprintf("Also, construct the contact graph under the cutoff %1.1e of adjusted-pvalue (%s)...", adjp.cutoff, as.character(Sys.time())), appendLF=T)
    }
    flag <- adjpval < adjp.cutoff
    adjmatrix <- flag
    adjmatrix[flag] <- zscore[flag]
    cgraph <- igraph::graph.adjacency(adjmatrix, mode="undirected", weighted=T, diag=F, add.colnames=NULL, add.rownames=NA)
    
    #########################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    result <- list(pval = pval,
                   adjpval = adjpval, 
                   cgraph  = cgraph)
    invisible(result)
}

#' Heterogenous Graph based Inference 
#' @title  Heterogenous Graph based Inference on Bipartite Network
#' @description Peforms random walk with restart with preferred seed sets. If seed sets are not given then the adjacencny
#' matrix is taken as the input as the input seed sets. THe restart parameter controls the random walk probability . This can be 
#' changed default is set to 0.8. Normalization of the matrix can be done by row,column,laplacian. For faster computation
#' Parallalization is implemented with multicores. Parallization is done using foreach package. 
#' @param ig : igraph object
#' @param normalise : normalise method 
#' @param setSeeds: vector or dataframe
#' @param restart: restart probability parameter
#' @param parallel: to execute in parallel either TRUE or FALSE
#' @param multcores: Number of cores to be used when running in parallel
#' @param Verbose: Verbose output
#' @name uNetwalk
#' @references  
#' \itemize{
#'   \item Kohler S, et al. Walking the Interactome for Prioritization of Candidate Disease Genes. American Journal of Human Genetics. 2008;82:949–958.
#'   \item Can, T., Çamoǧlu, O., and Singh, A.K. (2005). Analysis of protein-protein interaction networks using random walks. In BIOKDD '05: Proceedings of the 5th international workshop on Bioinformatics (New York, USA: Association for Computing Machinery). 61–68
#' }
#' @export

hgviNet <- function(g1,s1,s2,file=NULL,alpha=0.8,verbose=T) {
    
    startT <- Sys.time()
    # Stopping criteria
    
    if (!exists('s1') || !exists('s2')){
        stop("You must submit s1 and s2 matrices.\n")
    }
    
    if (class(g1) != "igraph"){
        stop("The function must apply to either 'igraph' object.\n")
    }
    
    if (!bipartite.mapping(g1)$res){
        stop("The function applies to bipartite graphs.\n")
    }
    
    
    
}

# 
# library(parallel)
# library(snow)
# cl <- makeCluster(detectCores())
# cl
# A <- matrix(rnorm(1000000), 1000)
# system.time(parMM(cl, A, A))
# system.time(A %*% A)

