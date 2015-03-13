#' jaccard.sim
#' @param
#' 


library(utils)
library(igraph)


jaccard.sim <- function(df){
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    
    sim <- a/(a + b + c)
    
    # changing diagonals to 1 and removing NaN's and NA
    diag(sim) <- 1
    s = as.matrix(sim)
    s[is.nan(s)] <- 0
    return (Matrix(s))
}

#' colNorm
#' @param
#' 

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

#' t.mat
#' @param
#' @param
#' @param
#' @param
#' @description


t.mat<-function(g1,s1,s2,normalise){
    g1 <- t(g1)
    seq<-s1          ##sequence similairty matrix normalised between 0 and 1 
    drugProt<-g1     ##drug target matrix   
    csim<-s2         ##drug similarity matrix normalised between 0 and 1
        
    new.drug_drug<-as.matrix(crossprod(drugProt))    
    
    new.prot_prot<-as.matrix(tcrossprod(drugProt))

    #calculate drug-drug similarity based on shared proteins based on jaccard similarity            
    norm_drug <- jaccard.sim(drugProt)

    # Jaccard similarity of two proteins based on shared compounds
    norm_prot <- jaccard.sim(t(drugProt))
    # Normalizing the matrices with equal weights
    drug.similarity.final<-0.5*(csim)+0.5*(norm_drug)
    prot.similarity.final<-0.5*(seq)+0.5*(norm_prot)    
 
    
    if(normalise == "laplace"){
        D1  <- Diagonal(x=(rowSums(drugProt))^(-0.5))
        D2  <- Diagonal(x=(colSums(drugProt))^(-0.5))   
        MTD <- D1 %*% g1 %*% D2
        D3  <- Diagonal(x=(rowSums(prot.similarity.final))^(-0.5))
        MTT <- D3 %*% s1 %*% D3
        D4  <- Diagonal(x=(rowSums(drug.similarity.final))^(-0.5))
        MDD  <- D4 %*% s2 %*% D4 
        M1<-cBind(MTT,t(MTD))
        M2<-cBind(MTD,MDD)
        M <- rBind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        # Returning the final matrix 
        return(Matrix(M))
    }
    
    
    if(normalise == "none"){
        MDD <- drug.similarity.final
        MTT <- prot.similarity.final
        MTD <- drugProt
        M1<-cBind(MTT,t(MTD))
        M2<-cBind(MTD,MDD)
        M <- rBind(M1,M2)
        M <- as.matrix(M)
        M[is.na(M)]<-0
        n =c(colnames(g1),rownames(g1))
        rownames(M) <- n
        colnames(M) <- n
        M <- colNorm(M)
        # Returning the final matrix 
        return(Matrix(M))
    }
}

#' rwr
#' @param
#' @param
#' @param
#' @param
#' @description

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
