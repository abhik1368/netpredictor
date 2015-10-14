##### Install the library
library(netpredictor)
library(igraph)


#### Load the data sets to perform prediction and analysis.
## Enzyme data

EN_ADJ <- read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_admat_dgc.txt",sep="\t",row.names = 1,header=TRUE)
EN_CSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_simmat_dc.txt",sep="\t",row.names=1,header=TRUE))
EN_PSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_simmat_dg.txt",sep='\t',row.names = 1,header = TRUE))

## Ion Channel
IC_ADJ <- read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_admat_dgc.txt",sep="\t",row.names = 1,header=TRUE)
IC_CSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_simmat_dc.txt",sep="\t",row.names = 1,header=TRUE))
IC_PSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_simmat_dg.txt",sep="\t",row.names=1,header=TRUE))


## GPCR
GPCR_ADJ <- read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_admat_dgc.txt",sep = '\t',row.names = 1,header=TRUE)
GPCR_CSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_simmat_dc.txt",sep='\t',row.names = 1,header=TRUE))
GPCR_PSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_simmat_dg.txt",sep='\t',row.names = 1,header=TRUE))

## NUCLEAR RECEPTOR
NR_ADJ <- read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_admat_dgc.txt",sep = '\t',row.names = 1,header=TRUE)
NR_CSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_simmat_dc.txt",sep = '\t',row.names = 1,header=TRUE))
NR_PSIM <- as.matrix(read.csv("http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_simmat_dg.txt",sep="\t",row.names = 1,header=TRUE))


## Performance check 
EN_ADJ <- t(EN_ADJ)
IC_ADJ <- t(IC_ADJ)
GPCR_ADJ <- t(GPCR_ADJ)
NR_ADJ <- t(NR_ADJ)
  
EN_PRED   <- nbiNet(EN_ADJ,alpha=0.5, lamda=0.5,  s1=EN_CSIM, s2=EN_PSIM,format = "matrix")  
IC_PRED   <- nbiNet(IC_ADJ,alpha=0.5, lamda=0.5,  s1=IC_CSIM, s2=IC_PSIM,format = "matrix")
GPCR_PRED <- nbiNet(GPCR_ADJ,alpha=0.5, lamda=0.5,  s1=GPCR_CSIM, s2=GPCR_PSIM,format = "matrix")  
NR_PRED   <- nbiNet(NR_ADJ,alpha=0.5, lamda=0.5,  s1=NR_CSIM, s2=NR_PSIM,format = "matrix")  




## Set a random seed for permutation analysis
set.seed(12345)
EN_PERM = list()
IC_PERM = list()
GPCR_PERM = list()
NR_PERM = list()
## Compute scores for 1000 permutations where you sample the matrix everytime

for ( i in 1:1000){
    
    message(sprintf("Computing for %s permuation" ,as.character(i)))

    EN_ADJ <- EN_ADJ[sample(nrow(EN_ADJ)),sample(ncol(EN_ADJ))]
    IC_ADJ <- IC_ADJ[sample(nrow(IC_ADJ)),sample(ncol(IC_ADJ))]
    GPCR_ADJ <- GPCR_ADJ[sample(nrow(GPCR_ADJ)),sample(ncol(GPCR_ADJ))]
    NR_ADJ <- NR_ADJ[sample(nrow(NR_ADJ)),sample(ncol(NR_ADJ))]
    
    
    EN_CSIM <- EN_CSIM[sample(nrow(EN_CSIM)),sample(ncol(EN_CSIM))]
    EN_PSIM <- EN_PSIM[sample(nrow(EN_PSIM)),sample(ncol(EN_PSIM))]
    
    IC_CSIM <- EN_CSIM[sample(nrow(IC_CSIM)),sample(ncol(IC_CSIM))]
    IC_PSIM <- EN_PSIM[sample(nrow(IC_PSIM)),sample(ncol(IC_PSIM))]

    GPCR_CSIM <- EN_CSIM[sample(nrow(GPCR_CSIM)),sample(ncol(GPCR_CSIM))]
    GPCR_PSIM <- EN_PSIM[sample(nrow(GPCR_PSIM)),sample(ncol(GPCR_PSIM))]
    

    NR_CSIM <- EN_CSIM[sample(nrow(NR_CSIM)),sample(ncol(NR_CSIM))]
    NR_PSIM <- EN_PSIM[sample(nrow(NR_PSIM)),sample(ncol(NR_PSIM))]
    
    
    EN_RESULT <- nbiNet(EN_ADJ,alpha=0.5, lamda=0.5, s1=EN_CSIM, s2=EN_PSIM,format = "matrix")  
    IC_RESULT <- nbiNet(IC_ADJ,alpha=0.5, lamda=0.5, s1=IC_CSIM, s2=IC_PSIM,format = "matrix")  
    GPCR_RESULT <- nbiNet(GPCR_ADJ,alpha=0.5, lamda=0.5,s1=GPCR_CSIM, s2=GPCR_PSIM,format = "matrix") 
    NR_RESULT <- nbiNet(NR_ADJ,alpha=0.5, lamda=0.5, s1=NR_CSIM, s2=NR_PSIM,format = "matrix")  
    
    EN_PERM[[i]] <- EN_RESULT
    IC_PERM[[i]] <- IC_RESULT
    GPCR_PERM[[i]] <- GPCR_RESULT
    NR_PERM[[i]] <- NR_RESULT
    
}

## Get the mean and standard deviation for each of the 1000 permuted matrix 

EN_PERM_MEAN <- apply(simplify2array(EN_PERM), 1:2, mean)
EN_PERM_SD <-  apply(simplify2array(EN_PERM), 1:2, sd)

IC_PERM_MEAN <- apply(simplify2array(IC_PERM), 1:2, mean)
IC_PERM_SD <-  apply(simplify2array(IC_PERM), 1:2, sd)

GPCR_PERM_MEAN <- apply(simplify2array(GPCR_PERM), 1:2, mean)
GPCR_PERM_SD <-  apply(simplify2array(GPCR_PERM), 1:2, sd)

NR_PERM_MEAN <- apply(simplify2array(NR_PERM), 1:2, mean)
NR_PERM_SD <-  apply(simplify2array(NR_PERM), 1:2, sd)


## Comput the Z-scores of matrix
EN_ZSCORE <- (EN_PRED  - EN_PERM_MEAN)/EN_PERM_SD
IC_ZSCORE <- (IC_PRED  - IC_PERM_MEAN)/IC_PERM_SD
GPCR_ZSCORE <- (GPCR_PRED  - GPCR_PERM_MEAN)/GPCR_PERM_SD
NR_ZSCORE <- (NR_PRED  - NR_PERM_MEAN)/NR_PERM_SD


## Compute the significance score for each of the different datasets
EN_SIG_NET<- 2*(pnorm(-abs(EN_ZSCORE)))
IC_SIG_NET<- 2*(pnorm(-abs(IC_ZSCORE)))
GPCR_SIG_NET<- 2*(pnorm(-abs(GPCR_ZSCORE)))
NR_SIG_NET<- 2*(pnorm(-abs(NR_ZSCORE)))



## Get the significant interactions where, P < 0.05 for each of the datasets

EN_SIG_NET[EN_SIG_NET <  0.05] <- 1
EN_SIG_NET[EN_SIG_NET != 1] <- 0

IC_SIG_NET[IC_SIG_NET <  0.05] <- 1
IC_SIG_NET[IC_SIG_NET != 1] <- 0

GPCR_SIG_NET[GPCR_SIG_NET <  0.05] <- 1
GPCR_SIG_NET[GPCR_SIG_NET != 1] <- 0

NR_SIG_NET[NR_SIG_NET <  0.05] <- 1
NR_SIG_NET[NR_SIG_NET != 1] <- 0


getNetwork <- function(pred,orig){
    
    ## Modify the edge type
    pdmat <- pred - orig
    edges <- as.data.frame(get.edgelist(graph.incidence(pdmat)))
    edges <- edges[order(edges$V1, edges$V2), ]
    edges$type <- "predicted"
    O_GRAPH <- as.data.frame(get.edgelist(graph.incidence(orig)))
    O_GRAPH  <- O_GRAPH [order(O_GRAPH $V1, O_GRAPH $V2), ]
    O_GRAPH$type <- "True"
    out <- rbind(edges,O_GRAPH)
    g <- graph.data.frame(out)
    g <- set.edge.attribute(g, "type", value=out$type)
    
    ## set tge vertex color 
    count_drugs <- length(unique(out$V1))
    count_target <- length(unique(out$V2))    
    drug_color <- rep("tomato", count_drugs)
    target_color <- rep("gold",count_target)
    V(g)$color<- c(drug_color,target_color)
    
    ## Removes loops and duplicated edges
    g <- simplify(g, remove.multiple=TRUE,remove.loops=TRUE)
    deg <- degree(g, mode="all")
    V(g)$size <- deg*0.15
    return (g)
}

g <- getNetwork(pred=EN_SIG_NET,orig=EN_ADJ)

V(g)$label.cex = 0.5
plot(g,edge.arrow.mode=0,layout=layout.fruchterman.reingold)

saveGML(g,"EN.gml","EN")


Enzyme <- melt(EN_SIG_NET)
IonChannel <- melt(IC_SIG_NET)
GPCR <- melt(GPCR_SIG_NET)
NR <- melt(NR_SIG_NET)

write.csv(Enzyme,"Enzyme_Results.csv")
write.csv(IonChannel,"IC_Results.csv")
write.csv(GPCR,"GPCR_Results.csv")
write.csv(NR,"NR_Results.csv")
