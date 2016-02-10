dCheckParallel <- function (multicores=NULL, verbose=T)
{
    
    # @import doMC
    # @import foreach
    
    flag_parallel <- F
    pkgs <- c("doMC","foreach")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            #suppressPackageStartupMessages(require(pkg, character.only=T))
            requireNamespace(pkg, quietly=T)
        })
        
        if(all(tmp)){
            doMC::registerDoMC()
            cores <- foreach::getDoParWorkers()
            if(is.null(multicores)){
                multicores <- max(1, ceiling(cores))
            }else if(is.na(multicores)){
                multicores <- max(1, ceiling(cores))
            }else if(multicores < 1 | multicores > 2*cores){
                multicores <- max(1, ceiling(cores))
            }else{
                multicores <- as.integer(multicores)
            }
            doMC::registerDoMC(multicores) # register the multicore parallel backend with the 'foreach' package
            
            if(verbose){
                message(sprintf("\tdo parallel computation using %d cores ...", multicores, as.character(Sys.time())), appendLF=T)
            }
            flag_parallel <- T
        }
        
    }
    
    return(flag_parallel)
}