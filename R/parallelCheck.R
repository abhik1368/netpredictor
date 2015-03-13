dCheckParallel <- function (multicores=NULL, verbose=T)
{
    
    # @import doMC
    # @import foreach
    
    if(0){
        flag_parallel <- F
        if (requireNamespace("doMC", quietly=T) & requireNamespace("foreach", quietly=T)) {
            if(1){
                doMC::registerDoMC()
                cores <- foreach::getDoParWorkers()
                if(is.null(multicores)){
                    multicores <- max(1, ceiling(cores*0.5))
                }else if(is.na(multicores)){
                    multicores <- max(1, ceiling(cores*0.5))
                }else if(multicores < 1 | multicores > cores){
                    multicores <- max(1, ceiling(cores*0.5))
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
    }
    
    if(1){
        flag_parallel <- F
        pkgs <- c("doMC","foreach")
        if(any(pkgs %in% rownames(installed.packages()))){
            tmp <- sapply(pkgs, function(pkg) {
                suppressPackageStartupMessages(require(pkg, character.only=T))
            })
            if(all(tmp)){
                doMC::registerDoMC()
                cores <- foreach::getDoParWorkers()
                print (cores)
                if(is.null(multicores)){
                    multicores <- max(1, ceiling(cores))
                }else if(is.na(multicores)){
                    multicores <- max(1, ceiling(cores))
                }else if(multicores < 1 | multicores > 2*cores){
                    multicores <- max(1, celing(cores))
                }else{
                    multicores <- as.integer(multicores)
                }
                doMC::registerDoMC(multicores) # register the multicore parallel backend with the 'foreach' package
                print (foreach::getDoParWorkers())
                if(verbose){
                    message(sprintf("\tdo parallel computation using %d cores ...", multicores, as.character(Sys.time())), appendLF=T)
                }
                flag_parallel <- T
            }
        }
    }
    
    return(flag_parallel)
}