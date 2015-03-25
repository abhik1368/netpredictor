#' get.candidate.graph
#' @title get.candidata.graph
#' @name get.candidate.graph
#' @param file: Input seeds file must be 2 column dataframe or matrix of seed nodes of a bipartite graph. 
#' @param affinity: Computed affinity matrix from Bipartite network.
#' @param top: retrieving top number of results.
#' @export

get.candidate.graph<- function(file=NULL,affinity=NULL,top=10,format=c("bipart","unipart")){
    
    ## function for using 'not in' in R 
    `%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
    
    if (match.arg(format)=="unipart"){
        if(is.null(file)){
            stop("A file is required..")
        }
        if (!exists("affinity")){
            stop("Affinity data/matrix is required.. ")
        }
        cnames <- colnames(file)
        pairs <- data.frame()
        for (cname in cnames){
            pnames <- rownames(file[file$cname==1,])
            ind <- which(rownames(affinity) %not in% pnames)
            candidates <- as.matrix(affinity[ind[!is.na(ind)],])
            sortedlist <- candidates[order(candidates[,1],decreasing=TRUE),]
            result <- as.matrix(sortedlist[1:top])
            pp <- cbind(result,cname)
            pairs <- rbind(pairs,pp)
        }
        print (pairs)
        return(pairs)
        
    }else if (match.arg(format)=="bipart"){
        if(is.null(file)){
            stop("A pair file is required.. \n")
        }
        if (!exists("affinity")){
            stop("Affinity data is required use bi.netwalk/netwalk to get it..\n")
        }
        
        drug.names <- as.character(unique(file$V2))
        pairs <- data.frame()
        
        for (i in 1:length(drug.names)){
            sub.fr <- file[file$V2==drug.names[i],]
            proteins <- as.character(sub.fr$V1)
            ind <- which(rownames(affinity) %not in% proteins)
            candidates <- as.matrix(affinity[,i][ind[!is.na(ind)]])
            sortedlist <- candidates[order(candidates[,1],decreasing=TRUE),]
            result <- as.matrix(sortedlist[1:top])
            pp <- cbind(result,drug.names[i])
            pairs <- rbind(pairs,pp)
        }
        
        colnames(pairs)[1]<-"weight"
        colnames(pairs)[2]<-"source"
        pairs$target <- rownames(PP)
        return(pairs)
   }
}

#' Export file to GML
#' @title save GML 
#' @description Export an igraph object to GML file which can be in Cytoscape or Gephi.
#' @param g: igraph object
#' @param fileName: exported file name
#' @param title: title of the exported file name
#' @export  

saveGML = function(g, fileName, title = "untitled") {
    attrToString = function(x) {
        m = mode(x)
        if(m == "numeric") {
            xc = sprintf("%.12f", x)
            
            xc[is.na(x)] = "NaN"
            xc[x == "Infinity"]= "Infinity"
            xc[x == "-Infinity"]= "-Infinity"
            x = xc
        } else {
            #Escape invalid characters
            x = gsub('"', "'", x)
            x = paste("\"", x , "\"", sep="")
        }
        x
    }
    
    illAttrChar = "[ -]"
    vAttrNames = list.vertex.attributes(g)
    vAttrNames = vAttrNames[vAttrNames != "id"]
    vAttr = lapply(vAttrNames, function(x) attrToString(get.vertex.attribute(g, x)))
    names(vAttr) = gsub(illAttrChar, "", vAttrNames)
    eAttrNames = list.edge.attributes(g)
    eAttrNames = eAttrNames[eAttrNames != "id"]
    eAttr = lapply(eAttrNames, function(x) attrToString(get.edge.attribute(g, x)))
    names(eAttr) = gsub(illAttrChar, "", eAttrNames)
    
    f = file(fileName, "w")
    cat("graph\n[", file=f)
    cat(" directed ", as.integer(is.directed(g)), "\n", file=f)
    for(i in seq_len(vcount(g))) {
        cat(" node\n [\n", file=f)
        cat("    id", i, "\n", file=f)
        for(n in names(vAttr)) {
            cat("   ", gsub("[\\._]", "", n), vAttr[[n]][i], "\n", file=f)
        }
        cat(" ]\n", file=f)
    }
    
    el = get.edgelist(g, names=FALSE)
    for (i in seq_len(nrow(el))) { 
        cat(" edge\n  [\n", file=f) 
        cat("  source", el[i,1], "\n", file=f) 
        cat("  target", el[i,2], "\n", file=f) 
        for(n in names(eAttr)) {
            cat("   ", gsub("[\\._]", "", n), eAttr[[n]][i], "\n", file=f)
        }
        cat(" ]\n", file=f) 
    }
    
    cat("]\n", file=f)
    cat("Title \"", title, '"', file=f, sep="")
    close(f)
}

#' Export file to GML
#' @title save GML 
#' @description Export an igraph object to gexf format graph file which can be in opened in Gephi.
#' @param g: igraph object
#' @param filepath: output graph name
#' @return Outputs a Gephi format gexf graph
#' @export

saveAsGEXF = function(g, filepath="output.gexf")
{
    # gexf nodes require two column data frame (id, label)
    if(is.null(igraph::V(g)$label))
        igraph::V(g)$label <- as.character(igraph::V(g)$name)
    
    # similarily if edges does not have weight, add default 1 weight
    if(is.null(igraph::E(g)$weight))
        igraph::E(g)$weight <- rep.int(1, igraph::ecount(g))
    
    nodes <- data.frame(cbind(V(g), V(g)$label))
    edges <- t(Vectorize(igraph::get.edge, vectorize.args='id')(g, 1:igraph::ecount(g)))
    
    # combine all node attributes into a matrix (and take care of & for xml)
    vAttrNames <- setdiff(igraph::list.vertex.attributes(g), "label") 
    nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&#038;",igraph::get.vertex.attribute(g, attr))))
    
    # combine all edge attributes into a matrix (and take care of & for xml)
    eAttrNames <- setdiff(igraph::list.edge.attributes(g), "weight") 
    edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&#038;",igraph::get.edge.attribute(g, attr))))
    
    # combine all graph attributes into a meta-data
    graphAtt <- sapply(igraph::list.graph.attributes(g), function(attr) sub("&", "&#038;",igraph::get.graph.attribute(g, attr)))
    
    # generate the gexf object
    output <- write.gexf(nodes, edges, 
                         edgesWeight=igraph::E(g)$weight,
                         edgesAtt = edgesAtt,
                         nodesAtt = nodesAtt,
                         meta=c(list(description="igraph -> gexf converted file", keywords="igraph, gexf, R, rgexf"), graphAtt))
    
    print(output, filepath, replace=T)
}