#'mergefile function is used to merge all xls files in the working path.
#'@export
#'@description
#'mergefile function is used to merge all xls files in the working path.
#'@param path the parameter is the working path and should be a string
#'@param format the way to get the results of the enrichment analysis,JAVA_GSEA,R_ClusterProfiler
#'@param type multiple pathways or multiple groups of one pathway,pathways,groups
#'@return a data.frame.
#'@examples
#'path="D:\\code\\r\\FigureYa13GSEA_Java_update"
#'result <- mergefile(path)

preprocess <- function(path,format,type="pathways") {
  setwd(path)
  if(format == "JAVA_GSEA"){
    if(type == "pathways"){
      fnames <- Sys.glob("*.xls")
      fdataset <- lapply(fnames,read.delim)
      names(fdataset) <- fnames
      if (requireNamespace("plyr", quietly = TRUE)) {
        result <- plyr::ldply(fdataset, data.frame)
      }
      result$pathway <- unlist(strsplit(result$.id,split = ".xls"))
      return(result)
    }
    if(type == "groups"){
      fnames <- Sys.glob("*.xls")
      fdataset <- lapply(fnames,read.delim)
      names(fdataset) <- fnames
      if (requireNamespace("plyr", quietly = TRUE)) {
        result <- plyr::ldply(fdataset, data.frame)
      }
      result$group <- unlist(strsplit(result$.id,split = ".xls"))
      return(result)
    }
  }
  if(format == "R_ClusterProfiler"){
    if(type == "pathways"){
      x <- read.csv("gsea_output.csv",header = TRUE,fileEncoding = "utf-8")
    }
    if(type == "groups"){
      fnames <- Sys.glob("*.csv")
      fdataset <- lapply(fnames,read.delim)
      names(fdataset) <- fnames
      if (requireNamespace("plyr", quietly = TRUE)) {
        result <- plyr::ldply(fdataset, data.frame)
      }
      return(result)
    }
  }

}
