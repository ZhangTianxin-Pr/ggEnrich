#'mergefile function is used to merge all xls files in the working path.
#'@export
#'@description
#'mergefile function is used to merge all xls files in the working path.
#'@param path the parameter is the working path and should be a string
#'@return a data.frame.
#'@examples
#'path="D:\\code\\r\\FigureYa13GSEA_Java_update"
#'result <- mergefile(path)

mergefile <- function(path) {
  setwd(path)
  fnames <- Sys.glob("*.xls")
  fdataset <- lapply(fnames,read.delim)
  names(fdataset) <- fnames
  result <- plyr::ldply(fdataset, data.frame)
  result$pathway <- unlist(strsplit(result$.id,split = ".xls"))
  return(result)
}
