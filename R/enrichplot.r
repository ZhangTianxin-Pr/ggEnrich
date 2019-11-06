#'enrichplot function is used to draw the enrichment plot.
#'@export
#'@description
#'enrichplot function is used to draw the enrichment plot.
#'@param file pathway file, the format is xls
#'@param type the type of plot, dot,line,barcode
#'@return a plot, the format is pdf
#'@examples
#'path="D:\\code\\r\\FigureYa13GSEA_Java_update"
#'file <- mergefile(path)
#'enrichplot <-enrichplot(file,type="dot")



enrichplot <-function(file,type){
  mycol <- c("red","navy","darkgreen","blueviolet","chocolate4")
  if(identical(type,"dot")){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
    ppoint <- ggplot2::ggplot(file, ggplot2::aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
      ggplot2::geom_point(shape=21) +
      ggplot2::scale_fill_manual(values = mycol) +
      ggplot2::labs(x = "", y = "Enrichment Score", title = "") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),limits =c(min(file$RUNNING.ES-0.02), max(file$RUNNING.ES+0.02))) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank(),axis.text.x = ggplot2::element_blank())+
      ggplot2::geom_hline(yintercept = 0)

    prank <- ggplot2::ggplot(file,ggplot2::aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+
      ggplot2::geom_tile()+
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = mycol) +
      ggplot2::labs(x = "treatment<--------------control", y = "", title = "") +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(panel.grid =ggplot2::element_blank()) +
      ggplot2::theme(panel.border = ggplot2::element_blank()) +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),axis.ticks.y = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank())+
      ggplot2::guides(color=FALSE)

    gA=ggplot2::ggplot_gtable(ggplot2::ggplot_build(ppoint))
    gB=ggplot2::ggplot_gtable(ggplot2::ggplot_build(prank))
    maxWidth = grid::unit.pmax(gA$widths, gB$widths)
    gA$widths <- as.list(maxWidth)
    gB$widths <- as.list(maxWidth)
    grid.newpage()
    grid.arrange(arrangeGrob(gA,gB,nrow=2,heights=c(.8,.3)))
    pdf<-pdf('prettyGSEApoint.pdf',width=8,height=4)
    dev.off()
    return(pdf)
    }else{install.packages("ggplot2", repos = "http://cran.r-project.org")}
  }


}
