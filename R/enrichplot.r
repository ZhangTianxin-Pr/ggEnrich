#'enrichplot function is used to draw the enrichment plot.
#'@export
#'@description
#'enrichplot function is used to draw the enrichment plot.
#'@param file pathway file, the format is xls
#'@param type the type of plot, dot,line,barcode
#'@param mycolor the type of plot, dot,line,barcode
#'@param plot_legend location of the legend,left,right
#'@return a plot, the format is pdf
#'@importFrom grDevices dev.off
#'@importFrom utils read.delim
#'@examples
#'path <- "D:\\code\\r\\FigureYa13GSEA_Java_update"
#'file <- mergefile(path)
#'enrichplot <-enrichplot(file,type="dot")
#'enrichplot <-enrichplot(file,type="line",plot_legend="left")
#'enrichplot <-enrichplot(file,type="barcode")



enrichplot <-function(file,type,mycolor=c("red","navy","darkgreen","blueviolet","chocolate4","black"),plot_legend="right"){
  if(type == "dot"){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      ppoint <- ggplot2::ggplot(file, ggplot2::aes(x=file$RANK.IN.GENE.LIST,y=file$RUNNING.ES,fill=pathway,group=pathway))+
      ggplot2::geom_point(shape=21) +
      ggplot2::scale_fill_manual(values = mycolor) +
      ggplot2::labs(x = "", y = "Enrichment Score", title = "") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),limits =c(min(file$RUNNING.ES-0.02), max(file$RUNNING.ES+0.02))) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank(),axis.text.x = ggplot2::element_blank())+
      ggplot2::geom_hline(yintercept = 0)
      if(plot_legend == "left"){
        ppoint <- ppoint +ggplot2::theme(legend.position=c(0,0),legend.justification = c(0,0)) +
        ggplot2::guides(fill=ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(legend.background = ggplot2::element_blank()) +
        ggplot2::theme(legend.key = ggplot2::element_blank())
        }
      prank <- ggplot2::ggplot(file,ggplot2::aes(file$RANK.IN.GENE.LIST,pathway,colour=pathway))+
      ggplot2::geom_tile()+
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = mycolor) +
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
    if (requireNamespace("grid", quietly = TRUE)) {
      grid::grid.newpage()
    }
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(gridExtra::arrangeGrob(gA,gB,nrow=2,heights=c(.8,.3)))
    }
    pdf<-pdf('prettyGSEApoint.pdf',width=8,height=4)
    dev.off()
    return(pdf)
    }
  }
  if(type == "line"){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      pline <- ggplot2::ggplot(file,ggplot2::aes(x=file$RANK.IN.GENE.LIST,y=file$RUNNING.ES,colour=pathway,group=pathway))+
        ggplot2::geom_line(size=1)+
        ggplot2::scale_color_manual(values = mycolor) +
        ggplot2::labs(x = "", y = "Enrichment Score", title = "") +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0),limits =c(min(file$RUNNING.ES-0.02), max(file$RUNNING.ES+0.02))) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid =ggplot2::element_blank()) +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),axis.ticks.x = ggplot2::element_blank(),axis.text.x = ggplot2::element_blank()) +
        ggplot2::geom_hline(yintercept = 0)
        if(plot_legend == "left"){
          pline<-pline+ggplot2::theme(legend.position=c(0,0),legend.justification = c(0,0)) +
          ggplot2::guides(colour=ggplot2::guide_legend(title = NULL)) +
          ggplot2::theme(legend.background = ggplot2::element_blank()) +
          ggplot2::theme(legend.key = ggplot2::element_blank())
        }

      prank <- ggplot2::ggplot(file,ggplot2::aes(file$RANK.IN.GENE.LIST,pathway,colour=pathway))+
        ggplot2::geom_tile()+
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = mycolor) +
        ggplot2::labs(x = "treatment<--------------control", y = "", title = "") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::theme(panel.grid =ggplot2::element_blank()) +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))+
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),axis.ticks.y = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank())+
        ggplot2::guides(color=FALSE)

      gA=ggplot2::ggplot_gtable(ggplot2::ggplot_build(pline))
      gB=ggplot2::ggplot_gtable(ggplot2::ggplot_build(prank))
      maxWidth = grid::unit.pmax(gA$widths, gB$widths)
      gA$widths <- as.list(maxWidth)
      gB$widths <- as.list(maxWidth)
      if (requireNamespace("grid", quietly = TRUE)) {
        grid::grid.newpage()
      }
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(gridExtra::arrangeGrob(gA,gB,nrow=2,heights=c(.8,.3)))
      }
      pdf<-pdf('prettyGSEAline.pdf',width=8,height=4)
      dev.off()
      return(pdf)
    }
  }

  if(type == "barcode"){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      pbar <- lapply(unique(file$pathway), function(ii) {
        dd <- file[file$pathway == ii,]
       ggplot2::ggplot(dd,ggplot2::aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES))+
        ggplot2::geom_bar(stat="identity",colour = "black")+
        ggplot2::labs(x = ii, y = "", title = ii) +
        ggplot2::theme_void() +
        ggplot2::xlim(0, max(file$RANK.IN.GENE.LIST))
      })
      t <- data.frame(a=as.numeric(1:1000),b=as.numeric(1:1000))
      pheat <- ggplot2::ggplot(t,ggplot2::aes(t$a,1,fill=t$b)) +
        ggplot2::geom_tile()+
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Decreased--------------Increased", y = "", title = "") +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 500)+
        ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::theme(panel.grid =ggplot2::element_blank()) +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::guides(fill=FALSE)
      num<-length(pbar)+1
      pbar[[num]] <- pheat
      if (requireNamespace("cowplot", quietly = TRUE)) {
        rel_heights<-c(rep(2,num-1),1)
        cowplot::plot_grid(plotlist=pbar, ncol=1,rel_heights = rel_heights)
      }
      ggplot2::ggsave("prettyGSEAbar.pdf")
    }
  }

}

