#'enrichplot function is used to draw the enrichment plot.
#'@export
#'@description
#'enrichplot function is used to draw the enrichment plot.
#'@param file pathway file, the format is xls
#'@param type the type of plot, dot,line,barcode
#'@param mycolor the type of plot, dot,line,barcode,clusterprofiler
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
      ppoint <- ggplot2::ggplot(file, ggplot2::aes(x=file$RANK.IN.GENE.LIST,
                                                   y=file$RUNNING.ES,
                                                   fill=pathway,group=pathway))+
      ggplot2::geom_point(shape=21) +
      ggplot2::scale_fill_manual(values = mycolor) +
      ggplot2::labs(x = "", y = "Enrichment Score", title = "") +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  limits =c(min(file$RUNNING.ES-0.02),
                                            max(file$RUNNING.ES+0.02))) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank())+
      ggplot2::geom_hline(yintercept = 0)
      if(plot_legend == "left"){
        ppoint <- ppoint +ggplot2::theme(legend.position=c(0,0),
                                         legend.justification = c(0,0)) +
        ggplot2::guides(fill=ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(legend.background = ggplot2::element_blank()) +
        ggplot2::theme(legend.key = ggplot2::element_blank())
        }
      prank <- ggplot2::ggplot(file,ggplot2::aes(file$RANK.IN.GENE.LIST,
                                                 pathway,
                                                 colour=pathway))+
      ggplot2::geom_tile()+
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = mycolor) +
      ggplot2::labs(x = "treatment<--------------control", y = "", title = "") +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::theme(panel.grid =ggplot2::element_blank()) +
      ggplot2::theme(panel.border = ggplot2::element_blank()) +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))+
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank())+
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
      pline <- ggplot2::ggplot(file,ggplot2::aes(x=file$RANK.IN.GENE.LIST,
                                                 y=file$RUNNING.ES,
                                                 colour=pathway,
                                                 group=pathway))+
        ggplot2::geom_line(size=1)+
        ggplot2::scale_color_manual(values = mycolor) +
        ggplot2::labs(x = "", y = "Enrichment Score", title = "") +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0),
                                    limits =c(min(file$RUNNING.ES-0.02),
                                              max(file$RUNNING.ES+0.02))) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid =ggplot2::element_blank()) +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black")) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank()) +
        ggplot2::geom_hline(yintercept = 0)
        if(plot_legend == "left"){
          pline<-pline+ggplot2::theme(legend.position=c(0,0),
                                      legend.justification = c(0,0)) +
          ggplot2::guides(colour=ggplot2::guide_legend(title = NULL)) +
          ggplot2::theme(legend.background = ggplot2::element_blank()) +
          ggplot2::theme(legend.key = ggplot2::element_blank())
        }

      prank <- ggplot2::ggplot(file,ggplot2::aes(file$RANK.IN.GENE.LIST,
                                                 pathway,
                                                 colour=pathway))+
        ggplot2::geom_tile()+
        ggplot2::theme_bw() +
        ggplot2::scale_color_manual(values = mycolor) +
        ggplot2::labs(x = "treatment<--------------control", y = "", title = "") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::theme(panel.grid =ggplot2::element_blank()) +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"))+
        ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank())+
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
        ggplot2::labs(x = "Decreased--------------Increased",
                      y = "", title = "") +
        ggplot2::scale_fill_gradient2(low = "blue",
                                      mid = "white",
                                      high = "red",midpoint = 500)+
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
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
  if(type == "clusterprofiler"){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p.res <- ggplot2::ggplot(file, ggplot2::aes_(x = ~x)) + xlab(NULL) +
        ggplot2::geom_line(ggplot2::aes_(y = ~runningScore,
                                         color= ~Description), size=1) +
        ggplot2::scale_color_manual(values = mycol) +
        ggplot2::geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
        ggplot2::ylab("Enrichment\n Score") +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::theme(legend.position = "top",
                       legend.title = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(fill = "transparent")) +
        ggplot2::theme(axis.text.y=ggplot2::element_text(size = 12, face = "bold"),
              axis.text.x=ggplot2::element_blank(),
              axis.ticks.x=ggplot2::element_blank(),
              axis.line.x=ggplot2::element_blank(),
              plot.margin=ggplot2::margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
      rel_heights <- c(1.5, .5, 1.5)
      i <- 0
      for (term in unique(file$Description)) {
        idx <- which(file$ymin != 0 & file$Description == term)
        file[idx, "ymin"] <- i
        file[idx, "ymax"] <- i + 1
        i <- i + 1
      }
      p2 <- ggplot2::ggplot(file, ggplot2::aes_(x = ~x)) +
        ggplot2::geom_linerange(ggplot2::aes_(ymin=~ymin,
                                              ymax=~ymax,
                                              color=~Description)) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggplot2::scale_color_manual(values = mycol)+
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::theme(legend.position = "none",
              plot.margin = ggplot2::margin(t=-.1, b=0,unit="cm"),
              axis.ticks = ggplot2::element_blank(),
              axis.text = ggplot2::element_blank(),
              axis.line.x = ggplot2::element_blank()) +
        ggplot2::scale_y_continuous(expand=c(0,0))
      df2 <- p.res$data
      df2$y <- p.res$data$geneList[df2$x]
      df2$gsym <- p.res$data$gsym[df2$x]
      selectgenes <- data.frame(gsym = selectedGeneID)
      selectgenes <- merge(selectgenes, df2, by = "gsym")
      selectgenes <- selectgenes[selectgenes$position == 1,]
      p.pos <- ggplot2::ggplot(selectgenes,
                               ggplot2::aes(x, y,
                                   fill = Description,
                                   color = Description,
                                   label = gsym)) +
        ggplot2::geom_segment(data=df2, ggplot2::aes_(x=~x, xend=~x, y=~y, yend=0),
                     color = "grey") +
        geom_bar(position = "dodge", stat = "identity") +
        scale_fill_manual(values = mycol, guide=FALSE) +
        scale_color_manual(values = mycol, guide=FALSE) +
        geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
        ylab("Ranked list\n metric") +
        xlab("Rank in ordered dataset") +
        theme_bw() +
        theme(axis.text.y=element_text(size = 12, face = "bold"),
              panel.grid = element_blank()) +
        geom_text_repel(data = selectgenes,
                        show.legend = FALSE,
                        direction = "x",
                        ylim = c(2, NA),
                        angle = 90,
                        size = 2.5, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines")) +
        theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
      plotlist <- list(p.res, p2, p.pos)
      n <- length(plotlist)
      plotlist[[n]] <- plotlist[[n]] +
        theme(axis.line.x = element_line(),
              axis.ticks.x = element_line(),
              axis.text.x = element_text(size = 12, face = "bold"))
      plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)
      pdf<-pdf('GSEA_multi_pathways.pdf',width=6,height=5)
      dev.off()
      return(pdf)
     }
    }
}

