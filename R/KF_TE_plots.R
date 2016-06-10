plotTE_FCs <- function(object, log2FC.cutoff = 1, cex=1, xlim, ylim, ..., cex.lab, cex.axis, cex.main, cex.sub, pch)
{
  list.of.packages <- "scales"
  if(list.of.packages %in% installed.packages()[,"Package"]){
    library("scales")
  }else{
    message('please install "scales" before run this function')
  }

  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")

  if (log2FC.cutoff <= 0 ) stop("log2FC.cutoff must be larger than 0")
  if (!missing(xlim)){
    if (length(xlim) != 2) stop("invalid xlim value, should have two elemtnts")
  }
  if (!missing(ylim)){
    if (length(ylim) != 2) stop("invalid ylim value, should have two elemtnts")
  }
  resultsTable <- object$resultsTable
  resultsTable <- resultsTable[complete.cases(resultsTable),]
  resultsTable <- resultsTable[order(resultsTable$pvalue_final, decreasing = TRUE),]
  resultsTable$colorscale <- rescale(-log10(resultsTable$pvalue_final))

  if (missing(xlim)){
    xmin <- min(resultsTable$mRNA_log2FC,na.rm=T) - 0.1
    xmax <- max(resultsTable$mRNA_log2FC,na.rm=T) + 0.1
    xlim <- c(xmin,xmax)
  }
  if (missing(ylim)){
    ymin <- min(resultsTable$log2FC_TE_final,na.rm=T) - 0.1
    ymax <- max(resultsTable$log2FC_TE_final,na.rm=T) + 0.1
    ylim <- c(ymin,ymax)
  }
  if (missing(cex.lab)) cex.lab <- 1.2
  if (missing(cex.axis)) cex.axis <- 1
  if (missing(cex.main)) cex.main <- 1.2
  if (missing(cex.sub)) cex.sub <- 1
  if (missing(pch)) pch <- 20

  #mRNA stable, RPF stable.
  variable <- which(abs(resultsTable$mRNA_log2FC - resultsTable$log2FC_TE_final) <= log2FC.cutoff |
                      (abs(resultsTable$mRNA_log2FC) < log2FC.cutoff & abs(resultsTable$log2FC_TE_final) < log2FC.cutoff))
  plot(resultsTable$mRNA_log2FC[variable],resultsTable$log2FC_TE_final[variable],pch=pch,col="gray90",
       xlim=xlim,xlab = "mRNA log2FC", ylab="TE log2FC",frame.plot=TRUE,cex=cex,
       ylim=ylim,
       cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, cex.sub=cex.sub, ...)

  #mRNA change, RPF stable.
  variable <- which( abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
                       abs(resultsTable$log2FC_TE_final) < log2FC.cutoff &
                       abs(resultsTable$mRNA_log2FC - resultsTable$log2FC_TE_final) >= log2FC.cutoff)
  colornames <- seq_gradient_pal("#BFBBFF","#0D5FFF")(resultsTable$colorscale[variable])
  points(resultsTable$mRNA_log2FC[variable],resultsTable$log2FC_TE_final[variable],
         pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- "transcription only"
  leg.col <- "#0D5FFF"

  #mRNA stable, RPF change.
  variable <- which( abs(resultsTable$mRNA_log2FC) < log2FC.cutoff &
                       abs(resultsTable$log2FC_TE_final) >= log2FC.cutoff &
                       abs(resultsTable$mRNA_log2FC - resultsTable$log2FC_TE_final) >= log2FC.cutoff )
  colornames <- seq_gradient_pal("#FFB69C","#FF2E06")(resultsTable$colorscale[variable])
  points(resultsTable$mRNA_log2FC[variable],resultsTable$log2FC_TE_final[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, "translation only")
  leg.col <- c(leg.col, "#FF2E06")

  #mRNA change, RPF change, homodirectional.
  variable <- which( sign(resultsTable$mRNA_log2FC) * sign(resultsTable$log2FC_TE_final) > 0 &
                       abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
                       abs(resultsTable$log2FC_TE_final) >= log2FC.cutoff &
                       abs(resultsTable$mRNA_log2FC - resultsTable$log2FC_TE_final) >= log2FC.cutoff )
  colornames <- seq_gradient_pal("#D1F7AD","#75E805")(resultsTable$colorscale[variable])
  points(resultsTable$mRNA_log2FC[variable],resultsTable$log2FC_TE_final[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, "homodirectional")
  leg.col <- c(leg.col, "#75E805")

  #mRNA change, RPF change, opposite change.
  variable <- which( sign(resultsTable$mRNA_log2FC) * sign(resultsTable$log2FC_TE_final) < 0 &
                       abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
                       abs(resultsTable$log2FC_TE_final) >= log2FC.cutoff &
                       abs(resultsTable$mRNA_log2FC - resultsTable$log2FC_TE_final) >= log2FC.cutoff )
  colornames <- seq_gradient_pal("#FFF1AF","#FFDE13")(resultsTable$colorscale[variable])
  points(resultsTable$mRNA_log2FC[variable],resultsTable$log2FC_TE_final[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, "opposite change")
  leg.col <- c(leg.col, "#FFDE13")

  abline(h= log2FC.cutoff, lty=2, col="gray")
  abline(v= log2FC.cutoff, lty=2, col="gray")
  abline(h= -log2FC.cutoff, lty=2, col="gray")
  abline(v= -log2FC.cutoff, lty=2, col="gray")

  abline(h= 0, lty=2, col="gray")
  abline(v= 0, lty=2, col="gray")


  legend("bottomright", legend=leg,pch=pch, col=leg.col, bty="n",cex=cex)
}

#=====================================================================================================================
#=====================================================================================================================

# KF Work 9.6.16
if(!exists("results_arabidopsis_xtail_analysis_17C_28C_raw_min50_100000")){
  load("~/Documents/Betty_Plant_Data/Arabidopsis/arabidopsis_xtail_analysis_17C_28C_attempt1.RData")
  load("~/Documents/Betty_Plant_Data/Arabidopsis/arabidopsis_xtail_analysis_17C_4C_attempt1.RData")
  plotTE_FCs(arabidopsis_xtail_analysis_17C_28C,main="17 C vs. 28 C")
}
