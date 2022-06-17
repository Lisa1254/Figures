#Modification of radarchart in R package fmsb to allow for custom axis and highlight of background grid to better demonstrate biological relevance
#https://rdrr.io/cran/fmsb/src/R/fmsb.R#sym-radarchart

#Added argument num_ax

radarchart_ca <- function(df, num_ax=c(18000,1000,100,10,1), caxislabels=NULL,
                          axistype=0, seg=4, pty=16, pcol=1:8, plty=1:6, plwd=1,
                          pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                          cglcol="navy", axislabcol="blue", title="", maxmin=TRUE,
                          na.itp=TRUE, centerzero=FALSE, vlabels=NULL, vlcex=NULL,
                          calcex=NULL, paxislabels=NULL, palcex=NULL, ...){
  #Checks for input of data
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }
  if (maxmin==FALSE) { # when the dataframe does not include max and min as the top 2 rows.
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  
  #Blank plot background
  # here the ... is reproduced to pass any additional plot arguments from the radarchart
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1) # define x-y coordinates without any plot
  
  #Grid Shape
  #Set up transformation factors for polar coordinates
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  
  #Variable for centre zero
  CGap <- ifelse(centerzero, 0, 1)
  
  #Concentric Grid (with some label information)
  for (i in 0:seg) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+CGap)/(seg+CGap), yy*(i+CGap)/(seg+CGap), lty=cglty, lwd=cglwd, border=cglcol)
    if (axistype==1|axistype==3) CAXISLABELS <- paste(i/seg*100,"(%)")
    if (axistype==4|axistype==5) CAXISLABELS <- sprintf("%3.2f",i/seg)
    if (!is.null(caxislabels)&(i<length(caxislabels))) CAXISLABELS <- caxislabels[i+1]
    if (axistype==1|axistype==3|axistype==4|axistype==5) {
      if (is.null(calcex)) text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol) else
        text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol, cex=calcex)
    }
  }
  
  #Grid Spokes
  if (centerzero) {
    arrows(0, 0, xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }else {
    arrows(xx/(seg+CGap), yy/(seg+CGap), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  
  #Peripheral of grid labels, if included according to the axis type
  PAXISLABELS <- df[1,1:n]
  if (!is.null(paxislabels)) PAXISLABELS <- paxislabels
  if (axistype==2|axistype==3|axistype==5) {
    if (is.null(palcex)) text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol) else
      text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol, cex=palcex)
  }
  
  #Variable name, the labels at the vertices
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) VLABELS <- vlabels
  if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
    text(xx*1.2, yy*1.2, VLABELS, cex=vlcex)
  
  #Set up for multiple datasets on same grid
  series <- length(df[[1]])
  SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  
  #Iterates for each row in dataframe, not including the max/min rows
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    for(j in 1:n){ #n is number of verticies
      section <- seg
      while (df[3, j] > num_ax[section]) {
        section <- section - 1
        #consider adding fail-safe for if a value is out of range, i.e. while also >=1
      }
      pf <- (num_ax[section] - df[3,j])/(num_ax[section]-num_ax[section+1])
      xxs[j] <- xx[j]*((pf + section)/(seg+1))
      yys[j] <- yy[j]*((pf + section)/(seg+1))
      
      
    }
    #pdensity makes lines crossing the data polygon
    if (is.null(pdensities)) {
      polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], col=pfcols[i-2])
    } else {
      polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], 
              density=pdensities[i-2], angle=pangles[i-2], col=pfcols[i-2])
    }
    
    #draw Point shape for data value
    points(xxs, yys, pch=ptys[i-2], col=pcols[i-2])
    
    
  }
  
}

