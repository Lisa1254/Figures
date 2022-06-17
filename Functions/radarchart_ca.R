#Modification of radarchart in R package fmsb to allow for custom axis and highlight of background grid to better demonstrate biological relevance
#https://rdrr.io/cran/fmsb/src/R/fmsb.R#sym-radarchart

#Added argument num_ax for values along grid polygon axes
#Since full axis is now provided, input df now only includes data rows. No more need for max/min (circumference/centre) rows of df
#If caxislabels=NULL, axes will be labelled with the same as the provided num_ax
#Removed some of the axistype options - now choice is either 0 for none, or 1 for labelled
#With additional axistype options removed, also removed arguments for paxislabel and palcex
#seg argument removed from function definition. Number of segments is now defined in function as the length of the provided num_ax - 1
#removed centrezero option
#Currently no support for imputing or handling NA values, as in original function; needs to be done to data before using function
#Outer 2 segments of grid will be highlighted in light gray, but intend to add customization for this option

radarchart_ca <- function(df, num_ax=c(18000,1000,100,10,1), caxislabels=NULL,
                          axistype=1, pty=16, pcol=1:8, plty=1:6, plwd=1,
                          pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                          cglcol="navy", axislabcol="blue", title="",
                          na.itp=TRUE, vlabels=NULL, vlcex=NULL,
                          calcex=NULL, ...){
  #Checks for input of data
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }

  
  #Blank plot background
  # here the ... is reproduced to pass any additional plot arguments from the radarchart
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1, ...) # define x-y coordinates without any plot
  
  #Grid Shape
  #Set up transformation factors for polar coordinates
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  
  #Define number of segments according to provided axis scale
  seg <- length(num_ax) - 1
  
  #Define axis labels if not provided
  if(is.null(caxislabels)) {
    caxislabels <- as.character(num_ax)
  }
  
  #Concentric Grid (with labels if including)
  for (i in c(seg, seg-1)) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+1)/(seg+1), yy*(i+1)/(seg+1), lty=cglty, lwd=cglwd, border=cglcol, col="gray94")
    if (axistype==1) {
      if (is.null(calcex)) text(-0.05, (i+1)/(seg+1), caxislabels[i+1], col=axislabcol) else
        text(-0.05, (i+1)/(seg+1), caxislabels[i+1], col=axislabcol, cex=calcex)
    }
  }
  seg_inner <- 0:(seg-2)
  seg_inner <- seg_inner[order(seg_inner, decreasing = TRUE)]
  for (i in seg_inner) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+1)/(seg+1), yy*(i+1)/(seg+1), lty=cglty, lwd=cglwd, border=cglcol, col="white")
    if (axistype==1) {
      if (is.null(calcex)) text(-0.05, (i+1)/(seg+1), caxislabels[i+1], col=axislabcol) else
        text(-0.05, (i+1)/(seg+1), caxislabels[i+1], col=axislabcol, cex=calcex)
    }
  }
  
  #Grid Spokes
  arrows(xx/(seg+1), yy/(seg+1), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)

  
  #Variable name, the labels at the vertices
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) VLABELS <- vlabels
  if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
    text(xx*1.2, yy*1.2, VLABELS, cex=vlcex)
  
  #Set up for multiple datasets on same grid
  SX <- length(df[[1]])
  #SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  
  #Iterates for each row in dataframe, not including the max/min rows
  for (i in 1:SX) { #series is number of datasets being plotted on same grid
    xxs <- xx
    yys <- yy
    for(j in 1:n){ #n is number of verticies
      section <- seg
      while (df[i, j] > num_ax[section]) {
        section <- section - 1
        #consider adding fail-safe for if a value is out of range, i.e. while also >=1
      }
      pf <- (num_ax[section] - df[i,j])/(num_ax[section]-num_ax[section+1])
      xxs[j] <- xx[j]*((pf + section)/(seg+1))
      yys[j] <- yy[j]*((pf + section)/(seg+1))
      
    }
    
    #pdensity makes lines crossing the data polygon
    if (is.null(pdensities)) {
      polygon(xxs, yys, lty=pltys[i], lwd=plwds[i], border=pcols[i], col=pfcols[i])
    } else {
      polygon(xxs, yys, lty=pltys[i], lwd=plwds[i], border=pcols[i], 
              density=pdensities[i], angle=pangles[i], col=pfcols[i])
    }
    
    #draw Point shape for data value
    points(xxs, yys, pch=ptys[i], col=pcols[i])
    
    
  }
  
}

