


plot.mean = function (x, ...) 
{
	t.mean = mean(x)
	t.sd = sd(x)
	plot(x, ...)
	abline(h=t.mean, col="blue")
	abline(h=t.mean+t.sd, col="gray80")
	abline(h=t.mean-t.sd, col="gray80")
	abline(h=t.mean+2*t.sd, col="gray70")
	abline(h=t.mean-2*t.sd, col="gray70")
	abline(h=t.mean+3*t.sd, col="gray60")
	abline(h=t.mean-3*t.sd, col="gray60")
}


print.mean = function (x, ...) {
	t.mean = mean(x)
	t.sd = sd(x)
	print(paste("Media = ",t.mean))
	print(paste("Desvio Padrao = ",t.sd))
	print(paste("68%  -> ",   t.mean-t.sd,   ' <-> ',t.mean+t.sd))
	print(paste("95%  -> ",   t.mean-2*t.sd, ' <-> ',t.mean+2*t.sd))
	print(paste("99.7%  -> ", t.mean-3*t.sd, ' <-> ',t.mean+3*t.sd))
}



hist.felipe = function (x, ...) {
	sw.result = shapiro.test(x)
	if (sw.result$p.value < 0.1) { 
	sw.result2 = paste("Distr. Nao Normal (Shapiro-Walk p-value:", round(sw.result$p.value, digits = 6), ")" ) 
	} else {
	sw.result2 = paste("Distr. Normal (Shapiro-Walk p-value:",round(sw.result$p.value, digits = 6),")") }
	hist(x, prob=TRUE, sub=sw.result2, ...)
	lines(density(x)) 
}


plot.rootunit = function (x, ...) {
	ru.result1 = adf.test(x)
	if (ru.result1$p.value < 0.05) { 
	ru.result1.txt = paste("no unit-root (ADF p-value:", round(ru.result1$p.value, digits = 6), ")" ) 
	} else {
	ru.result1.txt= paste("contains a unit-root (ADF p-value:",round(ru.result1$p.value, digits = 6),")") }

	plot(x, sub=ru.result1.txt, ...)
	abline(h=0, col="gray70")
}

plot.ets = function (x, ...) 
{
    if (x$components[3] == "N" & x$components[2] == "N") {
        plot(cbind(observed = x$x, Base = x$states[, 1]),main = "Decomposi??o", ...)
    }
    else if (x$components[3] == "N") {
        plot(cbind(Observado = x$x, Base = x$states[, 1],  Crescimento= x$states[, 
            "b"]), main = "Decomposi??o", ...)
    }
    else if (x$components[2] == "N") {
        plot(cbind(Observado = x$x, Base = x$states[, 1], Sazonalidade= x$states[, 
            "s1"]), main = "Decomposi??o", ...)
    }
    else {
        plot(cbind(Observado = x$x, Base= x$states[, 1], Crescimento= x$states[, 
            "b"], Sazonalidade= x$states[, "s1"]), main = "Decomposi??o", ...)
    }
}

plot.div = function(x,y) {
	x11()
	# mar = margem da figura 
	par(mfrow=c(x,y), mar=c(2.5,2.5,1.5,1.5), omi=c(0, 0, 0, 0) );
}

require(datasets)
require(grDevices) 
require(graphics)

## --- "math annotation" in plots :

######
# create tables of mathematical annotation functionality
######
make.table <- function(nr, nc) {
  savepar <- par(mar=rep(0, 4), pty="s")
  plot(c(0, nc*2 + 1), c(0, -(nr + 1)),
	  type="n", xlab="", ylab="", axes=FALSE)
  savepar
}

get.r <- function(i, nr) {
  i %% nr + 1
}

get.c <- function(i, nr) {
  i %/% nr + 1
}

draw.title.cell <- function(title, i, nr) {
  r <- get.r(i, nr)
  c <- get.c(i, nr)
  text(2*c - .5, -r, title)
  rect((2*(c - 1) + .5), -(r - .5), (2*c + .5), -(r + .5))
}

draw.plotmath.cell <- function(expr1, expr2, i, nr) {
  r <- get.r(i, nr)
  c <- get.c(i, nr)
  text((2*(c - 1) + 1), -r, expr1, col="black")
  text((2*c), -r, expr2, adj=c(.5,.5))
  rect((2*(c - 1) + .5), -(r - .5), (2*c + .5), -(r + .5), border="grey")
}

draw.plotmath.cell2 <- function(expr, i, nr, string = NULL) {
  r <- get.r(i, nr)
  c <- get.c(i, nr)
  if (is.null(string)) {
	 string <- deparse(expr)
	 string <- substr(string, 12, nchar(string) - 1)
  }
  text((2*(c - 1) + 1), -r, string, col="grey")
  text((2*c), -r, expr, adj=c(.5,.5))
  rect((2*(c - 1) + .5), -(r - .5), (2*c + .5), -(r + .5), border="grey")
}


# do livro

lag.plot1=function(series,max.lag=1,corr=TRUE,smooth=FALSE){ 
   name1=paste(deparse(substitute(series)),"(t-",sep="")
   name2=paste(deparse(substitute(series)),"(t)",sep="")
   data1=as.ts(series)
   max.lag=as.integer(max.lag)
   prow=ceiling(sqrt(max.lag))
   pcol=ceiling(max.lag/prow)
   a=stats::acf(series,max.lag,plot=FALSE)$acf[-1]
   old.par <- par(no.readonly = TRUE)
   par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
  for(h in 1:max.lag){                       
   plot(lag(series,-h), data1, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
    lines(lowess(ts.intersect(lag(series,-h),series)[,1],
                 ts.intersect(lag(series,-h),series)[,2]), col="red")
    if (corr==TRUE)
    legend("topright", legend=round(a[h], digits=2), text.col ="blue", bg="white", x.intersp=0)
    on.exit(par(old.par))
   }
}

lag.plot2=function(series1,series2,max.lag=0,corr=TRUE,smooth=FALSE){ 
   name1=paste(deparse(substitute(series1)),"(t-",sep="")
   name2=paste(deparse(substitute(series2)),"(t)",sep="")
   data1=as.ts(series1)
   data2=as.ts(series2)
   max.lag=as.integer(max.lag)
   m1=max.lag+1
   prow=ceiling(sqrt(m1))
   pcol=ceiling(m1/prow)
   a=ccf(series1,series2,max.lag,plot=FALSE)$acf
   old.par <- par(no.readonly = TRUE)
   par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
   for(h in 0:max.lag){                   
   plot(lag(series1,-h), series2, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
    lines(lowess(ts.intersect(lag(series1,-h),series2)[,1],
                 ts.intersect(lag(series1,-h),series2)[,2]), col="red")
    if (corr==TRUE)
    legend("topright", legend=round(a[m1-h], digits=2), text.col ="blue", bg="white", x.intersp=0)             
   on.exit(par(old.par))
   }
}

acf2=function(series,max.lag=NULL){
  num=length(series)
  if (is.null(max.lag)) max.lag=ceiling(10+sqrt(num))
  if (max.lag > (num-1)) stop("Number of lags exceeds number of observations")
  ACF=acf(series, max.lag, plot=FALSE)$acf[-1]
  PACF=pacf(series, max.lag, plot=FALSE)$acf
  LAG=1:max.lag/frequency(series)
  minA=min(ACF)
  minP=min(PACF)
  U=2/sqrt(num)
  L=-U
  minu=min(minA,minP,L)-.01
  old.par <- par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar = c(3,3,2,0.8),
    oma = c(1,1.2,1,1), mgp = c(1.5,0.6,0))
  plot(LAG, ACF, type="h",ylim=c(minu,1), 
    main=paste("Series: ",deparse(substitute(series))))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  plot(LAG, PACF, type="h",ylim=c(minu,1))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  on.exit(par(old.par))  
  ACF<-round(ACF,2); PACF<-round(PACF,2)    
  return(cbind(ACF, PACF)) 
  }
  
  
  
  
HTMLplot <- function (Caption = "", file = get(".HTML.file"), append = TRUE, 
    GraphDirectory = ".", GraphFileName = "", GraphSaveAs = "png", 
    GraphBorder = 1, Align = "center", Width = 600, Height = 600, 
    WidthHTML = NULL, HeightHTML = NULL, GraphPointSize = 12, 
    GraphBackGround = "white", GraphRes = 72, plotFunction = NULL, 
    ...) 
{
    if (exists("HTMLenv", where = ".GlobalEnv", mode = "environment")) {
        GraphDirectory = get(".HTML.outdir", envir = get("HTMLenv", 
            envir = .GlobalEnv))
    }
    cat("\n", file = file, append = append, ...)
    if (GraphFileName == "") {
        nowd <- date()
        GraphFileName <- paste("GRAPH_", substring(nowd, 5, 7), 
            substring(nowd, 9, 10), "_", substring(nowd, 12, 
                13), substring(nowd, 15, 16), substring(nowd, 
                18, 19), sep = "")
    }
    GraphFileName <- paste(GraphFileName, ".", GraphSaveAs, sep = "")
    AbsGraphFileName <- file.path(GraphDirectory, GraphFileName)
    if (GraphSaveAs == "png") {
        if (is.null(plotFunction)) 
            dev.print(device = png, file = AbsGraphFileName, 
                width = Width, height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
        else {
            if (exists("X11", env = .GlobalEnv) && Sys.info()["sysname"] != 
                "Windows" && Sys.info()["sysname"] != "Darwin") 
                bitmap(file = AbsGraphFileName, bg = GraphBackGround, 
                  res = GraphRes)
            else png(filename = AbsGraphFileName, width = Width, 
                height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
            plotFunction()
            dev.off()
        }
    }
    else if (GraphSaveAs %in% c("jpg", "jpeg")) {
        if (is.null(plotFunction)) 
            dev.print(device = jpeg, file = AbsGraphFileName, 
                width = Width, height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
        else {
            if (exists("X11", env = .GlobalEnv) && Sys.info()["sysname"] != 
                "Windows" && Sys.info()["sysname"] != "Darwin") 
                bitmap(filename = AbsGraphFileName, bg = GraphBackGround, 
                  res = GraphRes, type = "jpeg")
            else jpeg(filename = AbsGraphFileName, width = Width, 
                height = Height, pointsize = GraphPointSize, 
                bg = GraphBackGround)
            plotFunction()
            dev.off()
        }
    }
    else if (GraphSaveAs == "gif") {
        stop("Gif support was removed from base R because of patent restrictions. Use either jpg or png")
    }
    else stop("GraphSaveAs must be either jpg, png or gif")
    cat(paste("<img src='", GraphFileName, 
        "' border=", GraphBorder, if (!is.null(Width)) 
            paste(" width=", Width, sep = "")
        else "", if (!is.null(HeightHTML)) 
            paste(" height=", HeightHTML, sep = ""), if (!is.null(WidthHTML)) 
            paste(" width="), ">\n", sep = "", collapse = ""), 
        file = file, append = TRUE, sep = "")
    if (substitute(file) == ".HTML.file") 
        try(assign(".HTML.graph", TRUE, env = get("HTMLenv", 
            envir = .GlobalEnv)))
    invisible(return(TRUE))
}




