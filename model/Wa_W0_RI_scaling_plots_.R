#load data & libaries
library(stats); library(smatr)
setwd("/Users/dfalster/Ecology/_games/offspring size/data");

#returns vector from lo to hi with multiplication steps of incr
seq.log<-function(lo, hi, incr){temp<-NULL; while(lo < hi) {temp<-c(temp, lo); lo=lo*incr;}
					  if(max(temp)<hi) temp<-c(temp, hi); return (temp);} 
#log-log plot with legend and 1:1
	ll_plot=function(X,Y, XLAB, YLAB){
		plot(X, Y, log="xy", xlab=XLAB, ylab=YLAB)
		fit<-slope.test(log10(Y), log10(X), test.value = 1.0, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2) 
		r2<- cor(log10(Y), log10(X))^2; n = length(X)
  		leg<-c(paste("B = ",format(fit$b, digits =4)), paste("Int = ",format(Int, digits =4)),paste("n = ",n), paste("R2 = ",format(r2, digits =3)))
  		legend(x="topleft",legend=leg, lty=0, cex=1.0)}
#-----------------------------------------------------
#PLANTS
  Raw <- read.table("Plants.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
  Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$W0))

#MAMMALS
  Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  Dat <-subset(Raw, !is.na(Raw$W_newborn )&!is.na(Raw$W0))
  
  Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
  Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$Ea))
  Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI)&!is.na(Raw$Ea)&!is.na(Raw$W0))
  #data reduced to complete set for Wa, W0, E, RI


xlab="Adult body size(kg) [log scale]"; 
#Annual repro investment
  Y<- Dat$RI;  X<- Dat$Wa; ylab = "Annaul Reproductive Investment (kg) [log scale]";
#Lifetime repro investment
  Y<- Dat$RI*Dat$Ea;  X<- Dat$Wa; ylab = "Annaul Reproductive Investment (kg) [log scale]";
#Average adult lifespan
  Y<- Dat$Ea;  X<- Dat$Wa; ylab = "Reproductive Lifespan (yr) [log scale]";
#max repro lifespan 
  Y<- Dat$Ea_max;  X<- Dat$Wa; ylab = "Reproductive Lifespan (yr) [log scale]";
#W0 vs Wa
  Y<- Dat$W0;  X<- Dat$Wa; ylab = "Size at weaning (kg) [log scale]"; xlab = "Size adults (kg) [log scale]";
#W0/Wa vs Wa
  Y<- Dat$W0/Dat$Wa;  X<- Dat$Wa; ylab = "Ratio offspring:adult mass [log scale]";

#NEW BORN VS WEANING WEIGHT
  Dat <-subset(Raw, !is.na(Raw$W0)&!is.na(Raw$W_newborn))
  Y<-Dat$W_newborn; X<-Dat$W0; xlab="Size at Weaning(kg) [log scale]"; ylab = "Size at Birth(kg) [log scale]";
  
#ea_max vs Ea
  Dat<- subset(Raw, !is.na(Raw$Ea)&!is.na(Raw$Ea_max)&!is.na(Raw$Wa))
  X = Dat$Ea_max; Y = Dat$Ea; xlab = "Maximum lifespan (yrs)"; ylab ="Average lifespan (yrs)";

ll_plot(X,Y, xlab, ylab)
#--------------------------------------------------------------------------------
#publication plots
#2 row version 
postscript("fig3.eps", horizontal=FALSE, onefile=FALSE, height=6,width=10,pointsize=10)
par(mfrow=c(1,2))
#1. Annual reproductive output
	 Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
   Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
	 X= Dat$Wa; Y=Dat$RI;
	 Y_ax=seq.log(10^-6, 10^4, 10);	X_ax=seq.log(10^-6, 10^4, 10);
	 Y_lab =c(expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3), expression(10^4))	
	 X_lab =  Y_lab
   	
	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="darkgray", pch = 19, cex=1.2, xaxs="i", yaxs="i", family="helvetica")
		points(X, Y, type="p",pch = 1, cex=1.3,lwd =0.25)
		
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext("Annaul reproductive investment (kg)", side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
		mtext("Adult size (kg)", side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
	
		#fit slope
		fit<-slope.test(log10(Y), log10(X), test.value = 0.75, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=1.5, col ="darkgray") 

	#LOAD PLANT DATA
	Raw <- read.table("Plants.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
 	Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
      X= Dat$Wa; Y=Dat$RI;
	points(X, Y, type="p", pch=1, cex=1.3, lwd =0.25)
	fit<-slope.test(log10(Y), log10(X), test.value = 0.75, method = 'SMA')
	Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
	curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=1.5) 

	
#3. Lifetime RI (mammals)
	 Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  	 Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$Ea) &!is.na(Raw$RI))
	 X= Dat$Wa; Y=Dat$Ea*Dat$RI;
#	 Y_ax=seq.log(10^-3, 10^3, 10);	X_ax=seq.log(10^-5, 10^4, 10);
#	 Y_lab =c(expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3))	
#	 X_lab =c(expression(10^-5), expression(10^-4), Y_lab,expression(10^4))
		
	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="darkgray", pch = 19, cex=1.2, xaxs="i", yaxs="i", family="helvetica")
	  points(X, Y, type="p",pch = 1, cex=1.3, lwd = 0.25)
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext("Lifetime reproductive investment (kg)", side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
		mtext("Adult size (kg)", side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
	
		#fit slope
		fit<-slope.test(log10(Y), log10(X), test.value = 1.0, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=1.5, col ="darkgray") 
dev.off()

#--------------------------------------------------------------------------------
#3 row version 
postscript("fig3.eps", horizontal=FALSE, onefile=FALSE, height=10,width=6,pointsize=10, paper ="special" )
par(mfrow=c(3,1))
#1. Annual reproductive output
	 Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  	 Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
	 X= Dat$Wa; Y=Dat$RI;
	 Y_ax=seq.log(10^-6, 10^3, 10);	X_ax=seq.log(10^-5, 10^4, 10);
	 Y_lab =c(expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3))	
	 X_lab =c(Y_lab[2:10], expression(10^4))
		
	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="red", pch = 19, cex=1.0, xaxs="i", yaxs="i", family="helvetica")
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext("Annaul reproductive investment (kg)", side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
		#mtext("Adult size (kg) [log scale]", side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
	
		#fit slope
		fit<-slope.test(log10(Y), log10(X), test.value = 0.75, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2, col ="red") 

	#LOAD PLANT DATA
	Raw <- read.table("Plants.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
 	Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$RI))
      X= Dat$Wa; Y=Dat$RI;
	points(X, Y, type="p", pch=19, cex=1.0, col="blue")
	fit<-slope.test(log10(Y), log10(X), test.value = 0.75, method = 'SMA')
	Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
	curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2, col ="blue") 

#2. Reproductive lifespan (mammals)
	 Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  	 Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$Ea))
	 X= Dat$Wa; Y=Dat$Ea;
	 Y_ax=seq.log(10^-2, 10^2, 10);	X_ax=seq.log(10^-5, 10^4, 10);
	 Y_lab =c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2))	
	 X_lab =c(expression(10^-5), expression(10^-4), expression(10^-3), Y_lab, expression(10^2),expression(10^4))
		
	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="red", pch = 19, cex=1.0, xaxs="i", yaxs="i", family="helvetica")
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext("Average reproductive lifespan (kg)", side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
		#mtext("Adult size (kg) [log scale]", side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
	
		#fit slope
		fit<-slope.test(log10(Y), log10(X), test.value = 0.25, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2, col ="red") 
	
#3. Lifetime RI (mammals)
	 Raw <- read.table("Mammals.txt", header = T, sep="\t",strip.white = TRUE, blank.lines.skip = TRUE,)
  	 Dat <-subset(Raw, !is.na(Raw$Wa)&!is.na(Raw$Ea) &!is.na(Raw$RI))
	 X= Dat$Wa; Y=Dat$Ea*Dat$RI;
	 Y_ax=seq.log(10^-3, 10^3, 10);	X_ax=seq.log(10^-5, 10^4, 10);
	 Y_lab =c(expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2), expression(10^3))	
	 X_lab =c(expression(10^-5), expression(10^-4), Y_lab,expression(10^4))
		
	par(pty="s")
	plot(X, Y, type="p",log="xy", axes=F,ann=F, ylim=c(min(Y_ax), max(Y_ax)),xlim=c(min(X_ax), max(X_ax)),
		col="red", pch = 19, cex=1.0, xaxs="i", yaxs="i", family="helvetica")
		axis(2, at=Y_ax, labels=Y_lab, las=1, tck=0.012, , cex.axis=0.8)
		axis(4, at=Y_ax, labels=F,  tck=0.012)
		axis(1, at=X_ax, labels=X_lab, las=1, tck=0.012, cex.axis=0.8)
		axis(3, at=X_ax, labels=F, las=1, tck=0.012)
		mtext("Lifetime reproductive investment (kg)", side = 2, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
		mtext("Adult size (kg)", side = 1, line = 3.0, outer = F, at= NA, adj = 0.5, cex =1.1)	
	
		#fit slope
		fit<-slope.test(log10(Y), log10(X), test.value = 1.0, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2, col ="red") 
dev.off()
	
fit; Int; cor(log10(Y), log10(X))^2; length(X)
