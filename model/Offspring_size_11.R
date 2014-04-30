#Daniel Falster - offspring size model
#last modified Dec 2007

#----------------------------------------------------------
#SETUP: load libraries, set working directory
require(smatr)	#for SMA line-fitting
setwd("/Users/dfalster/Ecology/_Games/offspring size/model")  #setwd("F:\\_Games\\offspring size\\model")
rm(list=ls(all=TRUE)) 		#clear memory
source("ESS_functions_11.R")	#load sub routines

#----------------------------------------------------------
#load data
raw.mam <- read.table("../data/mammals.txt", header = T, sep="\t");
obs.mam <- subset(raw.mam, !is.na(raw.mam$Wa)&!is.na(raw.mam$W0));
raw.plt <- read.table("../data/Plants.txt", header = T, sep="\t")
obs.plt <- subset(raw.plt, !is.na(raw.plt$Wa)&!is.na(raw.plt$W0))

Wa.Range.mam =c(min(obs.mam$Wa), max(obs.mam$Wa));
Wa.Range.plt =c(min(obs.plt$Wa), max(obs.plt$Wa));
#----------------------------------------------------------
#Parameters
  q <- 0.2 	       #density independent mortality
  r <- 0.75	       #thinning mortality
  k <-10^-8
  alph_c <- 2; alph_k <-3  #for assymetric competition model
  RGR_switch = 1;  #1 = assume constant RGR
  g=0.8;		   #defines fraction of max size for maturity (only needed if RGR_switch !=1)

#Plants
	D=2; c_LRI = 0.895/35; b_LRI<- 0.6322944; s<- 0.00677 #optimal fit
#Mammals
    D=1; c_LRI = 0.895/20; b_LRI<- 0.9519; s=0.944; #optimal fit
#-------------------------------------------------------------------------------------
#Figures
#source("plots/Fitted_lines_version_1.R")
source("plots/Fitted_lines_version_2.R")
source("plots/Smith-Fretwell_curve_2.R")
source("plots/Competitive_asym.R")
source("plots/Fig_6-slope_change.R")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#MODEL ANALYSIS
#Check analytic solution against 2 numeric methods
plot(ESS_analytic(obs.mam$Wa, LRI(obs.mam$Wa)), ESS_numeric(obs.mam$Wa, LRI(obs.mam$Wa), 1e-5), log="xy", type='b');
curve(1*x, from = min(obs.mam$W0), to = max(obs.mam$W0), add = T, col='red')

plot(ESS_analytic(obs.mam$Wa, LRI(obs.mam$Wa)), ESS_numeric2(obs.mam$Wa, LRI(obs.mam$Wa), 1e-5, 600), log="xy", type='b');
curve(1*x, from = min(obs.mam$W0), to = max(obs.mam$W0), add = T, col='red')
#-------------------------------------------------------------------------------------
#Compare with constant and declining RGR=1
	D=1; RGR_switch =0;
	Wa = obs.mam$Wa;
	#compare predictions from two models
	plot(Wa, ESS_analytic(Wa, LRI(Wa)), type='l', log='xy');
	points(Wa, ESS_numeric(Wa, LRI(Wa), 1e-4), type='l', log='xy', col='red')

	#investigate change in slope as function B and competitive assymetry
	Wa = c(mean(obs.mam$Wa), max(obs.mam$Wa))
	plot(c(0,1),c(0,1), type='n', xlim= c(0.0, 1.0), ylim= c(-1, 0), xlab = 'Slope with constant RGR', ylab = 'Slope difference')
	for(alph_k in c(0.5,1,2,4, 8,16, 32, 64, 128, 256, 512, 1024))
		{
		X<-NULL; Y<-NULL; i=1;
		for(b_LRI in seq(1.0, 0.1, -0.01))
			{
			Y[i]= slope.test(log10(ESS_numeric(Wa, LRI(Wa), k)), log10(Wa), test.value = 1.0, method = 'SMA')$b;
			X[i] = (b_LRI-q)/(1-q);
			i=i+1;
			}
		points(X,Y-X, type='l')
		}
	points(X,0*X, type='l', lty= 'dashed')
 #-------------------------------------------------------------------------------------

#SOLVE ESS
#plot observed data
Wa=obs.mam$Wa; W0=obs.mam$W0;
	ll_plot(Wa, W0, "Adult size", "Offspring size")

#compare observed vs predicted offspring size
W0_pred= ESS_analytic(Wa, LRI(Wa));
	ll_plot(W0, W0_pred, "Observed", "Predicted")
	points(W0, W0, type='l');
	points(W0, W0_pred, type='p', col="blue");

#find optimum fit for slope
	W0_pred<-ESS_analytic(Wa, LRI(Wa));
	slope= slope.test(log10(W0_pred), log10(W0), test.value = 1.0, method = 'MA')$b
	count=1; mult=0.1;
	#choose to incraese or decrease slope
	flag=1; if(slope> 1){ flag=-1;}
	while(abs(slope-1) >10^-3 && count < 100)
		{
		temp=b_LRI;
		b_LRI = b_LRI *(1+flag*mult);
      	W0_pred<-ESS_analytic(Wa, LRI(Wa));
      	slope= slope.test(log10(W0_pred), log10(W0), test.value = 1.0, method = 'MA')$b
	    	print(c(count,temp, b_LRI,  slope, mult));
          	if((slope-1)*flag > 0) #increment too large  - overshoots optimum
			{mult = mult*0.5; b_LRI =temp;}
          	else  {count = 1+count;}
		}
	ll_plot(W0, W0_pred, "Observed", "Predicted")

#find optimum fit for intercept
	W0_pred<-ESS_analytic(Wa, LRI(Wa));
	slope= slope.test(log10(W0_pred), log10(W0), test.value = 1.0, method = 'MA')$b
	int<-10^(mean(log10(W0_pred))-slope*mean(log10(W0)))
	count=1; mult=0.5;
	#choose to incraese or decrease INTERCEPT
	flag=1; if(int > 1){ flag=-1;}
	while(abs(int-1) >10^-3 && count < 100)
		{
		temp=s;
		s = s *(1+flag*mult);
      	pred2<-ESS_analytic(Wa, LRI(Wa));
      	slope= slope.test(log10(pred2), log10(W0), test.value = 1.0, method = 'MA')$b
	    	int<-10^(mean(log10(pred2))-slope*mean(log10(W0)))
		print(c(count,temp, s,  int, mult));
          	if((int-1)*flag > 0) #increment too large  - overshoots optimum
			{mult = mult*0.5; s =temp;}
          	else  {W0_pred=pred2;}
		count = 1+count;
		}
	ll_plot(W0, W0_pred, "Observed", "Predicted")

#other properties at ESS
	W0_pred<-ESS_analytic(Wa, LRI(Wa));

	#Plot survival curves for ESS offsrping size species by species
	Obs<-obs.mam;
	pdf("survival.pdf", width=6, height=6, onefile=TRUE)
	for(i in 1:length(Obs$Wa))
		{
    	Plot.Surv(W0_pred[i], Wa[i], LRI(Wa[i]));
  		legend("bottomright", leg = c(
        		paste(as.character(Obs$Genus[i]), " ",as.character(Obs$species[i])),
				paste("Wa = ",format(Obs$Wa[i],digits=2), ", ", "Ea = ",format(Obs$Ea[i],digits=2)),
          		paste("W0 = ", format(Obs$W0[i],digits=2)," / ", format(W0_pred[i], digits =2)),
				paste("Wt = ", format(wt(W0_pred[i], Obs$Wa[i], LRI(Obs$Wa[i])),digits=2))
			));
     	 }
  	dev.off()

#------------------------------------------------------------------------------
#CHECKING THAT ESS IS INDEED ESS
W0_pred= ESS_analytic(Wa, LRI(Wa));
	#plot fitness function for given species
	i = 300;
	plot(function(x)R0(x, w0.Res = W0_pred[i], wa = Wa[i], LRI(Wa[i])), min(W0), max(W0), log="x", xlab = "R0(mut, res)", ylab = "R0")
	#make pairwise inavasibility plot
  	#??????

#------------------------------------------------------------------------------
#BASIC MODEL OUTPUT
#seed production v body size
	W0 = W0.init;
	plot(function(x)n(W0, LRI(x)),
		Wa.min, Wa.max, col = 1, lwd =1.5, log="xy",	#formatting
	      xlim = c(Wa.min*0.75, Wa.max*1.25),  xlab = "Adult size (kg) [log scale]", ylab= "Lifetime seed prod (life-1) [log scale]")
		curve(n(2*W0, LRI(x)), min(Obs$Wa), max(Obs$Wa),add = TRUE, col="blue")
#seed production v seed size
	WA = 10;
	plot(function(x)n(x, LRI(WA)),
		W0.min, W0.max, col = 1, lwd =1.5, log="xy",	#formatting
	      xlim = c(W0.min*0.75, W0.max*1.25),  xlab = "Offspring size (kg) [log scale]", ylab= "Lifetime seed prod (life-1) [log scale]")
		curve(n(x, LRI(WA/2)), min(Obs$Wa), max(Obs$Wa),add = TRUE, col="blue")
#Survival curves
	WA =10;
	Plot.Surv(0.01, WA, LRI(WA))

#Wt vs body size
	plot(function(x)wt(W0, x, LRI(x)),
		Wa.min, Wa.max, col = 1, lwd =1.5, log="xy",	#formatting
	      xlim = c(Wa.min*0.75, Wa.max*1.25),  xlab = "Adult size (kg) [log scale]", ylab= "Size at onset of self-thinning(kg)[log scale]")

#Wt vs W0
	#ll_plot(W0, wt(W0, Wa, lri), "Offspring size (kg)", "Size at onset thinning (kg)")
	ll_plot(W0, wt(W0, Wa, Lri_dat), "Offspring size (kg)", "Size at onset thinning (kg)")
	curve(1*x, min(Obs$W0), max(Obs$W0), add=TRUE)

#plot competition function
plot(function(x)alpha_ij(x, 1), 0.1, 10, log="x",  xlab="W1/W2", ylab = "a(ij) [log]")#, ylim = c(0,2))
	for(i in 1:8)
		{alph_k = i; curve(alpha_ij(x, 1),0.1,10, add=TRUE);}
