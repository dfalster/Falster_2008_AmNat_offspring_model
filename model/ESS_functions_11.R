#Dan Falster
#Functions for Optimal Offspring size model
#last modified Dec 4 2006

TOL = 10^-9
DERIV = 10^-6  #percentage difference when taking numerical derivatives

#returns vector from lo to hi with multiplication steps of incr
seq.log<-function(lo, hi, incr){temp<-NULL; while(lo < hi) {temp<-c(temp, lo); lo=lo*incr;}
					  if(max(temp)<hi) temp<-c(temp, hi); return (temp);} 

#offspring production - n estimates offspring production through scaling with body size, n_obs uses observed values for E and RI
	n <-function(w0, lri) {return(lri/w0);}    #Offspring size-number tradeoff
	RI<-function(wa)	{return(c_RI*wa^b_RI);}  #Rerproductive investment as function adult size
	RL<-function(wa)	{return(c_RL*wa^b_RL);}  #Rerproductive lifespan as function adult size
	LRI<-function(wa) {return(c_LRI*wa^b_LRI);}

#size at thinning onset as function offsprinf size, adult size, and lifetime reproductive investment LRI
	wt <- function(w0, wa, lri) { 
  		return((1/(s*lri)*wa^r*w0^(1-q)*D^(k/w0-q))^(1/(r-q)));} 
#surival functions
	s_est_phase <- function(w,w0){return (s*(w/w0)^(-k/w0));} 
	s_densindep_phase <- function(w,w0){return (s*D^(q-k/w0)*(w0/w)^(q));}
	s_ST_bound<- function(w, w0, wa, lri){ return ((w/wa)^(-r)/n(w0,lri));}
	#survival to adulthood - for resident only
	s_adult <-function(w0, wa,wt) { return(s*D^(q - k/w0)*wt^(r-q)*w0^q*wa^(-r)); }                        
	#survival to size w of mutant in env of residents 
	s_mut <-function(w, w0.mut, w0.res, wa, lri){
  		#standradise array dimensions for delaing with plotting etc
		LEN = max(length(w), length(w0.mut),length(w0.res), length(wa),length(lri));
		#print(LEN);
		if(length(w) ==LEN) {W = w;} else {W= rep(w[1], LEN);}
		if(length(w0.mut)==LEN) {W0.MUT = w0.mut;} else {W0.MUT= rep(w0.mut[1], LEN);}
		if(length(w0.res)==LEN) {W0.RES = w0.res;} else {W0.RES= rep(w0.res[1], LEN);}
		if(length(wa)==LEN) {WA = wa;} else {WA= rep(wa[1], LEN);}
		if(length(lri)==LEN) {LRI2= lri;} else {LRI2= rep(lri[1], LEN);}
		
		#solves for point of thinning onset
		WT <- wt(W0.RES, WA, LRI2);
		gam = gamma(W0.MUT, W0.RES, WA, WT);
		WT_mut = gam*WT*W0.MUT/W0.RES;
	  	surv = rep(0,times=LEN);
        for(i in 1:LEN)
           	{
           	#check to ensure offspring not larger than adult!
			if(W[i]>WA[i]) 
				{surv[i]=0;}  
			#establishment phase
			else if (W[i]<= D*W0.MUT[i])   
				{surv[i] = s_est_phase(W[i],W0.MUT[i]); } 
           	else 
           		{
		    	#excess juveniles - crash at onset of thinning 
		    	if (WT_mut[i]< D*W0.MUT[i]) 
				     {s_crash = s_ST_bound(D*W0.RES[i], W0.RES[i], WA[i], LRI2[i])/s_est_phase(D*W0.RES[i],W0.RES[i])   #magnitude of crash
               		 surv[i] = s*D^(- k/W0.MUT[i]) * (s_crash/alpha_ij(W0.MUT[i],W0.RES[i]))*(W[i]/D/W0.MUT[i])^(-r*alpha_ij(W0.MUT[i],W0.RES[i]));}
	            #too few juveniles - thinning starts after Wa --> density indept mortality only
 	       		else if (WT_mut[i]> WA[i]) 
					{surv[i] = s_densindep_phase(W[i],W0.MUT[i]);} 
           	   #normal case - thinning starts between D*w0 & Wa
           	   else 
			  		{if(W[i]<= WT_mut[i]) 
				  		{surv[i]=s_densindep_phase(W[i],W0.MUT[i]);}
			   		else                
				  		{
				  		surv[i] = s*D^(- k/W0.MUT[i]) * (gam[i]*WT[i]/D/W0.RES[i])^(-q)*(W[i]*W0.RES[i]/WT[i]/W0.MUT[i]/gam[i])^(-r*alpha_ij(gam[i]*W0.MUT[i],W0.RES[i]));}
				  	}
		   		}
           	}
	    return(surv); 
		}                                           

#adjustment function for RGR
gamma<-function(w0.mut, w0.res, wa, wt.res){
	gam= (( (wa/wt.res)^0.25*(1-(w0.res/w0.mut)^0.25) + (wa/w0.mut)^0.25 - g^0.25) / ( (wa/w0.res)^0.25 - g^0.25) )^4;
	if(RGR_switch == 1) 
		{gam=0*gam + 1;}
	return(gam);
	}

#for calculation of ESS
	#R0 function - fitness
	R0 <- function(w0.mut, w0.res, wa,lri) {
		return(n(w0.mut,lri)*s_mut(wa, w0.mut, w0.res, wa, lri));}
		
	#selection gradient
	dR0 <-function(w0, wa, lri) {
		return((R0(w0*(1+DERIV), w0, wa, lri) - R0(w0*(1-DERIV), w0, wa, lri))/(w0*2*DERIV));} 
	
    #solve for ESS with analyictal solution and adult size scaling
    ESS_analytic <-function(wa, lri){
	 	return(exp(2*(r-q)*(1-r)/r/alph_k /(q-1))*(s*lri/wa^q)^(1/(1-q)));}
	
	#solve for ESS by locating zero in selection gradient
	ESS_numeric<-function(wa, lri, min) {
		Sol = 0*wa;
		for(i in 1:length(wa))	#do for all values of wa if a vector 
			{Sol[i]= uniroot(dR0, wa[i]*c(min, 1), tol= TOL, wa = wa[i], lri = lri[i])$root;}
		return(Sol);}		
	
	#routine for solving ESS through impelementing canon equ
	ESS_numeric2 <-function(wa, lri, delta, max_counts)
		{res<- delta*wa
		for(i in 1:length(wa))	#do for all values of wa if a vector 
			{count <-0
			 mult = 0.5
	 		 der<-dR0(res[i],wa[i], lri[i])
	 		 #check to see if optimum is greater or lower seed size
	 		 if(der >0){ flag=1;}#increase seed mass
			 else  {flag=-1;}#decrease seed mass
			 #increment until opitmum located
			 while(abs(der) > TOL && count < max_counts)
			    {
			    W0_next = res[i]*(1+flag*mult);
          		    der<-dR0(W0_next,wa[i], lri[i]); 
 	        	    #print(  c(count, res, dR0(res[i], wa[i], lri[i]), mult)) 
          		    if(der*flag < 0  || W0_next > wa[i]) #increment too large  - overshoots optimum
				      {mult = mult*0.5}   #decrease increment
          		    else
              		{res[i] = W0_next} 
          		    count = 1+count;
         		    }
			print(count);
      		}
		return(res);
		}
		
#competition effect
	alpha_ij <-function (w1,w2){return (alph_c/(1+(w1/w2)^alph_k));}
	
#plot survival curves
#Create a plot of A-Ci curves showing Aj and Av
	Plot.Surv=function(w0, wa, lri){
 		WT <- wt(w0, wa, lri); 
		 plot(function(x)s_est_phase(x, w0=w0), w0, D*w0, lwd = 2, col = "red",xlab = "W (kg) [log]", xlim=c(w0/2,wa*2), ylim = c(0.0001,s*2), ylab= "Survival (0-1)", log="xy")
		 curve(s_densindep_phase(x, w0=w0), from = D*w0, to = max(min(WT, wa), D*w0), lwd = 2, col = "blue", add=TRUE)
	       curve(s_ST_bound(x, w0=w0, wa=wa, lri=lri), from = min(max(WT, D*w0),wa), to = wa, lwd = 2,col = "green", add=TRUE)
	       curve(s_mut(x, w0.mut=w0, w0.res =w0, wa =wa, lri =lri), from = w0, to = wa, lwd = 1,col = "black", add=TRUE)
 		 #points(c(W0, WA), c(s_est_phase(w0,w0), s_combined(wa, w0,w0,wa)), type = "p", cex = 1.5, )  
  		 #plot(c(W0, WT, WA), c(s_est_phase(W0,W0), s_thinning_phase(WT, W0,W0,WA), s_thinning_phase(WA, W0,W0,WA)), type = "p", cex = 1.5, log="xy")
  		legend("bottomleft", leg = c("establishment", "dens indep", "thinning", "all"), lwd=c(2,2,2,1), col=c("red", "blue", "green", "black"))
  		}


#log-log plot with legend and 1:1
	ll_plot=function(X,Y, XLAB, YLAB){
		plot(X, Y, log="xy", xlab=XLAB, ylab=YLAB)
		fit<-slope.test(log10(Y), log10(X), test.value = 1.0, method = 'SMA')
		Int<-10^(mean(log10(Y))-fit$b*mean(log10(X)))
		curve(Int*x^fit$b, min(X), max(X), add=TRUE, lwd=2) 
		r2<- cor(log10(Y), log10(X))^2; n = length(X)
  		leg<-c(paste("B = ",format(fit$b, digits =4)), paste("Int = ",format(Int, digits =4)),paste("n = ",n), paste("R2 = ",format(r2, digits =3)))
  		legend(x="topleft",legend=leg, lty=0, cex=1.0)
		}

#PIP 
require(lattice)  #for 3D plotting PIP 
PIP = function(low, up, wa){
	v= seq.log(low, up, (up/low)^0.05);
	N=length(v);
	x = matrix(rep(v, time=N), ncol= N, nrow=N)
	y = t(x);
	r= 0*x;
	for(i in 1:N){
	    for(j in 1:N){
		r[i,j] = max(0,R0(x[i,j], y[i,j], wa, LRI(wa))-1)
			}
		}
	wireframe(r ~ log10(x) * log10(y), 
          scales = list(arrows = FALSE),
	        perspective = F,
          drape = T, colorkey = TRUE,
          screen = list(z = 90, x =0, y=0))
	}
	
