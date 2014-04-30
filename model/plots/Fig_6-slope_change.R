#MS FIG 6

postscript("plots/fig6.eps", horizontal=FALSE, onefile=FALSE, 
height=6,width=6,pointsize=10, paper ="special" )

#investigate change in slope as function B and competitive assymetry 
   D=1; c_LRI = 0.895/20; b_LRI<- 0.9519; s=0.944; #optimal fit

	Wa = c(mean(obs.mam$Wa), max(obs.mam$Wa))
	
plot(c(0,1),c(0,1), type='n', xlim= c(0.0, 1.0), ylim= c(0, 1.0), lwd=1.0, axes=F,ann=F,xaxs="i", yaxs="i");
  box()
  axis(1, at=seq(0, 1, 0.2), labels=seq(0,1,0.2), las=1, tck=0.012)
  axis(1, at=seq(0, 1, 0.1), labels=F, tck=0.005)
  axis(3, at=seq(0, 1, 0.2), labels=F, tck=0.012); 
  axis(3, at=seq(0, 1, 0.1), labels=F, tck=0.005)
  axis(2, at=seq(0, 1, 0.2), labels=seq(0,1,0.2), las=1, tck=0.012)
  axis(2, at=seq(0, 1, 0.1), labels=F, tck=0.005)
  axis(4, at=seq(0, 1, 0.2), labels=F, tck=0.012); 
  axis(4, at=seq(0, 1, 0.1), labels=F, tck=0.005)
  mtext('Slope with constant RGR', side =1, line =3);
  mtext('Slope with declining RGR', side =2, line =3);
  
	
	for(alph_k in c(64, 16, 4)) 
		{	
		X<-NULL; Y<-NULL; i=1;
		for(b_LRI in seq(1.0, 0.2, -0.1))
			{
			Y[i]= slope.test(log10(ESS_numeric(Wa, LRI(Wa), k)), log10(Wa), test.value = 1.0, method = 'SMA')$b;
			X[i] = (b_LRI-q)/(1-q);
			i=i+1;
			}
		points(X,Y, type='l')
		}
	points(c(0,1),c(0,1), type='l', lty= 'dashed')
	
	
	for(alph_k in c(1)) 
		{	
		X<-NULL; Y<-NULL; i=1;
		for(b_LRI in seq(1.0, 0.4, -0.1))
			{
			Y[i]= slope.test(log10(ESS_numeric(Wa, LRI(Wa), k)), log10(Wa), test.value = 1.0, method = 'SMA')$b;
			X[i] = (b_LRI-q)/(1-q);
			i=i+1;
			}
		points(X,Y, type='l')
		}
	points(c(0,1),c(0,1), type='l', lty= 'dashed')

	for(alph_k in c(0.25)) 
		{	
		X<-NULL; Y<-NULL; i=1;
		for(b_LRI in seq(1.0, 0.6, -0.1))
			{
			Y[i]= slope.test(log10(ESS_numeric(Wa, LRI(Wa), k)), log10(Wa), test.value = 1.0, method = 'SMA')$b;
			X[i] = (b_LRI-q)/(1-q);
			i=i+1;
			}
		points(X,Y, type='l')
		}
	points(c(0,1),c(0,1), type='l', lty= 'dashed')

	
	dev.off()
